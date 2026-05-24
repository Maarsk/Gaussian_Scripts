import os
import numpy as np
import re
import tempfile

def write_xyz(filename, symbols, coords, comment=""):
    n = len(symbols)
    with open(filename, "w") as f:
        f.write(f"{n}\n")
        f.write(comment + "\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s}  {x:15.8f}  {y:15.8f}  {z:15.8f}\n")

def align(coord,elements,invert=False):
    atomic_symbols = {
    1: "H", 2: "He",
    3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
    19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
    27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
    35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce",
    59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
    67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta", 74: "W",
    75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb",
    83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
    99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg",
    107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl",
    115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
    }
    symbol_to_Z = {v: k for k, v in atomic_symbols.items()}
    z_element = np.zeros(len(elements))
    for i,atom in enumerate(elements):
        z_element[i] = symbol_to_Z[atom]
    if invert == True:
        I = np.array(([-1,0,0],[0,-1,0],[0,0,-1]))
        coord = coord @ I
    byp = coord[0:21]
    c = byp.mean(axis=0)
    q = byp - c
    C = q.T @ q
    vals, vecs = np.linalg.eigh(C)
    n = vecs[:, np.argmin(vals)]
    n = n / np.linalg.norm(n)
    for i,atom in enumerate(coord):
        I_n_minus = sum(z_element[i]*(atom-n)**2)
        I_n_plus = sum(z_element[i]*(atom+n)**2)
    if I_n_plus > I_n_minus:
        n = n
    print(I_n_plus, I_n_minus)
    r_ref = coord[0]
    x0 = r_ref - c
    x = x0 - np.dot(x0, n) * n
    x_hat = x / np.linalg.norm(x)

    y_hat = np.cross(n, x_hat)
    y_hat = y_hat / np.linalg.norm(y_hat)

    R = np.vstack([x_hat, y_hat, n])  # shape (3, 3)

    r_shift = coord - c
    coords_new = (R @ r_shift.T).T  # (N, 3)
    coords_new = coords_new - coords_new[0]
    if elements[1] == "N" and coords_new[1,1] < 0:
        print("N")
        R = np.vstack([x_hat, -y_hat, -n])
        r_shift = coord - c
        coords_new = (R @ r_shift.T).T  # (N, 3)
        coords_new = coords_new - coords_new[0]
    return coords_new

def read_xyz(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data_lines = lines[2:]
    coords = []
    elements = []
    for line in data_lines:
        parts = line.split()
        if len(parts) == 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
            elements.append(parts[0])
    return np.array(coords), np.array(elements)

def extract_gaussian_energy(file_path, variable):
    last_energy = []
    if variable == 0:
        energy_pattern = re.compile(
            r"SCF Done:\s+E\([A-Za-z0-9\-]+\)\s+=\s+(-?\d+\.\d+)"
        )
    elif variable == 1:
        energy_pattern = re.compile(r"Total Energy.*?=\s*(-?\d+\.\d+)")
    with open(file_path, 'r', errors='ignore') as f:
        for line in f:
            match = energy_pattern.search(line)
            if match:
                last_energy.append(float(match.group(1)))

    if last_energy is None:
        raise ValueError("No SCF energy found in file.")
    if last_energy == []:
        print("No SCF energy found in file.")
    return last_energy

def energy_pattern(file_path):
    pattern = re.compile(r"Total Energy.*?=\s*(-?\d+\.\d+)")
    with open(file_path, 'r', errors='ignore') as f:
        for line in f:
            if pattern.search(line):
                return 1
    return 0

def extract_opt_geometry(filename):
    atomic_symbols = {
        1: "H", 2: "He",
        3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
        11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar",
        19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe",
        27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
        35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
        43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
        51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce",
        59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy",
        67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta", 74: "W",
        75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb",
        83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
        91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
        99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf", 105: "Db", 106: "Sg",
        107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds", 111: "Rg", 112: "Cn", 113: "Nh", 114: "Fl",
        115: "Mc", 116: "Lv", 117: "Ts", 118: "Og"
    }

    pattern = re.compile(r"Standard orientation:")
    opt_pattern = re.compile(r"Optimization completed.")
    energies = extract_gaussian_energy(filename, energy_pattern(filename))

    geometries = []
    opt_geometries = []
    opt_energies = []
    opt = False
    with open(filename, 'r') as f:
        lines = f.readlines()
    j = 0
    i = 0
    while i < len(lines):
        if opt_pattern.search(lines[i]):
            try:
                opt_energies.append(energies[j-1])
            except Exception as error:
                print("filename")
                print(error)
                break
            opt = True
        if pattern.search(lines[i]):
            i += 5
            block = []
            while i < len(lines) and not re.match(r'\s*-{5,}\s*', lines[i]):
                try:
                    parts = lines[i].split()
                    if len(parts) == 6:
                        center, atomic_num, atomic_type, x, y, z = parts
                        block.append({
                            "center": int(center),
                            "atomic_symbol": atomic_symbols[int(atomic_num)],
                            "atomic_type": int(atomic_type),
                            "x": float(x),
                            "y": float(y),
                            "z": float(z)
                        })
                    i += 1
                except ValueError:
                    print(f"Error parsing line {i}: {lines[i]}")
                    break
            j += 1
            geometries.append(block)  # ← correct place
            if opt == True:
                opt_geometries.append(block)
                opt = False
            i += 1
        else:
            i += 1
    return opt_geometries, opt_energies

input_folder = os.path.join(os.getcwd(), 'Input', "out_to_info", "Renamed", "Mer")
import os
import tempfile

for filename in os.listdir(input_folder):
    if filename.endswith(".out"):
        print(f"Processing {filename}...")

        # Extract geometry + energies
        opt_geometries, opt_energies = extract_opt_geometry(
            os.path.join(input_folder, filename)
        )

        # Create temporary XYZ
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".xyz", delete=False) as tmp:
            temp_xyz_path = tmp.name

            i = 0
            for block in opt_geometries:
                tmp.write(str(len(block)) + "\n")

                try:
                    tmp.write(f"Energy: {opt_energies[i]}\n")
                except IndexError:
                    tmp.write(f"Energy: {opt_energies[i - 1]}\n")

                i += 1

                for atom in block:
                    tmp.write(
                        f"{atom['atomic_symbol']:<3}"
                        f"{atom['x']:>12.6f}"
                        f"{atom['y']:>12.6f}"
                        f"{atom['z']:>12.6f}\n"
                    )

            tmp.flush()
            os.fsync(tmp.fileno())  # ensure OS flush

        # Now the file is closed → safe to read
        coords, elements = read_xyz(temp_xyz_path)
        isomer = filename.split("_")[0]
        print(isomer)
        #conver to new geometry
        if isomer == "Fac":
            if elements[4] == "Cl":
                dic = np.array([1, 22, 25, 23, 26, 28, 27, 24, 13, 14, 15, 16, 17, 12, 8, 2, 7, 6, 5, 4, 3, 18, 19, 20, 21, 11, 10, 9])
                print("old geometry 0")
            elif elements[7] == "Cl":
                dic = np.array([1, 27, 28, 22, 23, 24, 25, 26, 13, 14, 15, 16, 17, 12, 8, 2, 7, 6, 5, 4, 3, 18, 19, 20, 21, 11, 10, 9])
                print("old geometry 1")
            elif elements[25] == "Cl":
                dic = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
        elif isomer == "Mer":
            if elements[4] == "Cl":
                dic = np.array([1, 22, 25, 28, 26, 23, 27, 24, 13, 14, 15, 16, 17, 12, 8, 2, 7, 6, 5, 4, 3, 18, 19, 20, 21, 11, 10, 9])
                print("old geometry 0")
            elif elements[7] == "Cl":
                dic = np.array([1, 27, 28, 22, 23, 24, 25, 26, 13, 14, 15, 16, 17, 12, 8, 2, 7, 6, 5, 4, 3, 18, 19, 20, 21, 11, 10, 9])
                print("old geometry 1")
            elif elements[25] == "Cl":
                dic = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])


        new_coords = np.zeros((len(coords), 3))
        new_elements = np.copy(elements)
        for i in range(len(coords)):
            new_coords[dic[i] - 1] = coords[i]
            new_elements[dic[i] - 1] = elements[i]
        write_xyz(os.path.join(os.getcwd(), "Input","out_to_info","Aligned", isomer,filename.strip(".out") + "_new_align.xyz"), new_elements,
                  align(new_coords, new_elements),comment=f"Energy: {opt_energies[0]}")
        #add align


