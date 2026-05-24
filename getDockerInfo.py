import os
import re
import numpy as np

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

def write_xyz(filename, symbols, coords, comment=""):
    n = len(symbols)
    with open(filename, "w") as f:
        f.write(f"{n}\n")
        f.write(comment + "\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s}  {x:15.8f}  {y:15.8f}  {z:15.8f}\n")

def read_xyz_internal(data_lines):
    coords = []
    elements = []
    for line in data_lines:
        parts = line.split()
        if len(parts) == 4:
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
            elements.append(parts[0])
    return np.array(coords), np.array(elements)

input_folder = os.path.join(os.getcwd(),"Input", "Docker")
filename = "Fac_5_DCM_XTB.docker.struc1.all.optimized.xyz"
energy_pattern = re.compile(r"Einter=\s*(-?\d+\.\d+)")
energies = []
with open(os.path.join(input_folder,filename), 'r') as f:
    lines = f.readlines()
    n_atoms = int(lines[0])
    n_structures = len(lines)//(n_atoms+2)
    for i in range(n_structures):
        match = energy_pattern.search(lines[i*(n_atoms+2)+1])
        if match:
            energies.append(float(match.group(1)))
        coords,elements = read_xyz_internal(lines[i*(n_atoms+2)+2:(i+1)*(n_atoms+2)])
        new_coords = align(coords,elements)
        filename_info = filename.split(".")
        write_xyz(os.path.join(input_folder,f"{filename_info[0]}_{filename_info[2]}_{i+1}.xyz"),elements,new_coords,comment=f"Interaction Energy: {energies[i]}")
with open(os.path.join(input_folder, f"{filename_info[0]}_{filename_info[2]}_energies.csv"), "w") as f2:
    f2.write("Structure,Energy\n")
    for i in range(len(energies)):
        f2.write(f"{filename_info[0]}_{filename_info[2]}_{i+1},{energies[i]}\n")