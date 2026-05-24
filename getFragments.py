import os
import sys
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
def write_xyz(filename, symbols, coords, comment=""):
    n = len(symbols)
    with open(filename, "w") as f:
        f.write(f"{n}\n")
        f.write(comment + "\n")
        for sym, (x, y, z) in zip(symbols, coords):
            f.write(f"{sym:2s}  {x:15.8f}  {y:15.8f}  {z:15.8f}\n")
filename = "test.xyz"
input_folder = os.getcwd()
out_folder = input_folder
file_path = os.path.join(input_folder, filename)
print(f"Elaborating: {file_path}")
matrice = read_xyz(file_path)[0]
elements = read_xyz(file_path)[1]
if elements[1] != "N":
    dic = np.array([1, 22, 25, 23, 26, 28, 27, 24, 13, 14, 15, 16, 17, 12, 8, 2, 7, 6, 5, 4, 3, 18, 19, 20, 21, 11, 10, 9])
    new_matrice = np.zeros((len(matrice), 3))
    new_elements = np.copy(elements)
    for i in range(len(matrice)):
        new_matrice[dic[i] - 1] = matrice[i]
        new_elements[dic[i] - 1] = elements[i]
    write_xyz(os.path.join(out_folder, filename.strip(".xyz") + "_new.xyz"), new_elements, new_matrice)
    write_xyz(os.path.join(out_folder, filename.strip(".xyz") + "_new_align.xyz"), new_elements,align(new_matrice, new_elements))
    new_matrice = (align(new_matrice, new_elements))
else:
    write_xyz(os.path.join(out_folder, filename.strip(".xyz") + "_align.xyz"), elements,align(matrice, elements))
    new_matrice = (align(matrice, elements))
    new_elements = np.copy(elements)

dic_frag = [[26,27],[21,22],[23,24],[25],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]]
dic_label = ["COax","COeq1","COeq2","Cl","Byp"]
n = new_matrice.shape[0]
all_idx = np.arange(n)
name, ext = os.path.splitext(filename)
for label, frag in zip(dic_label, dic_frag):
    mask = np.zeros(n, dtype=bool)
    mask[frag] = True
    frag_matrix = new_matrice[mask]
    rest_matrix = new_matrice[~mask]
    write_xyz(os.path.join(out_folder,name +"_"+ label + ".xyz"), new_elements[mask], frag_matrix)
    write_xyz(os.path.join(out_folder, name +"_"+ label + "_rest.xyz"), new_elements[~mask], rest_matrix)

