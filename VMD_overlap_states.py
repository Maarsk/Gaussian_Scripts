import os
import subprocess
import tempfile
from PIL import Image
import re
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import csv
import math
from math import acos
from itertools import combinations

start = time.perf_counter()


VMD = r"C:\Program Files\University of Illinois\VMD2\vmd.exe"
def available_stuff(stuff_id,stuff_list):
    final_stuff = []
    for i in stuff_id:
        final_stuff.append(stuff_list[i])
    return final_stuff

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

def legend(col_id_mol, labels, output):
    vmd_colors = {
    0:  "#0000ff",  # blue
    1:  "#ff0000",  # red
    2:  "#808080",  # gray
    3:  "#ffa500",  # orange
    4:  "#ffff00",
    5:  "#d2b48c",
    6:  "#c0c0c0",
    7:  "#008000",
    8:  "#ffffff",
    9:  "#ffc0cb",
    10: "#00ffff",
    11: "#800080",
    12: "#00ff00",
    13: "#e0b0ff",
    14: "#cc7722",
    15: "#99ccff",
    16: "#000000",
    17: "#ffff66",
    18: "#cccc00",
    19: "#66ff66",
    20: "#009900",
    21: "#66ffff",
    22: "#009999",
    23: "#6666ff",
    24: "#000099",
    25: "#8f00ff",
    26: "#b266ff",
    27: "#ff00ff",
    28: "#cc00cc",
    29: "#ff6666",
    30: "#990000",
    31: "#ff9933",
    32: "#cc6600"
}

# Build legend handles
    handles = [
    Patch(facecolor=vmd_colors[cid], edgecolor="black", label=f"{labels[i]}")
    for i, cid in enumerate(col_id_mol)
    ]

    # Create figure with only legend
    fig, ax = plt.subplots(figsize=(3, len(col_id_mol)*0.5))
    ax.axis("off")

    legend = ax.legend(
        handles=handles,
        loc="center left",
        frameon=False,
        fontsize=14,
        handlelength=1.5,
        borderpad=1,
    )

    # Save legend as PNG/JPG
    plt.savefig(output, dpi=300, bbox_inches="tight", pad_inches=0.2)
    plt.close()
    print(f"Legend saved to {output}")

def indices_for_group(files, key="functional"):
    out = []
    for f in files:
        info = parse_filename(f)
        if key == "functional":
            out.append(functional_to_index[info["functional"]])
        elif key == "solvent":
            out.append(solvent_to_index[info["solvent"]])
        elif key == "state":
            out.append(state_to_index[info["state"]])
        elif key == "isomer":
            out.append(isomer_to_index[info["isomer"]])
    return out

def run_pipeline(files, output_png, col_id_mol, camera_path = "VMD2/camera.vmd", molecule_path = "VMD2/molecule_CPK_colid.vmd"):
    output_png = os.path.abspath(output_png).replace("\\", "/")
    camera_vmd = os.path.abspath(camera_path).replace("\\", "/")
    molecules_vmd = os.path.abspath(molecule_path).replace("\\", "/")
    cube_files = files

    with tempfile.NamedTemporaryFile(delete=False, suffix=".vmd") as tmp:
        script_path = tmp.name.replace("\\", "/")

        # FIRST CUBE → molecule + isosurface
        for i in range(len(cube_files)):
            first = cube_files[i].replace("\\", "/")
            tmp.write(f'''
    color change rgb 0 0 0.5 1        
    set xyz_path "{first}"
    set col_id {col_id_mol[i]}
    source "{molecules_vmd}"
    '''.encode())
            # APPLY CAMERA

        tmp.write(f'''
    source "{camera_vmd}"
    '''.encode())

        # FINAL RENDER
        tmp.write(f'''
set outfile "{output_png}.tga"
display resize 3000 3000
axes location off
display depthcue off
display projection Orthographic
render TachyonInternal $outfile
quit
'''.encode())

    subprocess.run([VMD, "-dispdev", "text", "-e", script_path])

# --- DEFINITIONS -------------------------------------------------------------
FUNCTIONALS = ["B3LYP", "PBE0", "TPSSH", "M062X", "CAMB3LYP", "wB97XD"]
SOLVENTS    = ["Vacuum", "DCM", "THF"]
STATES      = ["S0", "TDFT", "TD", "TDA", "S1"]
ISOMERS     = ["Fac", "Mer", "Cis", "Trans"]

functional_to_index = {f: i for i, f in enumerate(FUNCTIONALS)}
solvent_to_index    = {s: i for i, s in enumerate(SOLVENTS)}
state_to_index      = {s: i for i, s in enumerate(STATES)}
isomer_to_index     = {i: j for j, i in enumerate(ISOMERS)}

# --- ALIASES -----------------------------------------------------------------
# Add any typos, alternative spellings, or synonyms here
ALIASES = {
    # Functionals
    "PBE1": "PBE0",
    "PBE": "PBE0",
    "CAM-B3LYP": "CAMB3LYP",
    "WB97XD": "wB97XD",
    "WB97X-D": "wB97XD",
    "TPPSH" : "TPSSh",

    # Solvents
    "GP": "Vacuum",
    "GasPhase": "Vacuum",
    "SCM-DCM": "DCM",
    "PCM-DCM": "DCM",

    # States
    "T-DFT": "TDFT",
    "T-TDA": "TDA",
    "S1toS0": "S1",
    # Isomers
    "fac": "Fac",
    "mer": "Mer",
    "cis": "Cis",
    "trans": "Trans",
}


# --- NORMALIZATION -----------------------------------------------------------

def normalize_token(token):
    """Normalize any token using ALIASES, otherwise return unchanged."""
    if token in ALIASES:
        return ALIASES[token]
    return token


# --- PARSER ------------------------------------------------------------------

def parse_filename(fname):
    """
    Extract solvent, state, isomer, functional from a filename,
    regardless of order, with alias correction.
    """
    base = os.path.splitext(os.path.basename(fname))[0]
    parts = base.split("_")

    solvent = None
    state = None
    isomer = None
    functional = None

    for p in parts:
        p_norm = normalize_token(p)

        if p_norm in FUNCTIONALS:
            functional = p_norm
        elif p_norm in SOLVENTS:
            solvent = p_norm
        elif p_norm in STATES:
            state = p_norm
        elif p_norm in ISOMERS:
            isomer = p_norm

    return {
        "filename": fname,
        "solvent": solvent,
        "state": state,
        "isomer": isomer,
        "functional": functional
    }


# --- GROUPING ---------------------------------------------------------------

def group_files(files):
    """
    Build the three grouping dictionaries:
    1. same functional + same solvent (different states)
    2. same solvent + same state (different functionals)
    3. same functional + same state (different solvents)
    """

    same_func_same_solv = {}
    same_solv_same_state = {}
    same_func_same_state = {}

    for f in files:
        info = parse_filename(f)

        key1 = (info["functional"], info["solvent"])
        key2 = (info["solvent"], info["state"])
        key3 = (info["functional"], info["state"])

        same_func_same_solv.setdefault(key1, []).append(f)
        same_solv_same_state.setdefault(key2, []).append(f)
        same_func_same_state.setdefault(key3, []).append(f)

    return same_func_same_solv, same_solv_same_state, same_func_same_state

def geometry_information(filename):
    def angle(a, b, c):
        r1 = ((matrice[a, 0] - matrice[b, 0]) ** 2 + (matrice[a, 1] - matrice[b, 1]) ** 2 + (
                    matrice[a, 2] - matrice[b, 2]) ** 2) ** 0.5
        r2 = ((matrice[b, 0] - matrice[c, 0]) ** 2 + (matrice[b, 1] - matrice[c, 1]) ** 2 + (
                    matrice[b, 2] - matrice[c, 2]) ** 2) ** 0.5
        r3 = ((matrice[a, 0] - matrice[c, 0]) ** 2 + (matrice[a, 1] - matrice[c, 1]) ** 2 + (
                    matrice[a, 2] - matrice[c, 2]) ** 2) ** 0.5
        theta = acos((r1 ** 2 + r2 ** 2 - r3 ** 2) / (2 * r1 * r2))
        return math.degrees(theta)
    matrice = read_xyz(filename)[0]
    dist_int = np.zeros(6)
    dic = np.array([25, 26, 1, 11, 23, 21])
    # New_geom
    dist = np.zeros(len(matrice))
    for i in range(len(matrice)):
        for j in range(3):
            dist[i] = dist[i] + (matrice[0, j] - matrice[i, j]) ** 2
        dist[i] = dist[i] ** 0.5
    comb = list(combinations(dic, 2))
    matrix_angles = np.zeros(len(comb))
    for i in range(len(comb)):
        matrix_angles[i] = angle(comb[i][0], 0, comb[i][1])
    for i in range(6):
        dist_int[i] = dist[dic[i]]
    return np.concatenate([dist_int, matrix_angles])
    name, ext = os.path.splitext(filename)
    out_file = out_folder_mkdir(name, 'Geom_Values_out', 'csv')
    with open(out_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header_angles)
        writer.writerow(np.round(matrix_angles, 1))
        writer.writerow(header_dist)
        writer.writerow(np.round(dist_int, 4))


if __name__ == "__main__":
    working_folder = "VMD2/Output/Mer/Renamed"
    files = []

    for filename in os.listdir(working_folder):
        if filename.endswith(".xyz"):
            fullpath = os.path.join(working_folder, filename)
            files.append(fullpath)
    g1, g2, g3 = group_files(files)
    header = np.array(
        ('Filename','d(Re-Cl/A', 'd(Re-Cax/A', 'd(Re-N1/A', 'd(Re-N2/A', 'd(Re-C1/A', 'd(Re-C2/A', 'theta(Cl-Re-Cax/deg)', 'theta(Cl-Re-N1/deg)', 'theta(Cl-Re-N2/deg)', 'theta(Cl-Re-C1/deg)',
         'theta(Cl-Re-C2/deg)', 'theta(Cax-Re-N1/deg)', 'theta(Cax-Re-N2/deg)', 'theta(Cax-Re-C1/deg)',
         'theta(Cax-Re-C2/deg)', 'theta(N1-Re-N2/deg)', 'theta(N1-Re-C1/deg)', 'theta(N1-Re-C2/deg)',
         'theta(N2-Re-C1/deg)', 'theta(N2-Re-C2/deg)', 'theta(C1-Re-C2/deg)'))
    for k, v in g1.items():
        try:
            geom_matrix = []
            files = []
            for file in v:
                files.append(((os.path.basename(file)).split(".")[0]).strip("_new_align"))
                vec = geometry_information(file)  # 1×21 vector
                geom_matrix.append(vec)
            print(files)

            geom_matrix = np.array(geom_matrix)  # → shape (n_files, 21)
            with open(os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_geominfo.csv"), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                distances = np.round(geom_matrix[:, :6], 4)
                angles = np.round(geom_matrix[:, 6:], 1)
                geom_matrix_rounded = np.hstack([distances, angles])

                # Combine filenames + geometry rows
                rows = [
                    [fname] + list(row)
                    for fname, row in zip(files, geom_matrix_rounded)
                ]

                writer.writerows(rows)
        except Exception as e:
            print("ERROR", v)
            print(e)
    for k, v in g2.items():
        try:
            geom_matrix = []
            files = []
            for file in v:
                files.append(((os.path.basename(file)).split(".")[0]).strip("_new_align"))
                vec = geometry_information(file)  # 1×21 vector
                geom_matrix.append(vec)
            print(files)

            geom_matrix = np.array(geom_matrix)  # → shape (n_files, 21)
            with open(os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_geominfo.csv"), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                distances = np.round(geom_matrix[:, :6], 4)
                angles = np.round(geom_matrix[:, 6:], 1)
                geom_matrix_rounded = np.hstack([distances, angles])

                # Combine filenames + geometry rows
                rows = [
                    [fname] + list(row)
                    for fname, row in zip(files, geom_matrix_rounded)
                ]

                writer.writerows(rows)
        except Exception as e:
            print("ERROR", v)
            print(e)
    for k, v in g3.items():
        try:
            geom_matrix = []
            files = []
            for file in v:
                files.append(((os.path.basename(file)).split(".")[0]).strip("_new_align"))
                vec = geometry_information(file)  # 1×21 vector
                geom_matrix.append(vec)
            print(files)

            geom_matrix = np.array(geom_matrix)  # → shape (n_files, 21)
            with open(os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_geominfo.csv"), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)
                distances = np.round(geom_matrix[:, :6], 4)
                angles = np.round(geom_matrix[:, 6:], 1)
                geom_matrix_rounded = np.hstack([distances, angles])

                # Combine filenames + geometry rows
                rows = [
                    [fname] + list(row)
                    for fname, row in zip(files, geom_matrix_rounded)
                ]

                writer.writerows(rows)
        except Exception as e:
            print("ERROR", v)
            print(e)
    for k, v in g1.items():
        try:
            ordered = sorted(indices_for_group(v, key="state"))
            legend(ordered, available_stuff(ordered, STATES),
                   os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_legend.png"))
        except:
            print("ERROR", v)

    print("\n=== Same solvent + same state ===")
    for k, v in g2.items():
        try:
            ordered = sorted(indices_for_group(v, key="functional"))
            legend(ordered, available_stuff(ordered, FUNCTIONALS),
                   os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_legend.png"))
        except:
            print("ERROR", v)

    print("\n=== Same functional + same state ===")
    for k, v in g3.items():
        try:
            ordered = sorted(indices_for_group(v, key="solvent"))
            legend(ordered, available_stuff(ordered, SOLVENTS),
                   os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_legend.png"))
        except:
            print("ERROR", v)
    print("\n=== Same functional + same solvent ===")
    for k, v in g1.items():
        # CPK
        print(k,v)
        output_name = os.path.join(working_folder,f"overlap_{k[0]}_{k[1]}_CPK")
        run_pipeline(v, output_name,indices_for_group(v, key="state"))
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_CPK_y")
        run_pipeline(v, output_name, indices_for_group(v, key="state"), "VMD2/camera_y.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        # Lines
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_lines")
        run_pipeline(v, output_name, indices_for_group(v, key="state"), "VMD2/camera.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_lines_y")
        run_pipeline(v, output_name, indices_for_group(v, key="state"), "VMD2/camera_y.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )

    print("\n=== Same solvent + same state ===")
    for k, v in g2.items():
        #CPK
        output_name = os.path.join(working_folder,f"overlap_{k[0]}_{k[1]}_CPK")
        run_pipeline(v, output_name,indices_for_group(v, key="functional"))
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_CPK_y")
        run_pipeline(v, output_name, indices_for_group(v, key="functional"), "VMD2/camera_y.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        #Lines
        output_name = os.path.join(working_folder,f"overlap_{k[0]}_{k[1]}_lines")
        run_pipeline(v, output_name,indices_for_group(v, key="functional"), "VMD2/camera.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_lines_y")
        run_pipeline(v, output_name, indices_for_group(v, key="functional"), "VMD2/camera_y.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )



    print("\n=== Same functional + same state ===")
    for k, v in g3.items():
        output_name = os.path.join(working_folder,f"overlap_{k[0]}_{k[1]}_CPK")
        run_pipeline(v, output_name,indices_for_group(v, key="solvent"))
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_CPK_y")
        run_pipeline(v, output_name, indices_for_group(v, key="solvent"), "VMD2/camera_y.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        #Lines
        output_name = os.path.join(working_folder,f"overlap_{k[0]}_{k[1]}_lines")
        run_pipeline(v, output_name,indices_for_group(v, key="solvent"), "VMD2/camera.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )
        output_name = os.path.join(working_folder, f"overlap_{k[0]}_{k[1]}_lines_y")
        run_pipeline(v, output_name, indices_for_group(v, key="solvent"), "VMD2/camera_y.vmd","VMD2/molecule_lines.vmd")
        Image.open(f"{output_name}.tga").save(
            f"{output_name}.png"
        )




# VMD color map (RGB approximations)


end = time.perf_counter()
print(f"Elapsed: {end - start:.6f} seconds")