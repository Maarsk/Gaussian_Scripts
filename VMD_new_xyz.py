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
    files = os.path.abspath(files).replace("\\", "/")
    cube_files = [files]

    with tempfile.NamedTemporaryFile(delete=False, suffix=".vmd") as tmp:
        script_path = tmp.name.replace("\\", "/")
        tmp.write(f'''
    color change rgb 0 0 0.5 1        
    set xyz_path "{files}"
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

if __name__ == "__main__":
    working_folder = "VMD2/Output/Docker/ACN"

    for filename in os.listdir(working_folder):
        if filename.endswith(".xyz"):
            fullpath = os.path.join(working_folder, filename)
            # CPK
            output_name = os.path.join(working_folder, f"{filename.split(".")[0]}_CPK")
            run_pipeline(fullpath, output_name, [0],camera_path = "VMD2/camera.vmd", molecule_path = "VMD2/molecule.vmd")
            Image.open(f"{output_name}.tga").save(
                f"{output_name}.png"
            )
            output_name = os.path.join(working_folder, f"{filename.split(".")[0]}_CPK_y")
            run_pipeline(fullpath, output_name, [0], "VMD2/camera_y.vmd",molecule_path = "VMD2/molecule.vmd")
            Image.open(f"{output_name}.tga").save(
                f"{output_name}.png"
            )


end = time.perf_counter()
print(f"Elapsed: {end - start:.6f} seconds")