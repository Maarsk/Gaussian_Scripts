import os
import numpy as np
import re
from Functions import read_xyz
files = []
input_folder = os.path.join(os.getcwd(),"Input","orca_xyz")
output_folder = os.path.join(os.getcwd(),"Output")

energy_pattern = re.compile(r"E\s+(-?\d+\.\d+)")

for filename in os.listdir(input_folder):
    info = filename.split(".")
    if len(info) == 3:
        num = info[-2]
        files.append((num,filename))
        name = info[0]
files.sort(key=lambda x: x[0])
energies = []
with open(os.path.join(output_folder, f"{name}_movie.xyz"), "w") as f1:
    for num,filename in files:
        with open(os.path.join(input_folder,filename),"r") as f:
            for line in f:
                if line.startswith("Coordinates from ORCA-job"):
                    match = energy_pattern.search(line)
                    if match:
                        energy = float(match.group(1))
                        energies.append(energy)
                        break
        coords,elements = read_xyz(os.path.join(input_folder,filename))
        f1.write(f"{len(elements)}\n")
        f1.write(f"Frame {num}, Energy {energy}\n")
        for i in range(len(elements)):
            f1.write(f"{elements[i]:<3} {coords[i][0]:>15.8f} {coords[i][1]:>15.8f} {coords[i][2]:>15.8f}\n")
with open(os.path.join(output_folder, f"{name}_energies.csv"), "w") as f2:
    f2.write("Frame,Energy\n")
    for i in range(len(energies)):
        f2.write(f"{i+1},{energies[i]}\n")