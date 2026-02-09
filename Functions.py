import csv
import os
import numpy as np

def out_folder_mkdir(filename,sub_folder,extention='txt'):
    if not os.path.exists('Output'):
        os.mkdir('Output')
    if not os.path.exists(f"Output/{sub_folder}"):
        os.mkdir(f"Output/{sub_folder}")
    return f"Output/{sub_folder}/{filename}.{extention}"
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