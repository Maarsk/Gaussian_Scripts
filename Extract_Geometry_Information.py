
from encodings.rot_13 import rot13
import math
from math import acos
import numpy as np
from itertools import combinations
import csv
import os
from Functions import out_folder_mkdir, read_xyz


def angle(a,b,c):
    r1 = ((matrice[a,0]-matrice[b,0])**2 + (matrice[a,1]-matrice[b,1])**2 + (matrice[a,2]-matrice[b,2])**2)**0.5
    r2 = ((matrice[b,0]-matrice[c,0])**2 + (matrice[b,1]-matrice[c,1])**2 + (matrice[b,2]-matrice[c,2])**2)**0.5
    r3 = ((matrice[a,0]-matrice[c,0])**2 + (matrice[a,1]-matrice[c,1])**2 + (matrice[a,2]-matrice[c,2])**2)**0.5
    theta = acos((r1**2+r2**2-r3**2)/(2*r1*r2))
    return math.degrees(theta)


input_folder = 'Input/xyz'
if not os.path.exists('Input'):
    os.mkdir('Input')
if not os.path.exists(input_folder):
    os.mkdir(input_folder)
    print("Input folder created.")
input('Put the .out files in Input/ex_out folder, then press Enter')
dist_int = np.zeros(6)
dic = np.array([4,6,13,15,7,1])
header_angles = np.array(('theta(Cl-Re-Cax/deg)', 'theta(Cl-Re-N1/deg)', 'theta(Cl-Re-N2/deg)', 'theta(Cl-Re-C1/deg)', 'theta(Cl-Re-C2/deg)', 'theta(Cax-Re-N1/deg)', 'theta(Cax-Re-N2/deg)', 'theta(Cax-Re-C1/deg)', 'theta(Cax-Re-C2/deg)', 'theta(N1-Re-N2/deg)',	'theta(N1-Re-C1/deg)', 'theta(N1-Re-C2/deg)', 'theta(N2-Re-C1/deg)', 'theta(N2-Re-C2/deg)', 'theta(C1-Re-C2/deg)'))
header_dist = np.array(('d(Re-Cl/A','d(Re-Cax/A','d(Re-N1/A','d(Re-N2/A','d(Re-C1/A','d(Re-C2/A'))


for filename in os.listdir(input_folder):
    if filename.lower().endswith(".xyz"):
        file_path = os.path.join(input_folder, filename)
        print(f"Elaborating: {file_path}")
        matrice = read_xyz(file_path)[0]
        dist = np.zeros(len(matrice))
        for i in range(len(matrice)):
            for j in range(3):
                dist[i] = dist[i]+(matrice[0,j]-matrice[i,j])**2
            dist[i]=dist[i]**0.5
        comb = list(combinations(dic,2))
        matrix_angles = np.zeros(len(comb))
        for i in range(len(comb)):
            matrix_angles[i] = angle(comb[i][0],0,comb[i][1])
        for i in range(6):
            dist_int[i]=dist[dic[i]]
        name, ext = os.path.splitext(filename)
        out_file = out_folder_mkdir(name,'Geom_Values_out','csv')
        with open(out_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header_angles)
            writer.writerow(np.round(matrix_angles,1))
            writer.writerow(header_dist)
            writer.writerow(np.round(dist_int,4))
            
input('Press enter to exit')