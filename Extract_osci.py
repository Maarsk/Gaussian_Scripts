#!/usr/bin/env python
# coding: utf-8


import re
import csv
import os

def out_folder_mkdir(filename):
    if not os.path.exists('Output'):
        os.mkdir('Output')
    if not os.path.exists('Output/Oscillations_out'):
        os.mkdir('Output/Oscillations_out')
    return f"Output/Oscillations_out/excited_states_{filename}"

state_header_re = re.compile(
    r"Excited State\s+(\d+):.*?([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+f=([\d.]+)"
)
transition_re = re.compile(r"(\d+)\s*->\s*(\d+)")

input_folder = 'Input/ex_out/'
if not os.path.exists('Input'):
    os.mkdir('Input')
if not os.path.exists(input_folder):
    os.mkdir(input_folder)
    print("Input folder created.")
print('Put the .out files in Input/ex_out folder, then press Enter')
input()

for filename in os.listdir(input_folder):
    if filename.lower().endswith(".out"):
        file_path = os.path.join(input_folder, filename)
        print(f"Processing {filename}...")
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()
                start_indices = [
                    i for i, line in enumerate(lines)
                    if "Excitation energies and oscillator strengths" in line
                ]

                if not start_indices:
                    raise ValueError("Excitation energies not found")

                start_idx = start_indices[-1]
                block = lines[start_idx:]

                results = []
                current_state = None

                for line in block:
                    header_match = state_header_re.search(line)
                    if header_match:
                        if current_state:
                            results.append(current_state)

                        current_state = {
                            "State": int(header_match.group(1)),
                            "Energy_eV": float(header_match.group(2)),
                            "Wavelength_nm": float(header_match.group(3)),
                            "OscillatorStrength_f": float(header_match.group(4)),
                            "Transitions": []
                        }

                    elif current_state:
                        t_match = transition_re.findall(line)
                        if t_match:
                            for (a, b) in t_match:
                                current_state["Transitions"].append(f"{a}->{b}")

                if current_state:
                    results.append(current_state)

                name, ext = os.path.splitext(filename)
                output_file = out_folder_mkdir(name)
                with open(output_file + '.csv', "w", newline="") as csvfile:
                    fieldnames = ["State", "Energy_eV", "Wavelength_nm", "OscillatorStrength_f", "Transitions"]
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for r in results:
                        writer.writerow({
                            "State": r["State"],
                            "Energy_eV": r["Energy_eV"],
                            "Wavelength_nm": r["Wavelength_nm"],
                            "OscillatorStrength_f": r["OscillatorStrength_f"],
                            "Transitions": ";".join(r["Transitions"])
                        })
        except Exception as e:
            print(f"Error processing {filename}: {e}")

input('Press enter to exit')

