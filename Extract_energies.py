
import csv
import re
import os

input_folder = os.path.join(os.getcwd(),"Input","out_to_energies")
output_folder = os.path.join(os.getcwd(),"Output")
if not os.path.exists(input_folder):
    os.mkdir(input_folder)
    print("Input folder created.")

if not os.path.exists(output_folder):
    os.mkdir(output_folder)
    print("Output folder created.")
    input("Press ENTER to continue...")

free_energy_pattern = re.compile(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)")
output_file = os.path.join(output_folder,"Energies.csv")

with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['File', 'SCF Energy (Hartree)', 'TD Energy (Hartree)', 'Free Energy (Hartree)'])
for filename in os.listdir(input_folder):
    if filename.lower().endswith(".out"):
        file_path = os.path.join(input_folder, filename)
        print(f"Processing {filename}...")

        with open(file_path, 'r', errors='ignore') as f:
            for line in f:
                energy_pattern_td = re.compile(r"Total Energy.*?=\s*(-?\d+\.\d+)")
                match_td = energy_pattern_td.search(line)
                energy_pattern_scf = re.compile(r"SCF Done:\s+E\(.*?\)\s+=\s+(-?\d+\.\d+)")
                match_scf = energy_pattern_td.search(line)
                if match_td:
                    print(f"Found energy pattern in {filename}")
                    #break
        last_energy_td = None
        last_free_energy = None
        last_energy_scf = None

        with open(file_path, 'r', errors='ignore') as f:
            for line in f:
                match_energy_td = energy_pattern_td.search(line)
                match_energy_scf = energy_pattern_scf.search(line)
                if match_energy_td:
                    last_energy_td = float(match_energy_td.group(1))

                if match_energy_scf:
                    last_energy_scf = float(match_energy_scf.group(1))

                match_free = free_energy_pattern.search(line)
                if match_free:
                    last_free_energy = float(match_free.group(1))

        with open(output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([filename, last_energy_scf, last_energy_td, last_free_energy])

print(f"Extraction complete! Results saved to {output_file}")
print(os.path.splitext(output_file)[0])
