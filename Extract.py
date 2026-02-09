#%%
import csv
import re
import os

input_folder = 'Input/energies_to_out/'

free_energy_pattern = re.compile(r"Sum of electronic and thermal Free Energies=\s+(-?\d+\.\d+)")

output_file = "energy.csv"

with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['File', 'SCF Energy (Hartree)', 'Free Energy (Hartree)'])
for filename in os.listdir(input_folder):
    if filename.lower().endswith(".out"):
        file_path = os.path.join(input_folder, filename)
        print(f"Processing {filename}...")
        if not os.path.exists(input_folder):
            os.mkdir(input_folder)
            print("Input folder created.")
            input("Press ENTER to continue...")
            exit()
        with open(file_path, 'r', errors='ignore') as f:
            for line in f:
                energy_pattern = re.compile(r"Total Energy.*?=\s*(-?\d+\.\d+)")
                match = energy_pattern.search(line)
                if match:
                    print(line)
                    break
                energy_pattern = re.compile(r"SCF Done:\s+E\(.*?\)\s+=\s+(-?\d+\.\d+)")
            print(energy_pattern)
        last_energy = None
        last_free_energy = None

        with open(file_path, 'r', errors='ignore') as f:
            for line in f:
                match_energy = energy_pattern.search(line)
                if match_energy:
                    last_energy = float(match_energy.group(1))
                    print(last_energy)

                match_free = free_energy_pattern.search(line)
                if match_free:
                    last_free_energy = float(match_free.group(1))

        with open(output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([filename, last_energy, last_free_energy])

print(f"\n✅ Extraction complete! Results saved to {output_file}")

