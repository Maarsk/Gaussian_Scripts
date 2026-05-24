import os
import subprocess
import tempfile
from PIL import Image
import re

VMD = r"C:\Program Files\University of Illinois\VMD2\vmd.exe"


def parse_filename(name):
    first_num = re.match(r"(\d+)", name).group(1)
    cleaned = name.replace(".cube", "").replace(".cub", "")
    middle = name.split("%")[1] if "%" in name else ""
    return first_num, cleaned, middle


def run_pipeline(working_folder, output_png):

    output_png = os.path.abspath(output_png).replace("\\", "/")
    mol_vmd = os.path.abspath("VMD/molecule_orbital_test.vmd").replace("\\", "/")
    iso_vmd = os.path.abspath("VMD/isosurface_test.vmd").replace("\\", "/")
    camera_vmd = os.path.abspath("VMD/camera.vmd").replace("\\", "/")

    cube_files = sorted(
        [f for f in os.listdir(working_folder) if f.endswith((".cube", ".cub"))]
    )

    with tempfile.NamedTemporaryFile(delete=False, suffix=".vmd") as tmp:
        script_path = tmp.name.replace("\\", "/")

        # FIRST CUBE → molecule + isosurface
        first = os.path.abspath(os.path.join(working_folder, cube_files[0])).replace("\\", "/")
        tmp.write(f'''
mol new "{first}" type cube waitfor all
molinfo top set center_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set rotate_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set scale_matrix  {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set global_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
set isovalue 0.03
set mode_orbital Transparent
set col_MO_plus 0
set col_MO_minus 1
source "{mol_vmd}"
'''.encode())

        # OTHER CUBES → isosurface only
        for cube in cube_files[1:]:
            cube_path = os.path.abspath(os.path.join(working_folder, cube)).replace("\\", "/")
            tmp.write(f'''
mol new "{cube_path}" type cube waitfor all
molinfo top set center_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set rotate_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set scale_matrix  {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
molinfo top set global_matrix {{{{1 0 0 0}} {{0 1 0 0}} {{0 0 1 0}} {{0 0 0 1}}}}
set isovalue 0.03
set mode_orbital Transparent
set col_MO_plus 0
set col_MO_minus 1
source "{iso_vmd}"
'''.encode())

        # APPLY CAMERA ONCE — AFTER ALL MOLECULES ARE LOADED
        tmp.write(f'''
source "{camera_vmd}"
'''.encode())

        # FINAL RENDER
        tmp.write(f'''
set outfile "{output_png}.tga"
render TachyonInternal $outfile
quit
'''.encode())

    subprocess.run([VMD, "-dispdev", "text", "-e", script_path])

if __name__ == "__main__":
    working_folder = os.path.abspath("VMD2/Output")

    # Build ONE image for ALL cubes
    output_name = "combined_isosurfaces"  # or whatever you want

    run_pipeline(
        working_folder=working_folder,
        output_png=os.path.join(working_folder, output_name)
    )

    # Convert TGA → PNG
    Image.open(f"{os.path.join(working_folder, output_name)}.tga").save(
        f"{os.path.join(working_folder, output_name)}.png"
    )

    print("[DONE] Combined image created.")

