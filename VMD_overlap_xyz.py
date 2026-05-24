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
    camera_vmd = os.path.abspath("VMD2/camera.vmd").replace("\\", "/")
    molecules_vmd = os.path.abspath("VMD2/molecule.vmd").replace("\\", "/")
    cube_files = sorted(
        [f for f in os.listdir(working_folder) if f.endswith((".xyz"))]
    )

    with tempfile.NamedTemporaryFile(delete=False, suffix=".vmd") as tmp:
        script_path = tmp.name.replace("\\", "/")

        # FIRST CUBE → molecule + isosurface
        for i in range(len(cube_files)):
            first = os.path.abspath(os.path.join(working_folder, cube_files[i])).replace("\\", "/")
            tmp.write(f'''
    set xyz_path "{first}"
    source "{molecules_vmd}"
    '''.encode())
            # APPLY CAMERA

        tmp.write(f'''
    set mol_id {i}
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

