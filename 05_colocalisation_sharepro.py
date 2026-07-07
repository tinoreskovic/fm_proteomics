import os

ids = ["ADAMTS15", "ADM", "AGER", "APOD", "CALB1", "CCL27", "CCN5", "CD300LG", "CLMP", "CRIM1", "DNER", "GHRL", "NEFL", "PRRT3", "PSPN", "ROBO1", "RTN4R", "SCGB3A1", "SCGB3A2", "SHBG", "WFIKKN2"]


base_dir = "./ShareProColoc"
script_path = os.environ.get("SHAREPRO_SCRIPT", "SharePro_coloc/src/SharePro/sharepro_coloc.py")
python_executable = os.environ.get("SHAREPRO_PYTHON", "python3")


for id in ids:
    z_file1 = os.path.join(base_dir, f"{id}.txt")
    z_file2 = os.path.join(base_dir, f"{id}_fat_mass.txt")
    ld_file = os.path.join(base_dir, f"{id}_ld.ld")
    save_path = os.path.join(base_dir, f"{id}")

    # Construct the SharePro command
    command = f'{python_executable} {script_path} --z "{z_file1}" "{z_file2}" --ld "{ld_file}" --save "{save_path}" --eps 0.0001'

    print(f"Running command: {command}")

    # Execute the command
    os.system(command)

