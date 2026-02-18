import os
import shutil

# Path to the main directory containing subfolders
source_dir = "/Volumes/GigiMurin/plip/sterol_lipids"

# Get the current directory
destination_dir = os.getcwd()

# Walk through all subfolders
for root, dirs, files in os.walk(source_dir):
    for file in files:
        if file.endswith(".pdb"):
            source_file = os.path.join(root, file)
            dest_file = os.path.join(destination_dir, file)
            
            # If a file with the same name exists, rename to avoid overwriting
            if os.path.exists(dest_file):
                base, ext = os.path.splitext(file)
                counter = 1
                while os.path.exists(os.path.join(destination_dir, f"{base}_{counter}{ext}")):
                    counter += 1
                dest_file = os.path.join(destination_dir, f"{base}_{counter}{ext}")
            
            shutil.copy2(source_file, dest_file)
            print(f"Copied: {source_file} -> {dest_file}")

print("All .pdb files copied.")