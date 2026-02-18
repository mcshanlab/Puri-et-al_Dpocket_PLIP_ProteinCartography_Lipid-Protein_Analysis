import os
import re

# Folder containing the files
folder = "."  # current directory, change if needed

for filename in os.listdir(folder):
    old_path = os.path.join(folder, filename)
    if os.path.isfile(old_path):
        # Find the part starting with BD and stopping before the first "_"
        match = re.search(r"(BD[^_]*)", filename)
        if match:
            new_name = match.group(1) + os.path.splitext(filename)[1]  # keep original extension
            new_path = os.path.join(folder, new_name)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_name}")

print("Renaming complete.")