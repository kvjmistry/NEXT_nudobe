import tarfile
import os
from pathlib import Path
import glob 
import sys

# USAGE: python extract.py <path to tars> <model name> <NME>
# python extract.py "/ospool/ap40/data/krishan.mistry/job/nudobe/Leptoquark_SM/" "Leptoquark" "SM"

# Function to extract files from tar and place in the specified folders
def extract_files(tar_path, path, model, nme):
    tar_base_name = Path(tar_path).stem
    file_id = tar_base_name.split('_')[-1]

    # Define the dynamic folder mapping based on the tar file's id
    folder_mapping = {
        f'{model}_{nme}_{file_id}.h5': 'nexus',
        f'{model}_{nme}_{file_id}_1mm_smear.h5' :  '1mm',
        f'{model}_{nme}_{file_id}_2mm_smear.h5' :  '2mm',
        f'{model}_{nme}_{file_id}_4mm_smear.h5' :  '4mm',
        f'{model}_{nme}_{file_id}_10mm_smear.h5': '10mm'
    }

    with tarfile.open(tar_path, 'r:*') as tar:
        for member in tar.getmembers():
            filename = Path(member.name).name
            if filename in folder_mapping:
                folder = path+folder_mapping[filename]
                dest_path = Path(folder) / filename
                os.makedirs(folder, exist_ok=True)
                with tar.extractfile(member) as source, open(dest_path, 'wb') as dest:
                    dest.write(source.read())
                print(f"Extracted {filename} to {folder}")


path  = sys.argv[1]
model = sys.argv[2]
nme   = sys.argv[3]
print("Path:",  path)
print("Model:", model)
print("NME:",   nme)
 

tar_files = glob.glob(path+"/*.tar")
print(tar_files)

# Process each tar file
for tar_file in tar_files:
    extract_files(tar_file, path, model, nme)
