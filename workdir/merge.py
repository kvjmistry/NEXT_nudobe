import tarfile
import os
from pathlib import Path
import glob 
import sys
import pandas as pd
import re

# USAGE: python merge.py <path to tars> <binsize>
# python extract.py "/ospool/ap40/data/krishan.mistry/job/nudobe/Leptoquark_SM/" 10mm

path  = sys.argv[1]
binsize=sys.argv[2]
print("Path:",     path)
print("Binsize:",  binsize)

files = glob.glob(path + "/" + binsize + "/*.h5")

# Given string
filename = os.path.basename(files[0])

# Use regex to remove the '_0' part and return the modified string
outfile = re.sub(r'_(\d+)_', '_', filename)

if (binsize == "nexus"):
    outfile = re.sub(r"_\d+", "_nexus", filename)

print("Filename Out:", outfile)

parts = []
hits = []

for f in files:
    print(f)

    part = pd.read_hdf(f, 'MC/particles')
    parts.append(part)

    # Same with the hits
    hit = pd.read_hdf(f, 'MC/hits')
    hits.append(hit)

parts = pd.concat(parts, ignore_index=True)
hits = pd.concat(hits, ignore_index=True)

print(parts)
print(hits)

# Open the HDF5 file in write mode
with pd.HDFStore(f"{path}/{outfile}", mode='w', complevel=5, complib='zlib') as store:
    # Write each DataFrame to the file with a unique key
    store.put('MC/particles',parts, format='table')
    store.put('MC/hits',hits, format='table')