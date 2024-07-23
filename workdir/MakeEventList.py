import pandas as pd
import sys

# USAGE: python MakeEventList.py <infile> <outfile>
# python MakeEventList.py "/ospool/ap40/data/krishan.mistry/job/nudobe/Leptoquark_SM/Leptoquark_SM_nexus.h5"  "Leptoquark_SM_events.txt"

infile  = sys.argv[1]
outfile = sys.argv[2]
print("Infile:",  infile)
print("Outfile:", outfile)

parts = pd.read_hdf(infile, "MC/particles")

events = sorted(parts.event_id.unique())

# Open the file in write mode
with open(outfile, 'w') as file:
    for index, item in enumerate(events):

        if (index == len(events)-1):
            file.write(f"{item}") 
        else:  
            file.write(f"{item}\n") 
