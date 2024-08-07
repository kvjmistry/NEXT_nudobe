#!/bin/bash 
Model="mbb"
NM="SM"
PRESSURE="1bar"
mkdir -p merged

# Define the smear values
smear_values=("1mm_smear" "2mm_smear" "4mm_smear" "10mm_smear")

# Loop over the smear values
for smear in "${smear_values[@]}"; do

    # Initialize the merged file
    echo "" > merged/${Model}_${NME}_${PRESSURE}_${smear}_reco_merged.txt

    # Loop over files in the current smear directory
    for f in $(ls ${Model}_${NME}_${PRESSURE}_${smear}); do
        # Append the contents of each file to the merged file
        cat ${Model}_${NME}_${PRESSURE}_${smear}/$f >> merged/${Model}_${NME}_${PRESSURE}_${smear}_reco_merged.txt
    done
done