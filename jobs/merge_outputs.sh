#!/bin/bash 
Model="mbb"
NME="SM"
PRESSURE="15bar"
mkdir -p merged

# Define the smear values
smear_values=("1mm_smear" "2mm_smear" "4mm_smear" "10mm_smear")

# Loop over the smear values
for smear in "${smear_values[@]}"; do


    echo "${smear}"
    PATHS=/ospool/ap40/data/krishan.mistry/job/nudobe/${Model}/${NME}/${PRESSURE}/reco_v1/${smear}/

    # Initialize the merged file
    echo "" > merged/${Model}_${NME}_${PRESSURE}_${smear}_reco_merged.txt

    # Loop over files in the current smear directory
    for f in $(ls ${PATHS}); do
        # Append the contents of each file to the merged file
        cat ${PATHS}/$f >> merged/${Model}_${NME}_${PRESSURE}_${smear}_reco_merged.txt
    done
done
