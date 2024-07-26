#!/bin/bash 
Model="Leptoquark_SM"
mkdir -p merged

for f in $(ls ${Model}_nexus | grep *.txt); do cat ${Model}_nexus/$f >> merged/${Model}_nexus_reco_merged.txt;done
for f in $(ls ${Model}_1mm_smear | grep *.txt); do cat ${Model}_1mm_smear/$f >> merged/${Model}_1mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_2mm_smear | grep *.txt); do cat ${Model}_2mm_smear/$f >> merged/${Model}_2mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_4mm_smear | grep *.txt); do cat ${Model}_4mm_smear/$f >> merged/${Model}_4mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_10mm_smear | grep *.txt); do cat ${Model}_10mm_smear/$f >> merged/${Model}_10mm_smear_reco_merged.txt;done