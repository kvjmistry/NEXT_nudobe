#!/bin/bash 
Model="mbb"
NM="SM"
PRESSURE="1bar"
mkdir -p merged

for f in $(ls ${Model}_${NME}_${PRESSURE}_nexus); do echo "" > merged/${Model}_${NME}_${PRESSURE}_nexus_reco_merged.txt;cat ${Model}_${NME}_${PRESSURE}_nexus/$f >> merged/${Model}_${NME}_${PRESSURE}_nexus_reco_merged.txt;done
for f in $(ls ${Model}_${NME}_${PRESSURE}_1mm_smear); do echo "" > merged/${Model}_${NME}_${PRESSURE}_1mm_smear_reco_merged.txt;cat ${Model}_${NME}_${PRESSURE}_1mm_smear/$f >> merged/${Model}_${NME}_${PRESSURE}_1mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_${NME}_${PRESSURE}_2mm_smear); do echo "" > merged/${Model}_${NME}_${PRESSURE}_2mm_smear_reco_merged.txt;cat ${Model}_${NME}_${PRESSURE}_2mm_smear/$f >> merged/${Model}_${NME}_${PRESSURE}_2mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_${NME}_${PRESSURE}_4mm_smear); do echo "" > merged/${Model}_${NME}_${PRESSURE}_4mm_smear_reco_merged.txt;cat ${Model}_${NME}_${PRESSURE}_4mm_smear/$f >> merged/${Model}_${NME}_${PRESSURE}_4mm_smear_reco_merged.txt;done
for f in $(ls ${Model}_${NME}_${PRESSURE}_10mm_smear); do echo "" > merged/${Model}_${NME}_${PRESSURE}_10mm_smear_reco_merged.txt;cat ${Model}_${NME}_${PRESSURE}_10mm_smear/$f >> merged/${Model}_${NME}_${PRESSURE}_10mm_smear_reco_merged.txt;done