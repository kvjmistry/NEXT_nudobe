# Script to untar and merge

model="Leptoquark"
NME="SM"
PRESS=1bar

python3 extract.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "${model}" "${NME}"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "nexus"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "1mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "2mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "4mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}/${NME}/${PRESS}/" "10mm"