# Script to untar and merge

model="Leptoquark"
NME="SM"

python3 extract.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "${model}" "${NME}"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "nexus"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "1mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "2mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "4mm"
python3 merge.py "/ospool/ap40/data/krishan.mistry/job/nudobe/${model}_${NME}/" "10mm"