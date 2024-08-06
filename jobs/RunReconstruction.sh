#!/bin/bash
#SBATCH -J GN # A single job name for the array
#SBATCH --nodes=1
#SBATCH --mem 4000 # Memory request (6Gb)
#SBATCH -t 0-24:00 # Maximum execution time (D-HH:MM)
#SBATCH -o log/GN_%A_%a.out # Standard output
#SBATCH -e log/GN_%A_%a.err # Standard error

start=`date +%s`

echo "The JOBID is ${SLURM_ARRAY_TASK_ID}" 
 
# Set the configurable variables
MODEL="mbb"
#MODEL="Leptoquark"

NME="SM"
#NME="QRPA"
#NME="IBM2"

PRESSURE=1bar
#PRESSURE=5bar
#PRESSURE=10bar
#PRESSURE=15bar

BINNING="1mm_smear"
#BINNING="2mm_smear"
#BINNING="4mm_smear"
#BINNING="10mm_smear"
#BINNING="nexus"

H5FILE="/home/argon/Projects/Krishan/NEXT_nudobe/files/${MODEL}/${NME}/${PRESSURE}/${MODEL}_${BINNING}.h5"
EVENTFILE=/home/argon/Projects/Krishan/NEXT_nudobe/files//${MODEL}/${NME}/${PRESSURE}/${MODEL}_events.txt

echo "Model is: ${MODEL}"
echo "NME is: ${NME}"
echo "BINNING is: ${BINNING}"
echo "H5File is: ${H5FILE}"
echo "EVENTFILE is: ${EVENTFILE}"


# Set the configurable variables
N=100 # The number of segments to run -- 1-103 slurm

# Create the directory
mkdir -p "${MODEL}_${NME}_${PRESSURE}_${BINNING}"
cd "${MODEL}_${NME}_${PRESSURE}_${BINNING}"

# Setup VENV so we have python
echo "Setting Up Python" 
source /home/argon/Projects/Krishan/venv/bin/activate

cp /home/argon/Projects/Krishan/NEXT_nudobe/workdir/TrackReconstruction_functions.py ./TrackReconstruction_functions_${SLURM_ARRAY_TASK_ID}.py

# Get the total number of lines in the file
total_lines=$(wc -l < ${EVENTFILE})

# Calculate the size of each segment
segment_size=$((total_lines / N))

# Calculate the starting line for the nth segment, assume slurm ids start from 1
start_line=$((segment_size * $(($SLURM_ARRAY_TASK_ID-1)) + 1))
echo "Start line is: ${start_line}"

# Calculate the ending line for the 5th segment
end_line=$((segment_size * $(($SLURM_ARRAY_TASK_ID))))
echo "End line is: ${end_line}"

# Extract the segment and save it to a new file so we can read this in for the job
sed -n "${start_line},${end_line}p" ${EVENTFILE} > segment_${SLURM_ARRAY_TASK_ID}.txt

# Run the reco
echo "Running Reco" 
# python3 /home/argon/Projects/Krishan/NEXT_nudobe/scripts/kinematics_reconstruction.py ${H5FILE} "segment_${SLURM_ARRAY_TASK_ID}.txt" "${MODEL}_${NME}_${BINNING}_${SLURM_ARRAY_TASK_ID}"  
python3 TrackReconstruction_functions_${SLURM_ARRAY_TASK_ID}.py ${H5FILE} "segment_${SLURM_ARRAY_TASK_ID}.txt" "${MODEL}_${NME}_${PRESSURE}_${BINNING}_${SLURM_ARRAY_TASK_ID}"  

rm segment_${SLURM_ARRAY_TASK_ID}.txt
rm TrackReconstruction_functions_${SLURM_ARRAY_TASK_ID}.py
ls -ltrh

echo; echo; echo;

echo "FINISHED....EXITING" 

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds 