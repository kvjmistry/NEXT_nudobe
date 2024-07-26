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
# Model_NME_binning
JOBNAME="Leptoquark_SM"

BINNING="1mm_smear"
#BINNING="2mm_smear"
#BINNING="4mm_smear"
#BINNING="10mm_smear"
#BINNING="nexus"

H5FILE="/home/argon/Projects/Krishan/NEXT_nudobe/files/${JOBNAME}_${BINNING}.h5"
EVENTFILE=/home/argon/Projects/Krishan/NEXT_nudobe/files/${JOBNAME}_events.txt

echo "JOBNAME is: ${JOBNAME}"
echo "BINNING is: ${BINNING}"
echo "H5File is: ${H5FILE}"
echo "EVENTFILE is: ${EVENTFILE}"


# Set the configurable variables
N=100 # The number of segments to run -- 1-103 slurm

# Create the directory
mkdir -p "${JOBNAME}_${BINNING}"
cd "${JOBNAME}_${BINNING}"

# Setup VENV so we have python
echo "Setting Up Python" 
source /home/argon/Projects/Krishan/venv/bin/activate

cp /home/argon/Projects/Krishan/NEXT_nudobe/workdir/TrackReconstruction_functions.py .

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
# python3 /home/argon/Projects/Krishan/NEXT_nudobe/scripts/kinematics_reconstruction.py ${H5FILE} "segment_${SLURM_ARRAY_TASK_ID}.txt" "${JOBNAME}_${BINNING}_${SLURM_ARRAY_TASK_ID}"  
python3 /home/argon/Projects/Krishan/NEXT_nudobe/workdir/TrackReconstruction.py ${H5FILE} "segment_${SLURM_ARRAY_TASK_ID}.txt" "${JOBNAME}_${BINNING}_${SLURM_ARRAY_TASK_ID}"  

rm segment_${SLURM_ARRAY_TASK_ID}.txt
ls -ltrh

echo; echo; echo;

echo "FINISHED....EXITING" 

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds 