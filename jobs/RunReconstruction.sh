#!/bin/bash
#SBATCH -J GN # A single job name for the array
#SBATCH --nodes=1
#SBATCH --mem 4000 # Memory request (6Gb)
#SBATCH -t 0-24:00 # Maximum execution time (D-HH:MM)
#SBATCH -o GN_%A_%a.out # Standard output
#SBATCH -e GN_%A_%a.err # Standard error

start=`date +%s`

echo "The JOBID is ${SLURM_ARRAY_TASK_ID}" | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt
 
# Set the configurable variables
# Model_NME_binning
JOBNAME="Leptoquark_SM_nexus"
INFILE="/home/argon/Projects/Krishan/NEXT_nudobe/files/"

# Set the configurable variables
N=100 # The number of segments to run

# Create the directory
mkdir -p $JOBNAME/jobid_"${SLURM_ARRAY_TASK_ID}"
cd $JOBNAME/jobid_"${SLURM_ARRAY_TASK_ID}"

# Setup IC so we have python
echo "Setting Up IC" 2>&1 | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt
source /home/argon/Projects/Krishan/IC/setup_IC.sh

# Get the total number of lines in the file
infile=/home/argon/Projects/Krishan/NEXT_nudobe/files/${JOBNAME}_events.txt
total_lines=$(wc -l < ${infile})

# Calculate the size of each segment
segment_size=$((total_lines / N))

# Calculate the starting line for the nth segment, assume slurm ids start from 1
start_line=$((segment_size * $(($SLURM_ARRAY_TASK_ID-2)) + 1))

# Calculate the ending line for the 5th segment
end_line=$((segment_size * $(($SLURM_ARRAY_TASK_ID-1))))

# Extract the segment and save it to a new file so we can read this in for the job
sed -n "${start_line},${end_line}p" ${infile} > segment_${JOBID}.txt

# Run the reco
echo "Running Reco" 2>&1 | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt
python3 /home/argon/Projects/Krishan/NEXT_nudobe/scripts/kinematics_reconstruction.py "segment_${JOBID}.txt" "${JOBNAME}"  2>&1 | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt

ls -ltrh  2>&1 | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt

echo; echo; echo;

echo "FINISHED....EXITING" 2>&1 | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds | tee -a log_"${SLURM_ARRAY_TASK_ID}".txt