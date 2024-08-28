#!/bin/bash

echo "Starting Job" 

JOBID=$1
echo "The JOBID number is: ${JOBID}" 

JOBNAME=$2
echo "The JOBNAME number is: ${JOBNAME}" 

echo "JOBID $1 running on `whoami`@`hostname`"

MODEL=$3
echo "Model name is: ${MODEL}"

NME=$4
echo "NME is: ${NME}"

PRESS=$5
echo "Pressure is: ${PRESS}"

BINNING=$6
echo "Binning is: ${BINNING}"

# Set the configurable variables
N_JOBS=$7
echo "Segement Size is: ${N_JOBS}"

start=`date +%s`

H5FILE=${MODEL}_${BINNING}.h5
EVENTFILE=${MODEL}_events.txt

echo "Model is: ${MODEL}"
echo "NME is: ${NME}"
echo "BINNING is: ${BINNING}"
echo "H5File is: ${H5FILE}"
echo "EVENTFILE is: ${EVENTFILE}"


# Get the total number of lines in the file
total_lines=$(wc -l < ${EVENTFILE})

# Calculate the size of each segment
segment_size=$((total_lines / N_JOBS))

# Calculate the starting line for the nth segment
start_line=$((segment_size * $JOBID + 1))
echo "Start line is: ${start_line}"

# Calculate the ending line for the 5th segment
end_line=$((segment_size * $(($JOBID+1))))
echo "End line is: ${end_line}"

# Extract the segment and save it to a new file so we can read this in for the job
sed -n "${start_line},${end_line}p" ${EVENTFILE} > segment_${JOBID}.txt

# Run the reco
echo "Running Reco" 
# python3 TrackReconstruction.py ${H5FILE} "segment_${JOBID}.txt" "${MODEL}_${NME}_${PRESS}_${BINNING}"  
python3 kinematics_reconstruction.py ${H5FILE} "segment_${JOBID}.txt" "${MODEL}_${NME}_${PRESS}_${BINNING}" 

ls -ltrh

echo; echo; echo;

rm segment_${JOBID}.txt

# Check for the exit file
if [ ! -e "${MODEL}_${NME}_${PRESS}_${BINNING}_reco.txt" ]; then
  echo "Error: File does not exist, returning with STATUS 1."
  exit 1
fi



echo "FINISHED....EXITING" 

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds 
exit 0