#!/bin/bash

echo "Starting Job" 

JOBID=$1
echo "The JOBID number is: ${JOBID}" 

JOBNAME=$2
echo "The JOBNAME number is: ${JOBNAME}" 

echo "JOBID $1 running on `whoami`@`hostname`"

MODEL=$3
echo "Model name is: ${MODEL}"

start=`date +%s`

# Setup nexus
echo "Setting Up NEXUS" 
source /software/nexus/setup_nexus.sh

# Set the configurable variables
N_EVENTS=100
CONFIG=ATPC_BSM.config.mac
INIT=ATPC_BSM.init.mac

echo "N_EVENTS: ${N_EVENTS}"

SEED=$((${N_EVENTS}*${JOBID} + ${N_EVENTS}))
echo "The seed number is: ${SEED}" 

# Change the config in the files
sed -i "s#.*random_seed.*#/nexus/random_seed ${SEED}#" ${CONFIG}
sed -i "s#.*start_id.*#/nexus/persistency/start_id ${SEED}#" ${CONFIG}
sed -i "s#.*output_file.*#/nexus/persistency/output_file ${MODEL}#" ${CONFIG}
sed -i "s#.*dist_file.*#/Generator/ElecPair/dist_file ${MODEL}.txt#" ${CONFIG}

# Print out the config and init files
cat ${INIT}
cat ${CONFIG}

# NEXUS
echo "Running NEXUS" 
nexus -n $N_EVENTS ${INIT}

# Loop over different bin sizes
for BIN in {1..5}
do
    python3 ${SCRIPT} ${MODEL} ${BIN}
done

ls -ltrh

tar -cvf *.h5 nudobe.tar

echo "FINISHED....EXITING" 

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds 