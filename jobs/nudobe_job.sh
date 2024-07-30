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

start=`date +%s`

# Setup nexus
echo "Setting Up NEXUS" 
source /software/nexus/setup_nexus.sh


# Re-source this stuff for now. There seems to be a G4 installation on some machines
export G4INSTALL=/software/geant4-v11.1.0/install;
export PATH=$G4INSTALL/bin:$PATH;
export LD_LIBRARY_PATH=$G4INSTALL/lib:$LD_LIBRARY_PATH;
cd $G4INSTALL/bin; source geant4.sh; cd -;


# Set the configurable variables
N_EVENTS=25
CONFIG=ATPC_BSM.config.mac
INIT=ATPC_BSM.init.mac

echo "N_EVENTS: ${N_EVENTS}"

SEED=$((${N_EVENTS}*${JOBID} + ${N_EVENTS}))
echo "The seed number is: ${SEED}" 

# Change the config in the files
sed -i "s#.*random_seed.*#/nexus/random_seed ${SEED}#" ${CONFIG}
sed -i "s#.*start_id.*#/nexus/persistency/start_id ${SEED}#" ${CONFIG}
sed -i "s#.*output_file.*#/nexus/persistency/output_file ${MODEL}_${JOBID}#" ${CONFIG}
sed -i "s#.*dist_file.*#/Generator/ElecPair/dist_file ${MODEL}_${NME}.txt#" ${CONFIG}
sed -i "s#.*gas_pressure.*#/Geometry/ATPC/gas_pressure ${PRESS} bar#" ${CONFIG}

# Print out the config and init files
cat ${INIT}
cat ${CONFIG}

# NEXUS
echo "Running NEXUS" 
nexus -n $N_EVENTS ${INIT}

# Compress the file
python3 CompressEvents.py ${MODEL}_${JOBID}

# Loop over different bin sizes
for BIN in {1,2,4,10}
do
    echo "Running with Bin: $BIN mm"
    python3 SmearEventsGeneral.py ${MODEL}_${JOBID} 0 0.4 0.4 ${BIN}
done

ls -ltrh

tar -cvf nudobe.tar *.h5

# Cleanup
rm *.h5
rm *.mac
rm *.txt
rm *.py

ls -ltrh

echo "FINISHED....EXITING" 

end=`date +%s`
let deltatime=end-start
let hours=deltatime/3600
let minutes=(deltatime/60)%60
let seconds=deltatime%60
printf "Time spent: %d:%02d:%02d\n" $hours $minutes $seconds 