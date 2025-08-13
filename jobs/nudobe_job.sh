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

# Set the configurable variables
CONFIG=ATPC_BSM.config.mac
INIT=ATPC_BSM.init.mac

echo "N_EVENTS: ${N_EVENTS}"

SEED=$((${JOBID} + 1))
echo "The seed number is: ${SEED}" 


if [ "$PRESS" -eq 1 ]; then
    
    N_EVENTS=120
    EID=$((${N_EVENTS}*${JOBID} + ${N_EVENTS}))
    echo "The EID number is: ${EID}" 

    # Change the config in the files
    sed -i "s#.*random_seed.*#/nexus/random_seed ${SEED}#" ${CONFIG}
    sed -i "s#.*start_id.*#/nexus/persistency/start_id ${EID}#" ${CONFIG}
    sed -i "s#.*output_file.*#/nexus/persistency/output_file ${MODEL}_${PRESS}bar#" ${CONFIG}
    sed -i "s#.*dist_file.*#/Generator/ElecPair/dist_file ${MODEL}_${NME}.txt#" ${CONFIG}
    sed -i "s#.*gas_pressure.*#/Geometry/ATPC/gas_pressure ${PRESS} bar#" ${CONFIG}
    sed -i "s#.*cube_size.*#/Geometry/ATPC/cube_size 6.182 m#" ${CONFIG}

    # Print out the config and init files
    cat ${INIT}
    cat ${CONFIG}

    # NEXUS
    echo "Running NEXUS" 
    nexus -n $N_EVENTS ${INIT}

    # Compress the file
    python3 CompressEvents.py ${MODEL}_${PRESS}bar

    # <Scale Factor> <CO2Percentage> <binsize> <pressure> <JOBID>
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 0 0.05 10 1.0 ${JOBID} # Just smearing
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1 0.05 10 1.0 ${JOBID} # Helium 10%
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1  0.1 10 1.0 ${JOBID} # 0.1 % CO2
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1 0.25 10 1.0 ${JOBID} # 0.25 % CO2
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1    5 10 1.0 ${JOBID} # 5.0 % CO2
    # python3 SmearEvents.py ${JOBNAME}_1bar 1    0 50 1.0 ${JOBID} # Pure Xe
    mv ${MODEL}_${PRESS}bar.h5 ${MODEL}_${PRESS}bar_nexus_${JOBID}.h5

elif [ "$PRESS" -eq 5 ]; then

    N_EVENTS=80
    EID=$((${N_EVENTS}*${JOBID} + ${N_EVENTS}))
    echo "The EID number is: ${EID}" 

    # Change the config in the files
    sed -i "s#.*random_seed.*#/nexus/random_seed ${SEED}#" ${CONFIG}
    sed -i "s#.*start_id.*#/nexus/persistency/start_id ${EID}#" ${CONFIG}
    sed -i "s#.*output_file.*#/nexus/persistency/output_file ${MODEL}_${PRESS}bar#" ${CONFIG}
    sed -i "s#.*dist_file.*#/Generator/ElecPair/dist_file ${MODEL}_${NME}.txt#" ${CONFIG}
    sed -i "s#.*gas_pressure.*#/Geometry/ATPC/gas_pressure ${PRESS} bar#" ${CONFIG}
    sed -i "s#.*cube_size.*#/Geometry/ATPC/cube_size 3.615 m#" ${CONFIG}

    # Print out the config and init files
    cat ${INIT}
    cat ${CONFIG}

    # NEXUS
    echo "Running NEXUS" 
    nexus -n $N_EVENTS ${INIT}

    # Compress the file
    python3 CompressEvents.py ${MODEL}_${PRESS}bar

    # <Scale Factor> <CO2Percentage> <binsize> <pressure> <JOBID>
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 0 0.05 10 5.0 ${JOBID} # Just smearing
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1 0.05 10 5.0 ${JOBID} # Helium 10%
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1  0.1 10 5.0 ${JOBID} # 0.1 % CO2
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1 0.25 10 5.0 ${JOBID} # 0.25 % CO2
    python3 SmearEvents.py ${MODEL}_${PRESS}bar 1    5 10 5.0 ${JOBID} # 5.0 % CO2
    # python3 SmearEvents.py ${JOBNAME}_1bar 1    0 50 1.0 ${JOBID} # Pure Xe
    mv ${MODEL}_${PRESS}bar.h5 ${MODEL}_${PRESS}bar_nexus_${JOBID}.h5

fi

ls -ltrh

echo "Taring the h5 files"
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