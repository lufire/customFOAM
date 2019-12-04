#!/bin/bash

#===================================================================#
# DEMrun script for ErgunTestMPI_restart testcase
# init ErgunTestMPI_restart
# Daniel Queteschiner - June 2014, DCS Computing GmbH
#===================================================================#

#- source CFDEM env vars
#. ~/.bashrc

#- include functions
#source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

parDEMrun()
{
    #--------------------------------------------------------------------------------#
    #- define variables
    logpath="$1"
    logfileName="$2"
    casePath="$3"
    headerText="$4"
    solverName="$5"
    nrProcs="$6"
    machineFileName="$7"
    debugMode="$8"
    #--------------------------------------------------------------------------------#

    if [ $debugMode == "on" ]; then
        debugMode="valgrind"
    elif [ $debugMode == "strict" ]; then
        #debugMode="valgrind --leak-check=full -v --trace-children=yes --track-origins=yes" 
        debugMode="valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes"  
    else
        debugMode=""
    fi

    #- clean up old log file
    rm $logpath/$logfileName

    #- change path
    cd $casePath/DEM

    #- header
    echo 2>&1 | tee -a $logpath/$logfileName
    echo "//   $headerText   //" 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- write path
    pwd 2>&1 | tee -a $logpath/$logfileName
    echo 2>&1 | tee -a $logpath/$logfileName

    #- run applictaion
    if [ $machineFileName == "none" ]; then
        mpirun -np $nrProcs $debugMode $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 2>&1 | tee -a $logpath/$logfileName
    else
        mpirun -machinefile $machineFileName -np $nrProcs $debugMode $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 2>&1 | tee -a $logpath/$logfileName
    fi

    #- keep terminal open (if started in new terminal)
    #read
}

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
logpath="$casePath"
caseName=${PWD##*/}
headerText="liggghts_init_$caseName"
logfileName="log_$headerText"
solverName="in.liggghts_init_dem"
nrProcs="4"
machineFileName="none"   # yourMachinefileName | none
debugMode="off"
#--------------------------------------------------------------------------------#

#- call function to run DEM case
#parDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode
cd $casePath/DEM
#mpirun -np 4 liggghts <in.liggghts_init_dem
mpirun -np $nrProcs $CFDEM_LIGGGHTS_SRC_DIR/$CFDEM_LIGGGHTS_LIB_NAME < $solverName 
