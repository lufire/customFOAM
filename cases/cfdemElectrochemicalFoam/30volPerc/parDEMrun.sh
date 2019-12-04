#!/bin/bash

#===================================================================#
# DEMrun script for ErgunTestMPI_restart testcase
# init ErgunTestMPI_restart
# Daniel Queteschiner - June 2014, DCS Computing GmbH
#===================================================================#

#- source CFDEM env vars
#. ~/.bashrc

#- include functions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

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
parDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode

