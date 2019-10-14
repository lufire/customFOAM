#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest
# Christoph Goniva - July. 2011, mod by Alice Hager - July 2011
#===================================================================#
refineMeshByCellSet()
{
   while [ $# -ge 1 ]
   do
      echo "creating cell set for primary zone - $1"
      cp system/topoSetDict.$1 system/topoSetDict
      topoSet > log.topoSet.$1 2>&1

      echo "refining primary zone - $1"
      refineMesh -dict system/refineMeshDict -overwrite > log.refineMesh.$1 2>&1
      shift
   done
}

. $WM_PROJECT_DIR/bin/tools/RunFunctions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
runOctave="false"
postproc="false"

# check if mesh was built
#if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
#    echo "mesh was built before - using old mesh"
#else
#    echo "mesh needs to be built"
#    cd $casePath/CFD
#    runApplication blockMesh
#    #runApplication snappyHexMesh -overwrite
#    #refineMeshByCellSet 1
#    #echo "refineMeshByCellSet finish"
#fi

#if [ -f "$casePath/DEM/post/data.liggghts_init_dem" ];  then
#    echo "LIGGGHTS init was run before - using existing restart file"
#else
#    cd $casePath
#    #- run DEM in new terminal
#    $casePath/parDEMrun.sh
#fi

#-------------------------------------------------------#
# adapt settings for init run
#mkdir $casePath/CFD/0/
#cp -R $casePath/CFD/0.org/* $casePath/CFD/0/
#cp -R $casePath/CFD/constant/initRun/* $casePath/CFD/constant/
#cp -R $casePath/CFD/system/initRun/* $casePath/CFD/system/
#-------------------------------------------------------#
#-------------------------------------------------------#

#-------------------------------------------------------#
# additional pressure jump
#cp system/topoSetDict.fixedJump system/topoSetDict
#runApplication topoSet
#runApplication createBaffles -overwrite
#rm -R $casePath/CFD/0/cellLevel
#rm -R $casePath/CFD/0/pointLevel
#-------------------------------------------------------#

#- run parallel CFD-DEM in new terminal
cd $casePath/CFD
#gnome-terminal --title='cfdemSolverPiso ErgunTestMPI_restart CFD'  -e "bash $casePath/parCFDDEMrun.sh" 
#bash $casePath/parCFDDEMrun.sh
gdb cfdemElectrochemicalFoam 
wait

#- wait until sim has finished then run octave
#echo "simulation finished? ...press enter to proceed"
#read
#rm -R $casePath/CFD/processor*
#cp -R $casePath/CFD/ $casePath/CFD_initRun/
#cd $casePath/CFD/
#ls -dt */ | tail -n +2 | xargs rm -rf
#cd ..
#cp -R $casePath/CFD_initRun/constant $casePath/CFD/constant
#cp -R $casePath/CFD_initRun/system $casePath/CFD/system
#
#
##-------------------------------------------------------#
## adapt settings for init or restart run
#cp -R $casePath/CFD/constant/restartRun/* $casePath/CFD/constant/
#cp -R $casePath/CFD/system/restartRun/* $casePath/CFD/system/
#
##- run parallel CFD-DEM in new terminal
#cd $casePath
##gnome-terminal --title='cfdemSolverPiso ErgunTestMPI_restart CFD'  -e "bash $casePath/parCFDDEMrun.sh" 
#bash $casePath/parCFDDEMrun.sh


#- wait until sim has finished then run octave
#echo "simulation finished? ...press enter to proceed"
#read
#rm -R $casePath/CFD/processor*
#-------------------------------------------------------#

if [ $runOctave == "true" ]
  then
    #- change path
    cd $casePath/Octave

    #- run octave
    octave processParticleContacts.m

    #- show plot 
    #evince cfdemSolverPiso_ErgunTestMPI.eps

    #- copy log file to test harness
    #cp $casePath/../$logfileName $testHarnessPath
    #cp cfdemSolverPiso_ErgunTestMPI.eps $testHarnessPath
fi

if [ $postproc == "true" ]
  then

    #- keep terminal open (if started in new terminal)
    #echo "simulation finished? ...press enter to proceed"
    #read

    #- get VTK data from CFD sim
    cd $casePath/CFD
    foamToVTK                                                   #- serial run of foamToVTK

    #- get VTK data from liggghts dump file
    cd $casePath/DEM/post
    python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_restart

    #source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh                       #- include functions
    #pseudoParallelRun "foamToVTK" $nrPostProcProcessors          #- pseudo parallel run of foamToVTK

    #- start paraview
    #paraview
fi
