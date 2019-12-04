#!/bin/sh
#cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
cd $casePath/CFD
foamToVTK
cd $casePath/DEM/post
python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_init
python -i $CFDEM_LPP_DIR/lpp.py  dump.liggghts_restart
cd $casePath

