#!/bin/sh
#cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
echo $casePath

if [ -d "$casePath/CFD_initRun/" ]; then
    rm -R $casePath/CFD
    mv $casePath/CFD_initRun $casePath/CFD
fi

cd $casePath/CFD
cleanCase
rm -R 0
cd ..

rm -R $casePath/DEM/post
mkdir $casePath/DEM/post
rm -R $casePath/DEM/restart
mkdir $casePath/DEM/restart

rm -R $casePath/log*

# ----------------------------------------------------------------- end-of-file
