#!/bin/sh
#cd ${0%/*} || exit 1    # run from this directory

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/CleanFunctions
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
echo $casePath

rm -R $casePath/log*
rm -R $casePath/CFD/processor*
rm -R $casePath/CFD/log*

# ----------------------------------------------------------------- end-of-file
