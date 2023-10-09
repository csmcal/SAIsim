#!/bin/bash
#
# SAIsim_EffectGridBatchTest.sh
# A wrapper for the SAIsim_EffectGridTest.py python script

# untar the Python installation, here version 3.8, packages separately
tar -xzf python38.tar.gz
tar -xzf packages.tar.gz

# make sure the script will use your Python installation, 
# and the working directory as its home location
export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$PWD

# run python script
# Assumes $1 = surEffect, $2 = repEffect, $3 = $PROCESS
args=("$@")
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}''
#echo ${ARGUMENTS}

python SAIsim_EffectGrid.py ${ARGUMENTS}