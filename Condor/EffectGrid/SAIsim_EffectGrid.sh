#!/bin/bash
#
# SAIsim_EffectGrid.sh
# A wrapper for the SAIsim_EffectGrid.py python script

# untar the python installation and associated packages
tar -xzf python310.tar.gz
tar -xzf packages.tar.gz

# make sure the script will use the python installation,
# and the working directory as its home location
export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$PWD

# run python script
# Assumes arguments = "$(survEffect) $reprEffect) $(encounterNum) $(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"
args=("$@")
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}' '${args[6]}' '${args[7]}''

echo "SAIsim_EffectGrid.py "${ARGUMENTS}$batchNum
python SAIsim_EffectGrid.py ${ARGUMENTS}