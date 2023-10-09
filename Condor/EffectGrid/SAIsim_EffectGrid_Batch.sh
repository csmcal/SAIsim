#!/bin/bash
#
# SAIsim_EffectGridBatch.sh
# A wrapper for running the SAIsim_EffectGrid.py python script 10 times

# untar the python installation
tar -xzf python310.tar.gz
tar -xzf packages.tar.gz


# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
export PYTHONPATH=$PWD/packages
# (and a good home directory?)
export HOME=$(pwd)

for batchNum in `seq 0 9`;
do
	# run python script
	##### Assumes $1 = surEffect, $2 = repEffect, $3 = $PROCESS
	# Assumes arguments = "$(survEffect) $reprEffect) $(encounterNum) $(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"
	args=("$@")
	ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}' '${args[6]}' '${args[7]}''

	echo "SAIsim_EffectGrid.py "${ARGUMENTS}$batchNum

	python SAIsim_EffectGrid.py ${ARGUMENTS}$batchNum
done
