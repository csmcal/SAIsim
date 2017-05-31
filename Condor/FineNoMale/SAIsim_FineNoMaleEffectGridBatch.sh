#!/bin/bash
#
# SAIsim_EffectGridBatch.sh
# A wrapper for running the SAIsim_EffectGrid.py python script 10 times

# untar the python installation
tar -xzf python.tar.gz

# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
# (and a good home directory?)
mkdir home
export HOME=$(pwd)/home

for locNum in `seq 0 9`;
do
	# run python script
	# Assumes $1 = surEffect, $2 = repEffect, $3 = $PROCESS
	args=("$@")
	ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}''

	# echo ${ARGUMENTS}$locNum

	python SAIsim_FineNoMaleEffectGrid.py ${ARGUMENTS}$locNum
done
