#!/bin/bash
#
# SAIsim_BigSmallBatch.sh
# A wrapper for running the SAIsim_BigSmall.py python script 10 times


# untar the python installation
tar -xzf python.tar.gz
# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
# (and a good home directory?)
mkdir home
export HOME=$(pwd)/home


# Old portion of the script for running replicates using the script, not the job submission
# for locNum in `seq 0 9`;
# do
# 	# run python script
# 	# Assumes $1 = surEffect1, $2 = repEffect1, $3 = surEffect2, $4 = repEffect2, $5 = spacing, $6 = $PROCESS
# 	args=("$@")
# 	ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}''

# 	# echo ${ARGUMENTS}$locNum

# 	python SAIsim_TwoSpacing.py ${ARGUMENTS}$locNum
# done


args=("$@")
# combine the batch and process number to the rep number
repNum=${args[6]}
while [ ${#repNum} -ne 2 ];
do
	repNum="0"$repNum
done
repNum=${args[5]}$repNum
# pass the arguments to the pyhton script
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '$repNum''
python SAIsim_TwoSpacingTest.py ${ARGUMENTS}