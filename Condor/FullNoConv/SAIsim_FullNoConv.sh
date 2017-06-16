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

# Make an output directory to compress for transfer back
DIRNAME='m'${args[0]}'i'${args[1]}'n'${args[3]}''
mkdir ${DIRNAME}
# mv SAIsim_Full.py ${DIRNAME}/
# mv SAIsim.py ${DIRNAME}/
cd ${DIRNAME}

# Run the script

# for locNum in `seq 0 9`;
# do
# 	# run python script
# 	# Assumes $1 = mutPower, $2 = invMutPower, $3 = convPower, $4 = $PROCESS
# 	args=("$@")
# 	ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}''

# 	# echo ${ARGUMENTS}$locNum

# 	python SAIsim_TwoSpacing.py ${ARGUMENTS}$locNum
# done

args=("$@")
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}''
python ../SAIsim_Full.py ${ARGUMENTS}

# Tar the result to be returned
cd ..
tar -czvf ${DIRNAME}.tar.gz ${DIRNAME}
