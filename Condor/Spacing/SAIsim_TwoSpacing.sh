#!/bin/bash
#
# SAIsim_BigSmallBatch.sh
# A wrapper for running the SAIsim_BigSmall.py python script 10 times

# untar the Python installation and packages separately
tar -xzf python310.tar.gz
tar -xzf packages.tar.gz

# make sure the script will use your Python installation, 
# and the working directory as its home location
export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$PWD

# run python script
# Assumes arguments = $(smallPos) $(smallSurvEffect) $(smallReprEffect)
# 					$(bigPos) $(bigSurvEffect) $(bigReprEffect)
# 					$(inv) $(breakpoint1) $(breakpoint2)
# 					$(lowFreq)
# 					$(encounterNum) $(popSize) $(numGens)
# 					$(noMaleCosts) $(noFemaleCosts) $(repNum)
args=("$@")
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}' '${args[6]}' '${args[7]}' '${args[8]}' '${args[9]}' '${args[10]}' '${args[11]}' '${args[12]}' '${args[13]}' '${args[14]}' '${args[15]}''

echo "SAIsim_TwoSpacing.py "${ARGUMENTS} #$batchNum
python SAIsim_TwoSpacing.py ${ARGUMENTS}


# Scratch for running replicates using the script, not the job submission
# for batchNum in `seq 0 9`;
# do
# 	# run python script
# 	# Assumes $1 = surEffect1, $2 = repEffect1, $3 = surEffect2, $4 = repEffect2, $5 = spacing, $6 = $PROCESS
# 	args=("$@")
# 	ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}''
#	echo "SAIsim_TwoSpacing.py "${ARGUMENTS}$batchNum
#	python SAIsim_TwoSpacing.py ${ARGUMENTS}
# done


# args=("$@")
# # combine the batch and process number to the rep number
# repNum=${args[6]}
# while [ ${#repNum} -ne 2 ];
# do
# 	repNum="0"$repNum
# done
# repNum=${args[5]}$repNum
# # pass the arguments to the pyhton script
# ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '$repNum''
# python SAIsim_TwoSpacing.py ${ARGUMENTS}


