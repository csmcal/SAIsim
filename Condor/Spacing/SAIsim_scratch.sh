#!/bin/bash
#
# SAIsim_BigSmallBatch.sh
# A wrapper for running the SAIsim_BigSmall.py python script 10 times

# untar the python installation
# tar -xzf python.tar.gz

# make sure the script will use the python installation
# export PATH=$(pwd)/python/bin:$PATH
# (and a good home directory?)
# mkdir home
# export HOME=$(pwd)/home

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
repNum=${args[6]}
echo $repNum
while [ ${#repNum} -ne 3 ];
do
	repNum="0"$repNum
done
echo $repNum
repNum=${args[5]}$repNum
echo $repNum
ARGUMENTS=''${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}''
echo $ARGUMENTS
ARGUMENTS1=${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}
echo $ARGUMENTS1
if [ "$ARGUMENTS" = "$ARGUMENTS1" ];
then
	echo 'Equal'
	echo "$ARGUMENTS"
else
	echo 'not equal'
	echo "$ARGUMENTS"
fi
ARGUMENTS2=${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '$repNum
echo $ARGUMENTS2

# python SAIsim_TwoSpacingInv.py ${ARGUMENTS}