#!/bin/bash
#
# SAIsim_Full_Arrangements.sh
# A bash wrapper for running the SAIsim_Full_Arrangements.py python script

# untar the python installation
tar -xzf python310.tar.gz
# untar the python packages
tar -xzf packages.tar.gz

# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$(pwd)

[ -d Checkpoints ] || mkdir Checkpoints

# Make an output directory to compress for transfer back
args=("$@")
# printf -v checkFile '../'${args[0]}
# printf -v checkTime "%06d" ${args[1]}
printf -v mut "%04.3e" ${args[2]}
printf -v inv "%04.3e" ${args[3]}
printf -v encounter "%03d" ${args[4]}
printf -v crossover "%04.3e" ${args[5]}
printf -v conversion "%04.3e" ${args[6]}
printf -v pop "%02.1e" ${args[7]}
printf -v gen "%02.1e" ${args[8]}
printf -v wf "%01d" ${args[9]}
printf -v nmc "%01d" ${args[10]}
printf -v nfc "%01d" ${args[11]}
printf -v replicate "%04d" ${args[12]}
DIRNAME='mut'$mut'inv'$inv'enc'$encounter'cro'$crossover'con'$conversion'pop'$pop'gen'$gen'wf'$wf'nmc'$nmc'nfc'$nfc'rep'$replicate''
[ -d "${DIRNAME}" ] || mkdir ${DIRNAME}
cd ${DIRNAME}

# Run the script
# Python script assumes the following arguments (indexed 1 higher than shell scripting)
#            $1 = checkpointFilePath, $2 = timeToCheckpoint,  
#            $3 = mutRate, $4 = invMutRate, $5 = encounterNum, 
#            $6 = crossoverRate, $7 = convRate, $8 = popSize, $9 = numGens,
#            $10 = wrightFisherRules, $11 = noMaleCost, $12 = noFemaleCost,
#            $13 = replicateNum
ARGUMENTS='../'${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}' '${args[5]}' '${args[6]}' '${args[7]}' '${args[8]}' '${args[9]}' '${args[10]}' '${args[11]}' '${args[12]}''
echo 'python ../SAIsim_Full_Arrangements.py '${ARGUMENTS}
python ../SAIsim_Full_Arrangements.py ${ARGUMENTS}

# Check for checkpoint exit from the python script
exit_status=$?
if [ "${exit_status}" -eq 85 ]; then
	exit 85
fi

# If not checkpointing, tar the result to be returned
cd ..
tar -czvf ${DIRNAME}.tar.gz ${DIRNAME}
