#!/bin/bash
#
# condor_pyscripts_checkpoint.sh
# A bash wrapper for running any python script with condor checkpointing by simply passing arguments

# untar the python installation
tar -xzf python310.tar.gz
# untar the python packages
tar -xzf packages.tar.gz

# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$(pwd)

[ -d Checkpoints ] || mkdir Checkpoints

# # Make an output directory to compress for transfer back
# args=("$@")
# # printf -v checkFile '../'${args[0]}
# # printf -v checkTime "%06d" ${args[1]}
# printf -v mut "%04.3e" ${args[2]}
# printf -v inv "%04.3e" ${args[3]}
# printf -v encounter "%03d" ${args[4]}
# printf -v crossover "%04.3e" ${args[5]}
# printf -v conversion "%04.3e" ${args[6]}
# printf -v pop "%02.1e" ${args[7]}
# printf -v gen "%02.1e" ${args[8]}
# printf -v wf "%01d" ${args[9]}
# printf -v nmc "%01d" ${args[10]}
# printf -v nfc "%01d" ${args[11]}
# printf -v replicate "%04d" ${args[12]}
# DIRNAME='mut'$mut'inv'$inv'enc'$encounter'cro'$crossover'con'$conversion'pop'$pop'gen'$gen'wf'$wf'nmc'$nmc'nfc'$nfc'rep'$replicate''
# [ -d "${DIRNAME}" ] || mkdir ${DIRNAME}
# cd ${DIRNAME}

# Run the script
ARGUMENTS="$@"
echo 'python '${ARGUMENTS}
python ${ARGUMENTS}

# Check for checkpoint exit from the python script
exit_status=$?
if [ "${exit_status}" -eq 85 ]; then
	exit 85
fi

# # If not checkpointing, tar the result to be returned
# cd ..
# tar -czvf ${DIRNAME}.tar.gz ${DIRNAME}
