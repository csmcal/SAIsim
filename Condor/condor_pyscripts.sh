#!/bin/bash
#
# condor_pyscripts.sh
# A bash wrapper for running any python script from an htcondor job by simply passing arguments

# untar the python installation
tar -xzf python310.tar.gz
# untar the python packages
tar -xzf packages.tar.gz

# make sure the script will use the python installation
export PATH=$(pwd)/python/bin:$PATH
export PYTHONPATH=$PWD/packages
export HOME=$(pwd)

# Run the script
ARGUMENTS="$@"
echo 'python '${ARGUMENTS}
python ${ARGUMENTS}
