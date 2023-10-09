#!/bin/bash
#
# condor_python310.sh
# A bash wrapper for running any python 3.10 script from an htcondor job by simply passing arguments

# Information collection used to diagnose memory overallocation errors when un-tarring:
# ls -lh $_CONDOR_SCRATCH_DIR/..

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
