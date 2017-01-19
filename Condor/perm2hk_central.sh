#!/bin/bash
#
# perm2hk_central.sh
# Runs on the submit node to run the submit 
# Also a wrapper for the perm_comb R script running r/qtl 
#
# Use as ./perm2hk_central.sh <numJobs> <inputFile> <phenList> <batchSize> <outFileID>
#


die () {
    echo >&2 "$@"
    exit 1
}

# Check that the right number of arguments is submitted
[ "$#" -ge 5 ] || die "at least 5 arguments (numJobs,input,phenList,numPermsPerJob,output) required, $# provided"

# Run the relevant r/qtl submission, with arguments
args=("$@")

#numAdditionalRuns=$((${args[0]} - 1))
#printf "$numAdditionalRuns\n"

# Assumes $1 = numJobs, $2 = inputFile, $3 = phenList, $4 = numPermsPerJob, $5 = outFileID
ARGUMENTS="arguments = \"\$(Process) ${args[1]} '${args[2]}' ${args[3]} ${args[4]}\""
#QUEUE="queue $numAdditionalRuns"
#com_out=$(condor_submit -append "$ARGUMENTS" -append "$QUEUE" perm2hk.sub)
com_out=$(condor_submit -append "$ARGUMENTS" perm2hk.sub -queue "${args[0]}")


# Wait until the condor queue is done, relies on predictable condor_submit output
CLUSTER=${com_out%.*}
CLUSTER=${CLUSTER##*er }
echo "$CLUSTER"
condor_wait perm2hk_$CLUSTER.log


# untar the R installation (here assumes untarred)
#tar -xzf R.tar.gz

# make sure the script will use the R installation
export PATH=$(pwd)/R/bin:$PATH

# run R combination script, Assumes $1 = numJobs, $5 = outFileID
R CMD BATCH --no-save --no-restore --quiet '--args '${args[0]}' '${args[4]}'' perm_comb.R comb_Rlog_$CLUSTER.Rout

# (alternate non-working calls:)
#Rscript --vanilla perm_comb.R ${args[0]} ${args[4]}
#R CMD BATCH --vanilla '--args numJobs='${args[0]}' outFileID='${args[4]}'' perm_comb.R RlogComb${args[0]}.Rout

# Remove separate RData files
# for pheno in `seq 0 $((${args[0]} - 1))`;
# do
# 	rm perm.${args[4]}.batch$jobNum.RData
# 	rm p2hkRlog_${args[4]}_$jobNum.Rout
#	
# done

# paste0("perm.${args[4]}.batch$jobNum.RData")

# Notify user of completion
echo "Condor Submission Cluster ${CLUSTER}; Type Perm2hk; ID ${args[4]} Completed" | mail -s "Condor ${CLUSTER} Completed" cmcallester@gmail.com


