#!/bin/bash
#
# perm2hk.sh
# A wrapper for the perm2hk R script running r/qtl 
#

# untar the R installation
tar -xzf R.tar.gz

# make sure the script will use the R installation
export PATH=$(pwd)/R/bin:$PATH

# run R
# Assumes $1 = $PROCESS, $2 = inputFile, $3 = phenList, $4 = numPerm, $5 = outFileID
args=("$@")
ARGUMENTS='--args '${args[0]}' '${args[1]}' '${args[2]}' '${args[3]}' '${args[4]}''
#echo ${ARGUMENTS}

#Rscript --vanilla perm2hk.R ${args[0]} ${args[1]} ${args[2]} ${args[3]} ${args[4]} &> Rlog${args[0]}.txt
# Since Rscript is broken by moving the install:
# R CMD BATCH --no-save --no-restore --quiet '--args jobNum='${args[0]}' inputFile='${args[1]}' pheVect='${args[2]}' numPerm='${args[3]}' outFileID='${args[4]}'' perm2hk.R Rlog${args[0]}.Rout &> Rlog${args[0]}.txt
R CMD BATCH --no-save --no-restore --quiet "$ARGUMENTS" perm2hk.R p2hkRlog_${args[4]}_${args[0]}.Rout