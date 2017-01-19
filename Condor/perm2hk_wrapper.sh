#!/bin/bash
#
# perm2hk_wrapper.sh
# Wraps the central script for a specific call:
#
# Does 1000 perms for every phenotype in goughF2_procSmRate.RData (the weights and simulated rates)

numJobs=100
numPerms=10
inputFile=goughF2_procSmRate.RData


for pheno in `seq 1 16`;
do
	outFileID=wk$pheno
	nohup ./perm2hk_central.sh $numJobs $inputFile $pheno $numPerms $outFileID &> nohupPerm2hk$outFileID.out&
	printf "$pheno $outFileID\n"
done


for pheno in `seq 17 32`;
do
	outFileID=rateSMwk$(($pheno - 16))
	nohup ./perm2hk_central.sh $numJobs $inputFile $pheno $numPerms $outFileID &> nohupPerm2hk$outFileID.out&
	printf "$pheno $outFileID\n"
done

#nohup ./perm2hk_central.sh $numJobs goughF2_genpr_v4.2.RData "4,5" 10 wk4-5 &> nohupPerm2hk$outFileID.out&
#nohup ./perm2hk_central.sh 2 goughF2_genpr_v4.2.RData 4 2 wk4test &> nohupPerm2hkwk4test.out&