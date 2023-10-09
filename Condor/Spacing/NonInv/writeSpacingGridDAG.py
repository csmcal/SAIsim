# Writes the Condor DAG submission file for specified effect grid
# Here the set of spacings between two mutations in a chromosome without inversions

import numpy as np

submitFile = 'SAIsim_TwoSpacing.sub'

smallSurvEffect = '0.86'
smallReprEffect = '0.11'

bigSurvEffect = '0.7'
bigReprEffect = '0.3'

#the script runs batches of 100, so n * 100 total
numBatches = 10

with open('SAIsim_TwoSpacingGrid.dag','w') as outfile:
	# Set the spacing, generally every 0.0125 or 0.025
	for spacing in np.around(np.linspace(0.0,0.75,num=61),8):
		# Add replicates for large sample sizes
		for batch in np.arange(numBatches):
			jobName = 'TS_'+'0'*(4-len(str(int(spacing*10000))))+str(int(spacing*10000))+"_"+str(batch)
			# print("Job "+jobName+":\n\tspacing\t"+str(spacing))
			jobLine = 'JOB '+jobName+' '+submitFile+'\n'
			varsLine = 'VARS '+jobName+' sur1=\"'+smallSurvEffect+'\" rep1=\"'+smallReprEffect+'\" sur2=\"'+\
				bigSurvEffect+'\" rep2=\"'+bigReprEffect+'\" spa=\"'+str(spacing)+'\" batchNum=\"'+str(batch)+'\"\n'
			outfile.write(jobLine)
			outfile.write(varsLine)

	outfile.write('CONFIG ../SAIsimDAG.config\n')