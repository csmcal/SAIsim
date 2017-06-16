# Writes the Condor DAG submission file for specified effect grid,
# Assumes the script runs 10 local instances

import numpy as np

submitFile = 'SAIsim_Full.sub'
outfile = open('SAIsim_FullParam.dag','w')

for mutPower in [-3,-4]:
	for invMutPower in [-3,-4]:
		for convPower in [-2]:
			jobName = 'm'+repr(-mutPower)+'i'+repr(-invMutFreq)+'c'+repr(-convPower)
			jobLine = 'JOB '+jobName+' '+submitFile+'\n'
			varsLine = 'VARS '+jobName+' mut=\"'+str(mutPower)+'\" inv=\"'+str(invMutFreq)\
				+'\" con=\"'+str(convPower)+'\"\n'
			outfile.write(jobLine)
			outfile.write(varsLine)

outfile.write('../CONFIG effectGridDAG.config\n')