# Writes the Condor DAG submission file for specified effect grid,
# Assumes the script runs 10 local instances

import numpy as np

submitFile = 'SAIsim_NeutralInversion.sub'
outfile = open('SAIsim_NeutralInversionGridSm.dag','w')

for invMutFreq in np.linspace(-2.5,-1,num=4):
	jobName = repr(-invMutFreq)
	jobLine = 'JOB '+jobName+' '+submitFile+'\n'
	varsLine = 'VARS '+jobName+' mut=\"'+str(invMutFreq)+'\"\n'
	outfile.write(jobLine)
	outfile.write(varsLine)

outfile.write('CONFIG ../effectGridDAG.config\n')