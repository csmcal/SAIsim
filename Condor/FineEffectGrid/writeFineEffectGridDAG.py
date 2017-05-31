# Writes the Condor DAG submission file for specified effect grid,
# Assumes the script runs 10 local instances

import numpy as np

submitFile = 'SAIsim_FineEffectGrid.sub'
outfile = open('SAIsim_FineEffectGrid.dag','w')

# for survEffect in np.linspace(0.0,1.0,num=21):
for survEffect in np.linspace(0.70,1.0,num=31):
	for reprEffect in np.linspace(0.0,0.20,num=21):
		jobName = str(int(survEffect*100))+'_'+str(int(reprEffect*100))
		jobLine = 'JOB '+jobName+' '+submitFile+'\n'
		varsLine = 'VARS '+jobName+' sur=\"'+str(survEffect)+'\" rep=\"'+str(reprEffect)+'\"\n'
		outfile.write(jobLine)
		outfile.write(varsLine)

outfile.write('CONFIG effectGridDAG.config\n')