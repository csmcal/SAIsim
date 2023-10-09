# Writes the Condor DAG submission file for specified effect grid,
# Assumes the sub, script keep track of replicate numbers

import numpy as np

submitFile = 'SAIsim_EffectGrid.sub'
outfile = open('SAIsim_EffectGrid.dag','w')

with open('SAIsim_EffectGrid.dag','w') as outfile:
	# Set the effect spacing
	for survEffect in np.around(np.linspace(0.05,1.0,num=20),8):
		for reprEffect in np.around(np.linspace(0.0,1.0,num=21),8):
			# jobName = 'EG_'+'{:.4f}'.format(survEffect)+'_'+\
			# 	'{:.4f}'.format(reprEffect)
			jobName = 'EG_'+'{:04d}'.format(int(survEffect*1000))+'_'+\
				'{:04d}'.format(int(reprEffect*1000))
			jobLine = 'JOB '+jobName+' '+submitFile+'\n'
			varsLine = 'VARS '+jobName+' sur=\"'+str(survEffect)+\
				'\" rep=\"'+str(reprEffect)+'\"\n'
			outfile.write(jobLine)
			outfile.write(varsLine)

	outfile.write('CONFIG ../SAIsimDAG.config\n')