# Writes the Condor DAG submission file for specified effect grid,
# Assumes the script runs 10 local instances

import numpy as np

submitFile = 'SAIsim_TwoSpacing.sub'

smallSurvEffect = '0.860'
smallReprEffect = '0.06'

bigSurvEffect = '0.7'
bigReprEffect = '0.15'

with open('SAIsim_TwoSpacingGrid.dag','w') as outfile:

	# for survEffect in np.linspace(0.0,1.0,num=21):
	for spacing in np.linspace(0.0,0.45,num=19):
	# for spacing in np.linspace(0.0,0.35,num=15):
		jobName = str(int(spacing*1000))
		jobLine = 'JOB '+jobName+' '+submitFile+'\n'
		varsLine = 'VARS '+jobName+' sur1=\"'+smallSurvEffect+'\" rep1=\"'+smallReprEffect+'\" sur2=\"'+\
			bigSurvEffect+'\" rep2=\"'+bigReprEffect+'\" spa=\"'+str(spacing)+'\"\n'
		outfile.write(jobLine)
		outfile.write(varsLine)

	outfile.write('CONFIG ../effectGridDAG.config\n')