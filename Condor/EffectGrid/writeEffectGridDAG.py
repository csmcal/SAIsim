# Writes the Condor DAG submission file for input specified effect grid, or default 

import argparse
import numpy as np


# EXAMPLE SUBMIT-DESCRIPTION

# SUBMIT-DESCRIPTION std_effect_grid_jobDesc {
# 	universe = vanilla
# 	log = Log/effectGrid_$(Cluster).log
# 	error = Log/effectGrid_$(Cluster)_$(Process).err
# 	executable = SAIsim_EffectGrid_Batch.sh
# 	output = Log/effectGrid_$(Cluster)_$(Process).out
# 	should_transfer_files = YES
# 	when_to_transfer_output = ON_EXIT
# 	transfer_input_files = SAIsim_EffectGrid.py, SAIsim.py, http://proxy.chtc.wisc.edu/SQUID/chtc/el8/python310.tar.gz, ../packages.tar.gz
# 	request_cpus = 1
# 	request_memory = 256MB
# 	request_disk = 640MB
# 	arguments = "$(survEffect) $(reprEffect) $(encounterNum) $(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"
# }



# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()


    # Add optional arguments to the parser:

    # I/O, file location, flag paramters
    parser.add_argument('--out', type=str, default='SAIsim_EffectGrid.dag',
                        help='The name of the DAG file to create')
    parser.add_argument('--sub', type=str, default='effect_grid',
                        help='The name of the condor job description written within the DAG file')
    parser.add_argument('--exec', type=str, default='SAIsim_EffectGrid_Batch.sh', #default='SAIsim_EffectGrid.sh',
                        help='The bash script to run within each compute node, for DAG-defined jobs')
    parser.add_argument('--jobFileReqs', type=str, nargs='+',
    					default=[
    						'SAIsim_EffectGrid.py', 
							'SAIsim.py', 
							'http://proxy.chtc.wisc.edu/SQUID/chtc/el8/python310.tar.gz', 
							'../packages.tar.gz'
						],
                        help='The files required for the script to run')
    parser.add_argument('--test', default=False, action='store_true',
                        help='A flag telling the script to run only jobs with large memory and disk requests, '+\
                        'and only a few queued runs per job')
    parser.add_argument('--rep_N', type=int, default=100,
                        help='The number of replicate simulations to run per paramter combination')

    # Parameters relevant to the jobs
    parser.add_argument('--enc_N', type=int, default=100,
                        help='The numbers of male encounters per female to run sims with')
    # parser.add_argument('--c', type=float, default=[2,5], nargs='+',
    #                     help='The gene conversion rate powers to run sims on: conversion rates/locus = exp(-c)')
    parser.add_argument('--N', type=float, default=1e3,
                        help='The population size to run sims on')
    parser.add_argument('--g', type=float, default=2e4,
                        help='The number of generations to run sims for')

    parser.add_argument('--min_S', type=float, default=0,
                        help='The lowest survival value simulated')
    parser.add_argument('--max_S', type=float, default=0.975,
                        help='The highest survival value simulated')
    parser.add_argument('--step_S', type=float, default=0.025,
                        help='The difference between successive survival values')


    parser.add_argument('--min_R', type=float, default=0.025,
                        help='The lowest mate competition value simulated')
    parser.add_argument('--max_R', type=float, default=1,
                        help='The highest mate competition value simulated')
    parser.add_argument('--step_R', type=float, default=0.025,
                        help='The difference between successive mate competition values')

    parser.add_argument('--no_female_costs', default=False, action='store_true',
                        help='A flag telling the script to run without survival costs for females')
    parser.add_argument('--no_male_costs', default=False, action='store_true',
                        help='A flag telling the script to run without survival costs for males')


    args = parser.parse_args()
    return args

def getMemDiskLines(popSize):
	if testing:
		return ['request_memory = 4GB','request_disk = 4GB']
	elif popSize < 8000:
		return ['request_memory = 640MB','request_disk = 640MB']
	else:
		return ['request_memory = 1GB','request_disk = 1GB']



def writeEffectGridDag(outFileName,jobTypeName,jobDescStr,sur_values,rep_values):
	with open(outFileName,'w') as outfile:
		outfile.write('SUBMIT-DESCRIPTION '+jobTypeName+' {\n'+jobDescStr+'}\n')
			# s0.0r0.025e100N1000g20000nmc0nfc1n0
		jobNameInfix = 'e'+repr(encNum)+'N'+repr(popSize)+'g'+repr(numGens)+\
						'nmc'+repr(int(nmc))+'nfc'+repr(int(nfc))
		jobVarsInfix = '\" encounterNum=\"'+repr(encNum)+'\" popSize=\"'+repr(popSize)+'\" numGens=\"'+repr(numGens)+\
						'\" noMaleCosts=\"'+repr(int(nmc))+'\" noFemaleCosts=\"'+repr(int(nfc))
		for s in sur_values:
			for r in rep_values:
				for repNum in range(numReps):
					jobName = 's'+repr(s)+'r'+repr(r)+jobNameInfix+'n'+repr(repNum)
					jobLine = 'JOB '+ jobName + ' ' + jobTypeName + '\n'
					outfile.write(jobLine)
					varsLine = 'VARS '+jobName+' survEffect=\"'+repr(s)+'\" reprEffect=\"'+repr(r)+\
						jobVarsInfix+'\" repNum=\"'+repr(repNum)+'\"\n'
					outfile.write(varsLine)
		outfile.write('CONFIG ../SAIsimDAG.config\n')
	return


# Main function, for managing paramter checks and writing the dag
def main():
	# Parse arguments
	args = parse_args()

	outFileName = args.out
	submitName = args.sub
	jobExecutable = args.exec
	jobInputFileList = args.jobFileReqs

	global testing
	testing = args.test

	global numReps
	global encNum
	global popSize
	global numGens
	global nfc
	global nmc
	numReps = args.rep_N
	encNum = args.enc_N
	popSize = int(args.N)
	numGens = int(args.g)
	nfc = args.no_female_costs
	nmc = args.no_male_costs

	numSigDigS = len(str(args.step_S))-2
	sur_values = np.round(np.arange(args.min_S,args.max_S+args.step_S,args.step_S),numSigDigS)
	numSigDigR = len(str(args.step_R))-2
	rep_values = np.round(np.arange(args.min_R,args.max_R+args.step_R,args.step_R),numSigDigR)

	if testing and numReps > 10:
			numReps = 1
			print("Detected > 10 reps given for a testing flagged DAG, defaulting to 1 rep")

	jobDescLines = [
		'universe = vanilla',
		'log = Log/effectGrid_$(Cluster).log',
		'error = Log/effectGrid_$(Cluster)_$(Process).err',
		'executable = '+str(jobExecutable),
		'output = Log/effectGrid_$(Cluster)_$(Process).out',
		'should_transfer_files = YES',
		'when_to_transfer_output = ON_EXIT',
		'transfer_input_files = '+', '.join(jobInputFileList),
		'request_cpus = 1'
	]

	argumentsLine = 'arguments = "$(survEffect) $(reprEffect) $(encounterNum) '+\
		'$(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"'



	# SUBMIT-DESCRIPTION std_effect_grid_jobDesc {
	# 	universe = vanilla
	# 	log = Log/effectGrid_$(Cluster).log
	# 	error = Log/effectGrid_$(Cluster)_$(Process).err
	# 	executable = SAIsim_EffectGrid_Batch.sh
	# 	output = Log/effectGrid_$(Cluster)_$(Process).out
	# 	should_transfer_files = YES
	# 	when_to_transfer_output = ON_EXIT
	# 	transfer_input_files = SAIsim_EffectGrid.py, SAIsim.py, http://proxy.chtc.wisc.edu/SQUID/chtc/el8/python310.tar.gz, ../packages.tar.gz
	# 	request_cpus = 1
	# 	request_memory = 256MB
	# 	request_disk = 640MB
	# 	arguments = "$(survEffect) $(reprEffect) $(encounterNum) $(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"
	# }

	jobTypeName = submitName+'_jobDesc'
	jobDescStr = '\t'+'\n\t'.join(jobDescLines)+'\n\t'+'\n\t'.join(getMemDiskLines(popSize))+\
		'\n\t'+argumentsLine+'\n'

	writeEffectGridDag(outFileName,jobTypeName,jobDescStr,sur_values,rep_values)

	return



if __name__ == "__main__":
	main()
