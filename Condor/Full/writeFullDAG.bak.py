# Writes the Condor DAG submission file on input or default 

import argparse

# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()

    # # Create a new group to store required arguments
    # required_args = parser.add_argument_group('required named arguments')

    # # Add required arguments to the parser
    # required_args.add_argument('--count_file', type=str, required=True,
    #                     help='The .csv file containing the estimated inverted and standard read counts, etc.')
    #                     # ../inv_freqs/freq_data.csv


    # Add optional arguments to the parser:

    # I/O, file location, flag paramters
    parser.add_argument('--out', type=str, default='SAIsim_Full.dag',
                        help='The name of the DAG file to create')
    parser.add_argument('--sub', type=str, #default='SAIsim_Full.sub',
                        help='The condor job submission file to run in the DAG jobs')
    parser.add_argument('--exec', type=str, default='SAIsim_Full.sh',
                        help='The condor job submission file to run in the DAG jobs')
    parser.add_argument('--jobFileReqs', type=str, nargs='+',
    					default=[
    						'SAIsim_Full.py', 
							'SAIsim.py', 
							'http://proxy.chtc.wisc.edu/SQUID/chtc/python39.tar.gz', 
							'../packages.tar.gz'
						],
                        help='The condor job submission file to run in the DAG jobs')
    parser.add_argument('--test', default=False, action='store_true',
                        help='A flag telling the script to run only jobs with large memory and disk requests, '+\
                        'and only a few queued runs per job')
    parser.add_argument('--rep', type=int, default=1000,
                        help='The number of replicate simulations to run per paramter combination')

    # Parameters relevant to the jobs
    parser.add_argument('--m', type=float, default=[3,5], nargs='+',
                        help='The SA allele mutation powers to run sims on: mutation rates = exp(-m)')
    parser.add_argument('--i', type=float, default=[4], nargs='+',
                        help='The inversion mutation powers to run sims on: mutation rates = exp(-i)')
    parser.add_argument('--e', type=int, default=[10,100], nargs='+',
                        help='The numbers of male encounters per female to run sims with')
    parser.add_argument('--r', type=float, default=[1], nargs='+',
                        help='The crossover map lengths of the simulated region to run sims on')
    parser.add_argument('--c', type=float, default=[2,5], nargs='+',
                        help='The gene conversion rate powers to run sims on: conversion rates/locus = exp(-c)')
    parser.add_argument('--N', type=int, default=[3], nargs='+',
                        help='The population size powers to run sims on: size = exp(N)')

    args = parser.parse_args()
    return args


memDict = {
	2:"120MB",
	3:"120MB",
	4:"1GB",
	"test":"4GB"
}

diskDict = {
	2:"640MB",
	3:"640MB",
	4:"1GB",
	"test":"4GB"
}


	# '\tuniverse = vanilla\n'+\
	# '\tlog = Log/full_$(Cluster).log\n'+\
	# '\terror = Log/full_$(Cluster)_$(Process).err\n'+\
	# '\texecutable = SAIsim_Full.sh\n'+\
	# '\targuments = \"$(mut) $(inv) $(enc) $(rec) $(con) $(pop) $(Process)\"\n'+\
	# '\toutput = Log/full_$(Cluster)_$(Process).out\n'+\
	# '\trequirements = (OpSys == \"LINUX\")\n'+\
	# '\tshould_transfer_files = YES\n'+\
	# '\twhen_to_transfer_output = ON_EXIT\n'+\
	# '\ttransfer_input_files = SAIsim_Full.py, SAIsim.py, http://proxy.chtc.wisc.edu/SQUID/chtc/python39.tar.gz, ../packages.tar.gz\n'+\
	# '\trequest_cpus = 1\n'

def writeFullModelDag(outFileName,jobDescStr,numReps):
	with open(outFileName,'w') as outfile:
		for popSizePower in sizePowers:
			for mutPower in mutPowers:
				for invMutPower in invPowers:
					for encounterNum in encNums:
						for recRate in recRates:
							for convPower in conPowers:
									jobPrefix = 'm'+repr(mutPower)+'i'+repr(invMutPower)+\
										'e'+repr(encounterNum)+'r'+repr(recRate)+'c'+repr(convPower)+\
										'N'+repr(popSizePower)
									for repNum in range(numReps):
										jobName = jobPrefix + 'n' + repr(repNum)
										jobTypeName = jobName + 'jobDesc'
										if testing:
											memReq = memDict["test"]
											diskReq = diskDict["test"]
										else:
											memReq = memDict[popSizePower]
											diskReq = diskDict[popSizePower]
										addedJobLines = [
											'request_memory = '+memReq,
											'request_disk = '+diskReq,
											# 'queue '+str(numReps)
											'arguments = \"'+repr(mutPower)+' '+repr(invMutPower)+\
												' '+repr(encounterNum)+' '+repr(recRate)+' '+repr(convPower)+\
												' '+repr(popSizePower)+' '+repr(repNum)+'\"'
											]
										thisJobDesc = jobDescStr+'\t'+'\n\t'.join(addedJobLines)+'\n'
										jobTypeDescLine = 'SUBMIT-DESCRIPTION '+jobTypeName+' {\n'+thisJobDesc+'}\n'
										outfile.write(jobTypeDescLine)
										jobLine = 'JOB '+ jobName + ' ' + jobTypeName + '\n'
										# varsLine = 'VARS '+jobName+' mut=\"'+str(mutPower)+'\" inv=\"'+str(invMutPower)+\
										# 	'\" enc=\"'+str(encounterNum)+'\" rec=\"'+str(recRate)+'\" con=\"'+str(convPower)+\
										# 	'\" pop=\"'+str(popSizePower)+'\" rep=\"'+str(repNum)+'\"\n'
										outfile.write(jobLine)
										# outfile.write(varsLine)
		outfile.write('CONFIG ../SAIsimDAG.config\n')
	return

# Main function, for managing paramter checks and writing the dag
def main():
	# Parse arguments
	args = parse_args()

	outFileName = args.out
	submitFile = args.sub
	jobExecutable = args.exec
	jobInputFileList = args.jobFileReqs

	global testing
	testing = args.test

	numReps = args.rep

	global mutPowers
	global invPowers
	global encNums
	global recRates
	global conPowers
	global sizePowers
	mutPowers = args.m
	invPowers = args.i
	encNums = args.e
	recRates = args.r
	conPowers = args.c
	sizePowers = args.N

	if testing and numReps > 10:
			numReps = 1
			print("Detected > 10 reps given for a testing flagged DAG, defaulting to 1 rep")

	jobDescLines = [
		'universe = vanilla',
		'log = Log/full_$(Cluster).log',
		'error = Log/full_$(Cluster)_$(Process).err',
		'executable = '+str(jobExecutable),
		# 'arguments = \"$(mut) $(inv) $(enc) $(rec) $(con) $(pop) $(rep)\"',
		'output = Log/full_$(Cluster)_$(Process).out',
		# 'requirements = (OpSys == \"LINUX\")',
		'should_transfer_files = YES',
		'when_to_transfer_output = ON_EXIT',
		'transfer_input_files = '+', '.join(jobInputFileList),
		'request_cpus = 1'
	]

	jobDescStr = '\t'+'\n\t'.join(jobDescLines)+'\n'

	writeFullModelDag(outFileName,jobDescStr,numReps)

	return



if __name__ == "__main__":
	main()


