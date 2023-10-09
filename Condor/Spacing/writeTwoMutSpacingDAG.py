
# Writes the Condor DAG submission file for specified two-mutation spacing simulation,
# Assumes the submission file, executable script keep track of replicate numbers

import numpy as np
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
    parser.add_argument('--out', type=str, default='SAIsim_TwoSpacing.dag',
                        help='The name of the DAG file to create')
    parser.add_argument('--sub', type=str, default='SAIsim_TwoSpacing',
                        help='The submit-description name prefix used within the DAG to define jobs types')
    parser.add_argument('--exec', type=str, default='SAIsim_TwoSpacing.sh',
                        help='The condor job submission file to run in the DAG jobs')
    parser.add_argument('--jobFileReqs', type=str, nargs='+',
    					default=[
    						'SAIsim_TwoSpacing.py', 
							'SAIsim.py', 
							'http://proxy.chtc.wisc.edu/SQUID/chtc/el8/python310.tar.gz', 
							'../packages.tar.gz'
						],
                        help='The condor job submission file to run in the DAG jobs')
    parser.add_argument('--test', default=False, action='store_true',
                        help='A flag telling the script to run only jobs with large memory and disk requests, '+\
                        'and only a few queued runs per job')
    parser.add_argument('--rep', type=int, default=1000,
                        help='The number of replicate simulations to run per paramter combination')

    # Parameters distinguishing the jobs
    ## Small mutation spacing parameters
    parser.add_argument('--small_min_pos', type=float, default=0.0,
                        help='The minimum position of the smaller mutation')
    parser.add_argument('--small_max_pos', type=float, default=1.0,
                        help='The maximum position of the smaller mutation')
    parser.add_argument('--small_pos_change', type=float, default=0.005,
                        help='The change in position of the smaller mutation between simulations')



    # Additional simulation parameters
    ## Big mutation parameters
    parser.add_argument('--big_pos', type=float, default=0.225,
                        help='The position of the larger mutation')
    parser.add_argument('--big_sur', type=float, default=0.7,
                        help='The survival effect of the larger mutation')
    parser.add_argument('--big_rep', type=float, default=0.3,
                        help='The reproductive effect of the larger mutation')

    ## Small mutation parameters
    parser.add_argument('--small_sur', type=float, default=0.86,
                        help='The survival effect of the smaller mutation')
    parser.add_argument('--small_rep', type=float, default=0.12,
                        help='The reproductive effect of the smaller mutation')

    ## Inversion parameters
    parser.add_argument('--inv', default=False, action='store_true',
    # parser.add_argument('--inv', type=float, default=None, nargs=2,
                        help='A flag telling the script to include an inversion ')
    parser.add_argument('--inv_bp_1', type=float, default=0.21,
                        help='The minimum position of the smaller mutation')
    parser.add_argument('--inv_bp_2', type=float, default=0.55,
                        help='The maximum position of the smaller mutation')

    ## The initial frequency of all variants
    parser.add_argument('--freq', type=float, default=0.5,
                        help='The frequency at which the two SNPs and the inversion will start')

    ## Other
    parser.add_argument('--no_male_costs', action='store_true',
                        help='Run simulations where males experience no survival costs')
    parser.add_argument('--no_female_costs', action='store_true',
                        help='Run simulations where females experience no survival costs')
    parser.add_argument('--enc_N', type=int, default=100,
                        help='The number of males encountered by each female')
    parser.add_argument('--N', type=float, default=1e3,
                        help='The population size; default is 1000')
    parser.add_argument('--g', type=float, default=2e4,
                        help='The number of generations to run sims for, ideally 20Ne')

    args = parser.parse_args()
    return args



def memDiskCheck(popSize,numGens):
	if testing:
		return ("4GB","4GB")
	else:
		if popSize <= 2000:
			return ("640MB","640MB")
		elif popSize <= 10000:
			return ("1GB","1GB")
		else:
			return ("2GB","2GB")






# #the script runs batches of 100, so n * 100 total
# numBatches = 10

# with open('SAIsim_TwoSpacingGrid.dag','w') as outfile:
# 	# Set the spacing, generally every 0.0125 or 0.025
# 	for spacing in np.around(np.linspace(0.0,0.75,num=61),8):
# 		# Add replicates for large sample sizes
# 		for batch in np.arange(numBatches):
# 			jobName = 'TS_'+'0'*(4-len(str(int(spacing*10000))))+str(int(spacing*10000))+"_"+str(batch)
# 			# print("Job "+jobName+":\n\tspacing\t"+str(spacing))
# 			jobLine = 'JOB '+jobName+' '+submitFile+'\n'
# 			varsLine = 'VARS '+jobName+' sur1=\"'+smallSurvEffect+'\" rep1=\"'+smallReprEffect+'\" sur2=\"'+\
# 				bigSurvEffect+'\" rep2=\"'+bigReprEffect+'\" spa=\"'+str(spacing)+'\" batchNum=\"'+str(batch)+'\"\n'
# 			outfile.write(jobLine)
# 			outfile.write(varsLine)

# 	outfile.write('CONFIG ../SAIsimDAG.config\n')



def writeTwoSpacingDag(outFileName,jobDescPrefix,jobDescStr,numReps):
	with open(outFileName,'w') as outfile:
		# Write the job description
		jobTypeName = jobDescPrefix + '_jobDesc'
		(memReq, diskReq) = memDiskCheck(popSize,numGens)
		addedJobLines = [
			'request_memory = '+memReq,
			'request_disk = '+diskReq,
			# 'queue '+str(numReps)
			'arguments = \"$(smallPos) $(smallSurvEffect) $(smallReprEffect) '+\
				'$(bigPos) $(bigSurvEffect) $(bigReprEffect) '+\
				'$(inv) $(breakpoint1) $(breakpoint2) '+\
				'$(initFreq) '+\
				'$(encounterNum) $(popSize) $(numGens) '+\
				'$(noMaleCosts) $(noFemaleCosts) $(repNum)\"'
			# 'arguments = \"'+\
			# 	repr(survEffect)+' '+repr(reprEffect)+\
			# 	' '+repr(encounterNum)+' '+repr(popSize)+' '+repr(numGens)+\
			# 	' '+repr(int(maleCosts))+' '+repr(int(femaleCosts))+\
			# 	' '+repr(repNum)+'\"'
			]
		thisJobDesc = jobDescStr+'\t'+'\n\t'.join(addedJobLines)+'\n'
		jobTypeDescLine = 'SUBMIT-DESCRIPTION '+jobTypeName+' {\n'+thisJobDesc+'}\n'
		outfile.write(jobTypeDescLine)

		# Write individual jobs
		positionRange = np.around(np.arange(smallMinPos,smallMaxPos+smallChangePos,smallChangePos),8)
		print(positionRange)
		for smallPos in positionRange:

			jobPrefix = 'sp'+repr(smallPos)+'ss'+repr(smallSurvEffect)+'sr'+repr(smallReprEffect)+\
				'bp'+repr(bigPos)+'bs'+repr(bigSurvEffect)+'br'+repr(bigReprEffect)+\
				'i'+repr(int(inv))+'i1'+repr(breakpoint1)+'i2'+repr(breakpoint2)+\
				'f'+repr(initFreq)+\
				'e'+repr(encounterNum)+'N'+repr(popSize)+'g'+repr(numGens)+\
				'nmc'+repr(int(noMaleCosts))+'nfc'+repr(int(noFemaleCosts))

			for repNum in range(numReps):
				jobName = jobPrefix + 'n' + repr(repNum)

				jobLine = 'JOB '+ jobName + ' ' + jobTypeName + '\n'
				outfile.write(jobLine)

				varsLine = 'VARS '+jobName+' '+\
					'smallPos=\"'+str(smallPos)+'\" smallSurvEffect=\"'+str(smallSurvEffect)+'\" smallReprEffect=\"'+str(smallReprEffect)+'\" '+\
					'bigPos=\"'+str(bigPos)+'\" bigSurvEffect=\"'+str(bigSurvEffect)+'\" bigReprEffect=\"'+str(bigReprEffect)+'\" '+\
					'inv=\"'+str(int(inv))+'\" breakpoint1=\"'+str(breakpoint1)+'\" breakpoint2=\"'+str(breakpoint2)+'\" '+\
					'initFreq=\"'+str(initFreq)+'\" '+\
					'encounterNum=\"'+str(encounterNum)+'\" popSize=\"'+str(popSize)+'\" numGens=\"'+str(numGens)+'\" '+\
					'noMaleCosts=\"'+str(int(noMaleCosts))+'\" noFemaleCosts=\"'+str(int(noFemaleCosts))+'\" repNum=\"'+str(repNum)+'\"\n'
				outfile.write(varsLine)

		outfile.write('CONFIG ../SAIsimDAG.config\n')
	return


# Main function, for managing paramter checks and writing the dag
def main():
	# Parse arguments
	args = parse_args()

	outFileName = args.out
	jobDescPrefix = args.sub
	jobExecutable = args.exec
	jobInputFileList = args.jobFileReqs

	global testing
	testing = args.test

	numReps = args.rep

	global smallMinPos
	global smallMaxPos
	global smallChangePos
	global smallSurvEffect
	global smallReprEffect
	smallMinPos = args.small_min_pos
	smallMaxPos = args.small_max_pos
	smallChangePos = args.small_pos_change
	smallSurvEffect = args.small_sur
	smallReprEffect = args.small_rep

	global bigPos
	global bigSurvEffect
	global bigReprEffect
	bigPos = args.big_pos
	bigSurvEffect = args.big_sur
	bigReprEffect = args.big_rep

	global inv
	global breakpoint1
	global breakpoint2
	inv = args.inv
	breakpoint1 = args.inv_bp_1
	breakpoint2 = args.inv_bp_2
	assert(breakpoint2>breakpoint1), "Inversion breakpoints should be given with the earlier break first"

	global initFreq
	initFreq = args.freq

	global noMaleCosts
	global noFemaleCosts
	global encounterNum
	global popSize
	global numGens
	noMaleCosts = args.no_male_costs
	noFemaleCosts = args.no_female_costs
	encounterNum = args.enc_N
	popSize = int(args.N)
	numGens = int(args.g)

	if testing and numReps > 10:
			numReps = 1
			print("Detected > 10 reps given for a testing flagged dag, defaulting to 1 rep")

	jobDescLines = [
		'universe = vanilla',
		'log = Log/twoSpacing_$(Cluster).log',
		'error = Log/twoSpacing_$(Cluster)_$(Process).err',
		'executable = '+str(jobExecutable),
		# 'arguments = \"$(mut) $(inv) $(enc) $(rec) $(con) $(pop) $(rep)\"',
		'output = Log/twoSpacing_$(Cluster)_$(Process).out',
		# 'requirements = (OpSys == \"LINUX\")',
		'should_transfer_files = YES',
		'when_to_transfer_output = ON_EXIT',
		'transfer_input_files = '+', '.join(jobInputFileList),
		'request_cpus = 1'
	]

	jobDescStr = '\t'+'\n\t'.join(jobDescLines)+'\n'

	writeTwoSpacingDag(outFileName,jobDescPrefix,jobDescStr,numReps)

	return



if __name__ == "__main__":
	main()
