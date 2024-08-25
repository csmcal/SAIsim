
# Writes the Condor DAG submission file for specified two-mutation spacing simulation,
# Assumes the submission file, executable script keep track of replicate numbers

import numpy as np
import argparse

# Define and collect input arguments
def parse_args():
    parser = argparse.ArgumentParser()

    # I/O, file location, flag paramters
    parser.add_argument('--nomp', default=False, action='store_true',
                        help='A flag telling the script to run without parallel processing')

    # Parameters distinguishing the jobs
    ## Basic replications
    parser.add_argument('--rep', type=int, default=1000,
                        help='The number of replicate simulations to run per paramter combination')

    ## Survival effect value to use
    parser.add_argument('--sur', type=float, default=1.0,
                        help='The survival effect of the mutation simulated')

    ## Reproductive quality spacing parameters
    parser.add_argument('--rep_min', type=float, default=0.0,
                        help='The minimum reproductive quality value')
    parser.add_argument('--rep_max', type=float, default=2.0,
                        help='The maximum reproductive quality value')
    parser.add_argument('--rep_step', type=float, default=0.1,
                        help='The change in reproductive quality value')

    ## Frequency spacing parameters
    parser.add_argument('--freq_min', type=float, default=0.05,
                        help='The minimum mutation frequency')
    parser.add_argument('--freq_max', type=float, default=1.0,
                        help='The maximum mutation frequency')
    parser.add_argument('--freq_step', type=float, default=0.05,
                        help='The change in mutation frequency')

    ## Spacing parameters for proportion of heterozygote numbers out of the maximum possible at the allele frequency
    parser.add_argument('--ph_min', type=float, default=1.0,
                        help='The minimum proportion of heterozygotes out of the maximum possible at the allele frequency')
    parser.add_argument('--ph_max', type=float, default=1.0,
                        help='The maximum proportion of heterozygotes out of the maximum possible at the allele frequency')
    parser.add_argument('--ph_step', type=float, default=1.0,
                        help='The change in proportion of heterozygotes out of the maximum possible at the allele frequency')

    ## Encounter numbers
    parser.add_argument('--enc_N', type=int, nargs='+', default=[2,3,4,10,100],
                        help='The number of males encountered by each female')


    # Additional simulation parameters
    parser.add_argument('--no_male_costs', action='store_true',
                        help='Run simulations where males experience no survival costs')
    parser.add_argument('--no_female_costs', action='store_true',
                        help='Run simulations where females experience no survival costs')
    parser.add_argument('--N', type=float, default=1e3,
                        help='The population size; default is 1000')
    # parser.add_argument('--N', type=float, nargs='+', default=[1e3,1e4,1e5],
    #                     help='The population size; default is 1000')

    args = parser.parse_args()
    return args


# For running all sims locally, either sequentially or
# by founding a multiprocessing pool to asynchronously call the simulator
# Sims print their results, rather than returning them
def runSingMutPatGridPool():

	# Prepare the multiprocessing pool, if multiprocessing is enabled
	if run_mp:
		import multiprocessing as mp
		num_proc = mp.cpu_count()//2
		pool = mp.Pool(num_proc)
		# async_objects = []

	# Run each combination of parameters
	from SAIsim_SingMutPaternityGrid import runSim
	repSet = np.around(np.linspace(rep_min,rep_max,
		int(np.around((rep_max-rep_min)/rep_step,8))+1),8)
	print("Reproductive qualities")
	print(repSet)
	freqSet = np.around(np.linspace(freq_min,freq_max,
		int(np.around((freq_max-freq_min)/freq_step,8))+1),8)
	print("Allele frequencies")
	print(freqSet)
	propHetSet = np.around(np.linspace(ph_min,ph_max,
		int(np.around((ph_max-ph_min)/ph_step,8))+1),8)
	print("Proportions of heterozygotes")
	print(propHetSet)
	print("Encounter numbers")
	print(encNums)
	for repEffect in repSet:
		for freq in freqSet:
			for propHet in propHetSet:
				for encN in encNums:
					for repN in np.arange(numReps):
						if run_mp:
							# async_objects += [pool.apply_async(runSim,
							# 	args=(mutPos,surEffect,repEffect,freq,propHet,encN,
							# 		popSize,noMaleCosts,noFemaleCosts,repN))]
							pool.apply_async(runSim,args=(mutPos,surEffect,repEffect,
								freq,propHet,encN,popSize,noMaleCosts,noFemaleCosts,repN))
						else:
							runSim(mutPos,surEffect,repEffect,freq,propHet,encN,
								popSize,noMaleCosts,noFemaleCosts,repN)

	if run_mp:
		# # Wait for the results
		# for proc in async_objects:
		# 	proc.wait()
		# Close the process pool
		pool.close()
		# block until all tasks are complete and processes close
		pool.join()

	return


# Main function, for managing paramter checks
def main():
	# Parse arguments
	args = parse_args()

	global run_mp
	run_mp = not args.nomp

	global numReps
	numReps = args.rep

	global freq_min
	global freq_max
	global freq_step
	freq_min = args.freq_min
	freq_max = args.freq_max
	freq_step = args.freq_step

	global ph_min
	global ph_max
	global ph_step
	ph_min = args.ph_min
	ph_max = args.ph_max
	ph_step = args.ph_step

	global rep_min
	global rep_max
	global rep_step
	rep_min = args.rep_min
	rep_max = args.rep_max
	rep_step = args.rep_step

	global encNums
	encNums = args.enc_N

	global surEffect
	surEffect = args.sur
	global popSize
	popSize = int(args.N)
	global mutPos
	mutPos = 0.5
	global noMaleCosts
	global noFemaleCosts
	noMaleCosts = args.no_male_costs
	noFemaleCosts = args.no_female_costs

	runSingMutPatGridPool()

	return



if __name__ == "__main__":
	main()

