
# Script created for running with Condor,
# 	to take in frequency and reproductive quality values and model the relative mating success,
# 	of a single mutation in an otherwise undifferentiable background for a single mate competition
# Meant to be run with a number of simulation variations:
# 	variable encounter number
# 	variable population size
# 	variable departure from Hard-Weinberg Equilibrium

import sys
import SAIsim as sim

# Main function, for managing paramter checks and writing the dag
def main():
	# Parse arguments
	
	# Assumes arguments = "$(surEffect) $(repEffect) $(frequency) $(propHet) 
	# 		$(noMaleCosts) $(noFemaleCosts) $(encounterNum) $(popSize) $(repNum)"
	surEffect = float(sys.argv[1])
	repEffect = float(sys.argv[2])
	frequency = float(sys.argv[3])
	propHet = float(sys.argv[4])
	noMaleCosts = bool(int(sys.argv[5]))
	noFemaleCosts = bool(int(sys.argv[6]))
	encounterNum = int(sys.argv[7])
	popSize = int(sys.argv[8])
	replicateNum = int(sys.argv[9])
	mutPos = 0.5

	
	runSim(mutPos,surEffect,repEffect,frequency,propHet,encounterNum,popSize,
		noMaleCosts,noFemaleCosts,replicateNum)
	return


# For generating a SAIsim population object from a single mutation specification
# Population starts with the mutation in the middle of the only chromosome arm
# Relies on some parameter defaults in SAIpop
def genPopSingMutEffects(mutPos,surEffect,repEffect,frequency,propHet,encounterNum,popSize,
		noMaleCosts,noFemaleCosts):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#   size # of mutations as the diploid population has 2N alleles
	mutation = [mutPos,surEffect,repEffect,0,int(2*popSize*frequency)]
	# Parameters as described in general de novo simulator script,
	#  most parameters will not be used or take effect
	mutRate = .1
	mutRateInv = .1
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1
	choiceNoiseSD = 1 #USED
	invRecBuffer = .1
	# Generate the population
	(genomes,sexes,record) = sim.genGenomesSexes(popSize,[mutation],[])
	pop = sim.SAIpop(popSize, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, 
		conversionRate, recombRate, encounterNum, choiceNoiseSD, invRecBuffer, 
		willMutate = False, willMutInv = False, willConvert = False, willRecombine = False, 
		noMaleCost = noMaleCosts, noFemaleCost = noFemaleCosts, removeFixed=True, #verbose=True,
		genomes=genomes, sexes=sexes, record=record)
	return pop

# Run the Simulation
def runSim(mutPos,surEffect,repEffect,frequency,propHet,encounterNum,popSize,
		noMaleCosts,noFemaleCosts,replicateNum,verbose=False):
	# print((mutPos,surEffect,repEffect,frequency,propHet,encounterNum,popSize,
	# 	noMaleCosts,noFemaleCosts,replicateNum,verbose))

	pop = genPopSingMutEffects(mutPos,surEffect,repEffect,frequency,propHet,
		encounterNum,popSize,noMaleCosts,noFemaleCosts)
	pop.step(recReprodParents=True)

	# Record Relevant Output
	# ASSUMES MUTATION ID = 0
	surEffStr = repr(surEffect)+'0'*(5-len(repr(surEffect)))
	repEffStr = repr(repEffect)+'0'*(5-len(repr(repEffect)))
	freqStr = '{:.3f}'.format(frequency)
	propHetStr = '{:.3f}'.format(propHet)
	encNumStr = '{:03d}'.format(encounterNum)
	popSizeStr = '{:03d}'.format(popSize)
	repNumStr = '{:03d}'.format(replicateNum)
	mutFileName = 's'+surEffStr+'r'+repEffStr+'f'+freqStr+'d'+propHetStr+\
		'nmc'+repr(int(noMaleCosts))+'nfc'+repr(int(noFemaleCosts))+\
		'm'+encNumStr+'N'+popSizeStr+'n'+repNumStr+'.txt'
	pop.writeMutation(mutFileName,0,fillRecord = True)
	# print(mutFileName)
	if verbose:
		pop.printRecord()
	return



if __name__ == "__main__":
	main()



