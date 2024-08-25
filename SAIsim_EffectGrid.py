
# Script created for running with Condor,
# 	to take in sur, rep values and model the change in frequency of a single mutation
# 	in an otherwise undifferentiable background
# Meant to be run with a number of simulation variations:
# 	no male survival costs
# 	no female survival costs
# 	variable encounter number
# 	variable population size

import sys
import SAIsim as sim

# Assumes arguments = "$(survEffect) $reprEffect) $(encounterNum) $(popSize) $(numGens) $(noMaleCosts) $(noFemaleCosts) $(repNum)"
surEffect = float(sys.argv[1])
repEffect = float(sys.argv[2])
encounterNum = int(sys.argv[3])
popSize = int(sys.argv[4])
numGens = int(sys.argv[5])
noMaleCosts = bool(int(sys.argv[6]))
noFemaleCosts = bool(int(sys.argv[7]))
replicateNum = int(sys.argv[8])
mutPos = 0.5
recordEveryN = popSize


# For generating a SAIsim population object from a single mutation specification
# Population starts with the mutation at 50% frequency in the middle of the only chromosome arm
# Relies on some parameter defaults in SAIpop
def genPopSingMutEffects(mutPos,surEffect,repEffect,encounterNum,popSize,noMaleCosts,noFemaleCosts):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#   size # of mutations as the diploid population has 2N alleles
	mutation = [mutPos,surEffect,repEffect,0,popSize]
	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = .1
	mutRateInv = .1
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1
	# encounterNum = 100 #USED
	choiceNoiseSD = 1 #USED
	invRecBuffer = .1
	# Generate the population
	(genomes,sexes,record) = sim.genGenomesSexes(popSize,[mutation],[])
	pop = sim.SAIpop(popSize, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,
		encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = False,
		willConvert = False, willRecombine = False, noMaleCost = noMaleCosts, noFemaleCost = noFemaleCosts,
		removeFixed=True, #verbose=True,
		genomes=genomes, sexes=sexes, record=record)
	return pop

# Run the Simulation
pop = genPopSingMutEffects(mutPos,surEffect,repEffect,encounterNum,popSize,noMaleCosts,noFemaleCosts)
# pop.recordNthInNGens(recordEveryN,numGens)
pop.recordNthInNStopWhenFixed(recordEveryN,numGens)

# Record Relevant Output
# ASSUMES MUTATION ID = 0, effect values only to two decimal places, =< three digit replicate ID
surEffStr = repr(surEffect)+'0'*(5-len(repr(surEffect)))
repEffStr = repr(repEffect)+'0'*(5-len(repr(repEffect)))
# repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
repNumStr = '{:03d}'.format(replicateNum)
mutFileName = 's'+surEffStr+'r'+repEffStr+'n'+repNumStr+'.txt'
pop.writeMutation(mutFileName,0,fillRecord = True)
# pop.printRecord()



