
# Script created for running with Condor, to take in effects and model the change in frequency
#   in an otherwise undifferentiable background

import sys
import SAIsim as sim

surEffect = float(sys.argv[1])
repEffect = float(sys.argv[2])
replicateNum = int(sys.argv[3])
mutPos = 0.5
size = 1000
numGens = 20000
# size = 1000
# numGens = 100
recordEveryN = 10


# For generating a SAIsim population object from a single mutation specification
# Population starts with the mutation at 50% frequency in the middle of the only chromosome arm
# Relies on some parameter defaults in SAIpop
def genPopSingMutEffects(mutPos,surEffect,repEffect,size):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#   size # of mutations as the diploid population has 2N alleles
	mutation = [mutPos,surEffect,repEffect,0,size]
	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = .1
	mutRateInv = .1
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1 
	encounterNum = 100 #USED
	choiceNoiseSD = .5 #USED
	invRecBuffer = .1
	# Generate the population
	(genomes,sexes,record) = sim.genGenomesSexes(size,[mutation],[])
	pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
		encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = False,\
		willConvert = False, willRecombine = False, genomes=genomes, sexes=sexes, record=record)
	return pop

# Run the Simulation
pop = genPopSingMutEffects(mutPos,surEffect,repEffect,size)
# pop.recordNthInNGens(recordEveryN,numGens)
pop.recordNthInNStopWhenFixed(recordEveryN,numGens)

# Record Relevant Output
# ASSUMES MUTATION ID = 0, effect values only to two decimal places, =< three digit replicate ID
surEffStr = repr(surEffect)+'0'*(4-len(repr(surEffect)))
repEffStr = repr(repEffect)+'0'*(4-len(repr(repEffect)))
repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
mutFileName = 's'+surEffStr+'r'+repEffStr+'n'+repNumStr+'.txt'
pop.writeMutation(mutFileName,0)



