
# Script created for running with Condor, to take in effects of two alleles
#   and model the change in frequency in an otherwise undifferentiable background

import sys
import SAIsim as sim

surEffectS = float(sys.argv[1])
repEffectS = float(sys.argv[2])
surEffectB = float(sys.argv[3])
repEffectB = float(sys.argv[4])
# Expects something [0,0.9] in Morgans
spacing = float(sys.argv[5])
replicateNum = int(sys.argv[6])
mutPos1 = 0.05
mutPos2 = 0.05 + spacing
size = 1000
numGens = 20000
# numGens = 3000
# size = 1000
# numGens = 100
recordEveryN = 10


# For generating a SAIsim population object from a two mutation specification
# Population starts with both at 50% frequency in the middle of separate chromosome arms
# Relies on some parameter defaults in SAIpop
def genPopTwoMutEffects(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#   size # of mutations as the diploid population has 2N alleles
	mutationS = [mutPos1,surEffectS,repEffectS,0,size]
	mutationB = [mutPos2,surEffectB,repEffectB,0,size]
	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = .1
	mutRateInv = .1
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1
	encounterNum = 100 #USED
	# encounterNum = 20 #USED
	choiceNoiseSD = .5 #USED
	invRecBuffer = .1
	# Generate the population
	(genomes,sexes,record) = sim.genGenomesSexes(size,[mutationS,mutationB],[])
	pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
		encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = False,\
		willConvert = False, genomes=genomes, sexes=sexes, record=record)
	return pop

# Run the Simulation
pop = genPopTwoMutEffects(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size)
# pop.recordNthInNGens(recordEveryN,numGens)
pop.recordNthInNStopWhenFixed(recordEveryN,numGens)

# Record Relevant Output
# ASSUMES MUTATION ID = 0, effect values only to two decimal places, =< four digit replicate ID
surEffStrS = repr(surEffectS)+'0'*(5-len(repr(surEffectS)))
repEffStrS = repr(repEffectS)+'0'*(5-len(repr(repEffectS)))
surEffStrB = repr(surEffectB)+'0'*(5-len(repr(surEffectB)))
repEffStrB = repr(repEffectB)+'0'*(5-len(repr(repEffectB)))
spacing = repr(spacing)+'0'*(5-len(repr(spacing)))
repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
# mutFileName = 'sS'+surEffStrS+'rS'+repEffStrS+'sB'+surEffStrB+'rB'+repEffStrB\
# 	+'spacing'+spacing+'n'+repNumStr
mutFileName = 'spacing'+spacing+'n'+repNumStr
# mutFileName = 'spacing'+spacing+'n'+repNumStr
# pop.writeMutation(mutFileName+'Small.txt',0)
# pop.writeMutation(mutFileName+'Big.txt',1)
pop.writeAllMutFreqTable(mutFileName+'Mut.txt')
# pop.printGenomes()


