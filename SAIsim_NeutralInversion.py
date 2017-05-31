
# Script created for running with Condor, to take in effect magnitudes and model the change in frequency
#   of multiple variably linked alleles in an otherwise undifferentiable background

import sys
import SAIsim as sim
import numpy as np

invMutPower = float(sys.argv[1])
replicateNum = int(sys.argv[2])
size = 1000
numGens = 20000
# numGens = 100
recordEveryN = 10



# For generating a SAIsim population object
# Relies on some parameter defaults in SAIpop ()
def genInvMutPop(invMutPower,size):
	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = .1
	mutRateInv = np.power(10,invMutPower)
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1
	encounterNum = 100 #USED
	# encounterNum = 20 #USED
	choiceNoiseSD = .5 #USED
	invRecBuffer = .1
	# Generate the population
	pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
		encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = True,\
		willConvert = False)
	return pop

# Run the Simulation
pop = genInvMutPop(invMutPower,size)
pop.recordNthInNGens(recordEveryN,numGens)
# pop.recordNthInNStopWhenFixed(recordEveryN,numGens)

# Record Relevant Output
invMutStr = repr(invMutPower)+'0'*(4-len(repr(invMutPower)))
repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
invFileName = 'iMP'+invMutStr+'n'+repNumStr

pop.writeInvCharTable(invFileName+'Char.txt')
pop.writeAllInvFreqTable(invFileName+'Count.txt')
# pop.printGenomes()


