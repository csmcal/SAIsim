
# Script created for running with Condor, to take in mutation rates and run the full SAIsim model

import sys
import SAIsim as sim
import numpy as np

mutPower = float(sys.argv[1])
invMutPower = float(sys.argv[2])
convPower = float(sys.argv[3])
replicateNum = int(sys.argv[4])
size = 1000
# numGens = 20000
numGens = 100
recordEveryN = 10



# For generating a SAIsim population object
# Relies on some parameter defaults in SAIpop ()
def genFullPop(mutPower,invMutPower,convPower,size):
	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = np.power(10,mutPower)
	mutRateInv = np.power(10,invMutPower)
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = np.power(10,convPower)
	recombRate = 1
	encounterNum = 100 #USED
	# encounterNum = 20 #USED
	choiceNoiseSD = .5 #USED
	invRecBuffer = .1
	# Generate the population
	pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
		encounterNum, choiceNoiseSD, invRecBuffer)
	# pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
	# 	encounterNum, choiceNoiseSD, invRecBuffer, willConvert = False)
	return pop

# Run the Simulation
pop = genFullPop(mutPower,invMutPower,convPower,size)
pop.recordNthInNGens(recordEveryN,numGens)

# Record Relevant Output
mutStr = repr(mutPower)+'0'*(4-len(repr(mutPower)))
invMutStr = repr(invMutPower)+'0'*(4-len(repr(invMutPower)))
convStr = repr(convPower)+'0'*(4-len(repr(convPower)))
repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
# filePrefix = 'mp'+mutStr+'imp'+invMutStr+'cp'+convStr+'n'+repNumStr
filePrefix = 'n'+repNumStr

pop.writeRecordTables(filePrefix)


