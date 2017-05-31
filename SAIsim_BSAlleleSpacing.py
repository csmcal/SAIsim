
# Script created for running with Condor, to take in effect magnitudes and model the change in frequency
#   of multiple variably linked alleles in an otherwise undifferentiable background

import sys
import SAIsim as sim

surEffectS = float(sys.argv[1])
repEffectS = float(sys.argv[2])
numMut = int(sys.argv[3])
spacing = float(sys.argv[4])
surEffectB = float(sys.argv[5])
repEffectB = float(sys.argv[6])
replicateNum = int(sys.argv[7])
size = 1000
numGens = 20000
# numGens = 1000
recordEveryN = 10

# For generating a set of equivalent 
def genMutationList(surEffect,repEffect,numMut,spacing):
	currentSpacing = 0.05
	mutList = []
	mutHap = []
	indivSurEffect = np.power(surEffect,(1/float(numMut)))
	indivRepEffect = repEffect/float(numMut)
	for i in range(numMut):
		mutList += [[currentSpacing,indivSurEffect,indivRepEffect,0]]
		mutHap += [[currentSpacing,indivSurEffect,indivRepEffect,i+1]]
	return (mutList,mutHap)


# For generating a SAIsim population object from a two mutation specification
# One mutation gets broken into spacing and 
# Population starts with both at 50% frequency in the middle of separate chromosome arms
# Relies on some parameter defaults in SAIpop
def genPopTwoChromHap(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#    set the # of mutations = size as the diploid population has 2N alleles
	mutationS = [0.5,surEffectS,repEffectS,1]
	(mutList,mutHap) = genMutationList(surEffect,repEffect,numMut,spacing)
	mutationB = [mutPos2,surEffectB,repEffectB,1]
	mutList = [[mutPos1,surEffectS,repEffectS,0],[mutPos2,surEffectB,repEffectB,0]]

	doubleHap = []

	inversion = [0.02,0.38,0]
	# invList = [inversion]
	invList = [inversion,inversion,inversion,inversion]
	inversion2 = [0.02,0.38,1]
	inversion3 = [0.02,0.38,2]
	inversion4 = [0.02,0.38,3]

	doubleInvHap = [[[mutationS,mutationB],[inversion]]]
	doubleHap = [[[mutationS,mutationB],[]]]
	# sInvHap = [[[mutationS],[inversion]]]
	sInvHap = [[[mutationS],[inversion2]]]
	sHap = [[[mutationS],[]]]
	# bInvHap = [[[mutationB],[inversion]]]
	bInvHap = [[[mutationB],[inversion3]]]
	bHap = [[[mutationB],[]]]
	# emptyInvHap = [[[],[inversion]]]
	emptyInvHap = [[[],[inversion4]]]
	# emptyHap = [[[],[]]]
	invClassNum = size//50
	# classNum = size//4
	numPerHap = size//4
	# hapList = [[doubleInvHap,numPerHap],[doubleHap,numPerHap],[sInvHap,numPerHap],\
	# 	[sHap,numPerHap],[bInvHap,numPerHap],[bHap,numPerHap],[emptyInvHap,numPerHap]]
	hapList = [[doubleInvHap,invClassNum],[doubleHap,numPerHap],[sInvHap,invClassNum],\
		[sHap,numPerHap],[bInvHap,invClassNum],[bHap,numPerHap],[emptyInvHap,invClassNum]]

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
	(genomes,sexes,record) = sim.genGenoSexFromWholeGenHap(size,mutList,invList,hapList)
	pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
		encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = False,\
		willConvert = False, genomes=genomes, sexes=sexes, record=record)
	return pop

# Run the Simulation
pop = genPopTwoMutInvHap(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size)
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
mutFileName = 'sS'+surEffStrS+'rS'+repEffStrS+'sB'+surEffStrB+'rB'+repEffStrB\
	+'spa'+spacing+'n'+repNumStr
# mutFileName = 'spacing'+spacing+'n'+repNumStr
# pop.writeMutation(mutFileName+'Small.txt',0)
# pop.writeMutation(mutFileName+'Big.txt',1)
pop.writeAllMutFreqTable(mutFileName+'Mut.txt')
pop.writeAllInvFreqTable(mutFileName+'Inv.txt')
for invID in range(4):
	filename = mutFileName+'Inv'+str(invID)+'.txt'
	pop.writeInversion(filename,invID)
# pop.printGenomes()


