
# Script created for running with Condor, to take in effects of two alleles
#   and model the change in frequency in an otherwise undifferentiable background

# May be tested as:
#  python SAIsim_TwoSpacing.py 0.525 0.86 0.11 0.2125 0.7 0.3 1 0.2125 0.55 0.5 100 500 4000 0 0 0
#  time python SAIsim_TwoSpacing.py 0.5375 0.86 0.12 0.225 0.7 0.3 1 0.21 0.55 0.5 100 1000 1000 0 0 1

import sys
import SAIsim as sim

# Expects something [0,0.9] in Morgans
posS = float(sys.argv[1])
surEffectS = float(sys.argv[2])
repEffectS = float(sys.argv[3])
# Expects something [0,0.9] in Morgans
posB = float(sys.argv[4])
surEffectB = float(sys.argv[5])
repEffectB = float(sys.argv[6])

hasInv = bool(int(sys.argv[7]))
breakpoint1 = float(sys.argv[8])
breakpoint2 = float(sys.argv[9])

initFreq = float(sys.argv[10])

encounterNum = int(sys.argv[11])
popSize = int(sys.argv[12])
numGens = int(sys.argv[13])
noMaleCosts = bool(int(sys.argv[14]))
noFemaleCosts = bool(int(sys.argv[15]))
replicateNum = int(sys.argv[16])

startBNearEquil = True

recordEveryN = popSize

# mutPos1 = 0.2125
# mutPos2 = 0.05 + spacing
# size = 1000
# numGens = 20*size
# recordEveryN = 10
# lowInvStart = False


# For generating a SAIsim population object from a two mutation specification
#   and potentially an inversion
# Population starts with all at a specified initFreq frequency
# Relies on some parameter defaults in SAIpop
def genPopTwoMutInvHap(posS,surEffectS,repEffectS,posB,surEffectB,repEffectB,
		hasInv,breakpoint1,breakpoint2,initFreq,encounterNum,popSize,
		noMaleCosts,noFemaleCosts):
	# Initialize the mutation for the mutation list for genGenomesSexes
	mutationS = [posS,surEffectS,repEffectS,0]
	mutationB = [posB,surEffectB,repEffectB,1]
	mutList = [[posS,surEffectS,repEffectS,0],[posB,surEffectB,repEffectB,0]]

	# Initialize the inversion for the inversion list for genGenomesSexes
	inversion = [breakpoint1,breakpoint2,0]
	invList = []
	if hasInv:
		invList = [inversion]

	# Initialize the heplotype definitions, to be passed to the population object
	doubleHap = [[[mutationS,mutationB],[]]]
	if posS > posB:
		doubleHap = [[[mutationB,mutationS],[]]]
	sHap = [[[mutationS],[]]]
	bHap = [[[mutationB],[]]]

	if hasInv:
		doubleInvHap = [[[mutationS,mutationB],[inversion]]]
		if posS > posB:
			doubleInvHap = [[[mutationB,mutationS],[inversion]]]
		sInvHap = [[[mutationS],[inversion]]]
		bInvHap = [[[mutationB],[inversion]]]
		emptyInvHap = [[[],[inversion]]]

		num3var = int(initFreq**3*2*popSize)
		num2var = int(initFreq**2*(1-initFreq)*2*popSize)
		num1var = int(initFreq*((1-initFreq)**2)*2*popSize)

		hapList = [[doubleInvHap,num3var],
			[doubleHap,num2var],
			[sInvHap,num2var],
			[sHap,num1var],
			[bInvHap,num2var],
			[bHap,num1var],
			[emptyInvHap,num1var]]

		hapDef = [[[0,1],[0]]]

		if startBNearEquil:
			numInvAndS = int(initFreq**2*.5*2*popSize)
			numInvOrS = int(initFreq*(1-initFreq)*.5*2*popSize)
			numBOnly = int(0.5*((1-initFreq)**2)*2*popSize)

			hapList = [[doubleInvHap,numInvAndS],
				[doubleHap,numInvOrS],
				[sInvHap,numInvAndS],
				[sHap,numInvOrS],
				[bInvHap,numInvOrS],
				[bHap,numBOnly],
				[emptyInvHap,numInvOrS]]

	else:
		num2var = int(initFreq**2*2*popSize)
		num1var = int(initFreq*(1-initFreq)*2*popSize)
		hapList = [[doubleHap,num2var],
			[sHap,num1var],
			[bHap,num1var]]

		hapDef = [[[0,1],[]]]


		if startBNearEquil:
			numS = int(initFreq*.5*2*popSize)
			numBOnly = int(0.5*(1-initFreq)*2*popSize)

			hapList = [[doubleHap,numS],
				[sHap,numS],
				[bHap,numBOnly]]

	print(hapDef)
	print(hapList)

	# Parameters as described in general de novo simulator script, most parameters will not be used
	mutRate = .1
	mutRateInv = .1
	mutEffectDiffSD = .2
	minInvLen = .1
	conversionRate = 10.0**-2.0
	recombRate = 1
	choiceNoiseSD = 1 #USED
	invRecBuffer = .1

	# Generate the population
	(genomes,sexes,record) = sim.genGenoSexFromWholeGenHap(popSize,mutList,invList,hapList)
	pop = sim.SAIpop(popSize, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, 
		conversionRate, recombRate, encounterNum, choiceNoiseSD, invRecBuffer, 
		willMutate=False, willMutInv=False, willConvert=False, 
		noMaleCost=noMaleCosts, noFemaleCost=noFemaleCosts, 
		removeFixed=True, recordHaplotypes=hapDef,
		genomes=genomes, sexes=sexes, record=record)
	return pop


# Run the Simulation
pop = genPopTwoMutInvHap(posS,surEffectS,repEffectS,posB,surEffectB,repEffectB,
		hasInv,breakpoint1,breakpoint2,initFreq,encounterNum,popSize,
		noMaleCosts,noFemaleCosts)
pop.recordNthInNStopWhenFixed(recordEveryN,numGens)

endStr = 'Gen '+str(pop.record[0][-1])+'\n'
endAlleleSize = 2*pop.record[5][-1] # In case the last record is an adult stage
if hasInv:
	if -1 == pop.record[3][0][4]:
		endStr += 'Inv Count\t'+str(pop.record[4][0][0][-1])+'\tFreq\t'+str(pop.record[4][0][0][-1]/endAlleleSize)+'\n'
	else: # must have been lost
		endStr += 'Inv Unrecorded\n'
if -1 == pop.record[1][0][5]:
	endStr += 'Mut S Count\t'+str(pop.record[2][0][0][-1])+'\tFreq\t'+str(pop.record[2][0][0][-1]/endAlleleSize)+'\n'
else: # must have been lost
	endStr += 'Mut S Unrecorded\n'
if -1 == pop.record[1][1][5]:
	endStr += 'Mut B Count\t'+str(pop.record[2][1][0][-1])+'\tFreq\t'+str(pop.record[2][1][0][-1]/endAlleleSize)+'\n'
else: # must have been lost
	endStr += 'Mut B Unrecorded'
print(endStr)

# Record Relevant Output
# ASSUMES MUTATION ID = 0, effect values only to two decimal places, =< four digit replicate ID
surEffStrS = repr(surEffectS)+'0'*(5-len(repr(surEffectS)))
repEffStrS = repr(repEffectS)+'0'*(5-len(repr(repEffectS)))
surEffStrB = repr(surEffectB)+'0'*(5-len(repr(surEffectB)))
repEffStrB = repr(repEffectB)+'0'*(5-len(repr(repEffectB)))
posStrS = '{:.4f}'.format(posS)
posStrB = '{:.4f}'.format(posB)
hasInvStr = str(int(hasInv))
initFreqStr = '{:.3f}'.format(initFreq)
repNumStr = '{:03d}'.format(replicateNum)


filePrefix = 'pos'+posStrS+'inv'+hasInvStr+'freq'+initFreqStr+'n'+repNumStr

# Write the output files
pop.writeAllMutFreqTable(filePrefix+'Mut.txt')
pop.writeInversion(filePrefix+'Inv.txt',0,fillRecord=True) # Should silently return false if Inv not present
pop.writeHaplotypes(filePrefix+'Hap.txt',0,fillRecord=True)


