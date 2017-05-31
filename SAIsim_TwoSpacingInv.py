
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
# numGens = 100
# numGens = 100
recordEveryN = 10


# For generating a SAIsim population object from a two mutation specification
# Population starts with both at 50% frequency in the middle of separate chromosome arms
# Relies on some parameter defaults in SAIpop
def genPopTwoMutInvHap(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size):
	# Initialize the mutation for the mutation list for genGenomesSexes,
	#   size # of mutations as the diploid population has 2N alleles
	mutationS = [mutPos1,surEffectS,repEffectS,0]
	mutationB = [mutPos2,surEffectB,repEffectB,1]
	mutList = [[mutPos1,surEffectS,repEffectS,0],[mutPos2,surEffectB,repEffectB,0]]

	inversion = [0.02,0.38,0]
	invList = [inversion]
	# invList = [inversion,inversion,inversion,inversion]
	# inversion2 = [0.02,0.38,1]
	# inversion3 = [0.02,0.38,2]
	# inversion4 = [0.02,0.38,3]

	doubleInvHap = [[[mutationS,mutationB],[inversion]]]
	doubleHap = [[[mutationS,mutationB],[]]]
	sInvHap = [[[mutationS],[inversion]]]
	# sInvHap = [[[mutationS],[inversion2]]]
	sHap = [[[mutationS],[]]]
	bInvHap = [[[mutationB],[inversion]]]
	# bInvHap = [[[mutationB],[inversion3]]]
	bHap = [[[mutationB],[]]]
	emptyInvHap = [[[],[inversion]]]
	# emptyInvHap = [[[],[inversion4]]]
	# emptyHap = [[[],[]]]
	invClassNum = size//50
	# classNum = size//4
	# numPerHap = size//4
	numPerHap = size//8
	hapList = [[doubleInvHap,numPerHap],[doubleHap,numPerHap],[sInvHap,numPerHap],\
		[sHap,numPerHap],[bInvHap,numPerHap],[bHap,numPerHap],[emptyInvHap,numPerHap]]
	hapCounts = [[numPerHap] * 8]
	# hapList = [[doubleInvHap,invClassNum],[doubleHap,numPerHap],[sInvHap,invClassNum],\
	# 	[sHap,numPerHap],[bInvHap,invClassNum],[bHap,numPerHap],[emptyInvHap,invClassNum]]
	# print(hapList)

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
	return pop, hapCounts
	# return pop

def getHapCounts(pop):
	hapCounts = [0]*8
	for indiv in pop.males + pop.females:
		for chrom in indiv.genome:
				for hom in chrom:
					numInv = len(hom[1])
					numMut = len(hom[0])
					# print((numInv,numMut))
					if numInv == 1:
						if numMut == 2:
							hapCounts[0] += 1
						elif numMut == 0:
							hapCounts[3] += 1
						elif numMut == 1:
							if hom[0][0][3] == 0:
								hapCounts[1] += 1
							elif hom[0][0][3] == 1:
								hapCounts[2] += 1
							else:
								print("Unexpected mutation: "+str(hom[0][0]))
								print(hom)
						else:
							print("Too many mutations in Hap accounting")
							print(hom)
					elif numInv == 0:
						if numMut == 2:
							hapCounts[4] += 1
						elif numMut == 0:
							hapCounts[7] += 1
						elif numMut == 1:
							if hom[0][0][3] == 0:
								hapCounts[5] += 1
							elif hom[0][0][3] == 1:
								hapCounts[6] += 1
							else:
								print("Unexpected mutation: "+str(hom[0][0]))
								print(hom)
						else:
							print("Too many mutations in Hap accounting")
							print(hom)
					else:
						print("Too many inversions in Hap accounting")
						print(hom)
					# for mut in hom[0]:
					# for inv in hom[1]:
	return hapCounts

# For recording haplotype data as well:
def recordNthInNGensWithHap(pop, setSize, numGens, hapCounts):
	# Be wary of this, it may not simulate the full number of generations,
	#   but really you wouldn't necessarily record the final step anyway
	numSets = int(numGens/setSize)
	for i in range(numSets):
		pop.stepNGens(setSize-1)
		pop.recordStep()
		hapCounts += [getHapCounts(pop)]
	return hapCounts


def writeAllHapCounts(counts,pop,filename):
	with open(filename, 'w') as outfile:
		outfile.write('Gen\tDouInv\tSInv\tBInv\tInv\tDouble\tS\tB\tWT\n')
		for g in range(len(counts)):
			outfile.write(str(pop.record[0][g])+"\t"+"\t".join([str(x) for x in counts[g]])+"\n")

def writeFinalHapCounts(counts,filename):
	outfile = open(filename, 'w')
	outfile.write('DouInv\tSInv\tBInv\tInv\tDouble\tS\tB\tWT\n')
	outfile.write("\t".join([str(x) for x in counts]))


# Run the Simulation
pop, hapCounts = genPopTwoMutInvHap(mutPos1,surEffectS,repEffectS,mutPos2,surEffectB,repEffectB,size)
# pop.recordNthInNGens(recordEveryN,numGens)
# pop.recordNthInNStopWhenFixed(recordEveryN,numGens)
hapCounts = recordNthInNGensWithHap(pop, recordEveryN, numGens, hapCounts)

# Record Relevant Output
# ASSUMES MUTATION ID = 0, effect values only to two decimal places, =< four digit replicate ID
surEffStrS = repr(surEffectS)+'0'*(5-len(repr(surEffectS)))
repEffStrS = repr(repEffectS)+'0'*(5-len(repr(repEffectS)))
surEffStrB = repr(surEffectB)+'0'*(5-len(repr(surEffectB)))
repEffStrB = repr(repEffectB)+'0'*(5-len(repr(repEffectB)))
spacing = repr(spacing)+'0'*(5-len(repr(spacing)))
repNumStr = '0'*(3-len(str(replicateNum)))+str(replicateNum)
mutFileName = 'sS'+surEffStrS+'rS'+repEffStrS+'sB'+surEffStrB+'rB'+repEffStrB\
	+'spacing'+spacing+'n'+repNumStr
# mutFileName = 'spacing'+spacing+'n'+repNumStr
# pop.writeMutation(mutFileName+'Small.txt',0)
# pop.writeMutation(mutFileName+'Big.txt',1)
pop.writeAllMutFreqTable(mutFileName+'Mut.txt')
# pop.writeAllInvFreqTable(mutFileName+'Inv.txt')
# for invID in range(4):
# 	filename = mutFileName+'Inv'+str(invID)+'.txt'
# 	pop.writeInversion(filename,invID)
pop.writeInversion(mutFileName+'Inv.txt',0)
# pop.printGenomes()

# writeFinalHapCounts(getHapCounts(pop),mutFileName+'Hap.txt')
writeAllHapCounts(hapCounts,pop,mutFileName+'Hap.txt')

