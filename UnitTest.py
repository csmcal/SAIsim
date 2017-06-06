
# Simulator script for testing SAI forward simulations with particular parameters

import SAIsim as sim

# Limited size for testing
size = 10
# given 10^-8 SAmut/base/gen and 
mutRate = .1
# Approximated by 
mutRateInv = .1
# Deviation of normally distributed offset of mutation effects from [x,1-x], with x = uniform[0,1)
mutEffectDiffSD = .2
# Reasonable minimum length given that 1 distance unit on the chromosome = 1 Morgan
minInvLen = .1
# Chromosome length in M, Drosophila arms are 50 cM
lenChrom = 1.0
# The expected # of crossovers per chromosome arm
# parameter = 1 gives 1 Morgan chromosome arm
recombRate = 1 
# Must be <= pop size
encounterNum = 20
# SD for the normally distributed additive noise to male quality in the female choice
choiceNoiseSD = .5
# The probability per instance of heterozygosity that a gene undergoes conversion in a heterozygous gamete,
#  which one is converted is selected at random, treated independently from number of crossover events
conversionRate = 10.0**-2.0
# Inversion record buffer represents how far from inversion edges to examine mutation effects
# CONSIDER making this a parameter of the record function (Doesn't really work as recording is simultaneous)
invRecBuffer = .1

# Testing Recombination
isFly = True
willConvert = False
willRecombine = True

mut1 = [0.05,.7,.2,0]
mut2 = [0.40,.8,.1,1]

inv1 = [0.02,0.98,0]
inv2 = [0.02,0.98,1]
inv3 = [0.04,0.10,3]
inv4 = [0.34,0.45,4]

def testSameInvRec():
	print('Testing recombination between the same inversion:')
	testGenome = [[[[mut1],[inv1]],[[mut2],[inv1]]]]
	print('Test Genome:')
	print(testGenome)
	for i in range(10):
		fly = sim.individual('F',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
			lenChrom,isFly,willConvert,willRecombine,testGenome)
		print(fly.genGamete())

def testSameInvDiffLabelRec():
	print('Testing recombination between differently labeled inversions:')
	testGenome = [[[[mut1],[inv1]],[[mut2],[inv2]]]]
	print('Test Genome:')
	print(testGenome)
	for i in range(10):
		fly = sim.individual('F',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
			lenChrom,isFly,willConvert,willRecombine,testGenome)
		print(fly.genGamete())

# Run the tests
testSameInvRec()
testSameInvDiffLabelRec()



# # For genGenomesSexes
# initMutList = [[0.5,0.95,0.1,0,40]]
# initInvList = [[0.01,0.06,0,60]]
# genomes = [[[[[],[]],[[],[]]]]]
# sexes = ['F']

# # For genGenoSexFromWholeGenHap
# mutList = [[0.5,0.95,0.1,0],[0.04,0.45,0.45,0]]
# invList = [[0.01,0.06,0]]
# mut0 = [0.5,0.95,0.1,0]
# mut1 = [0.04,0.45,0.45,1]
# inv0 = [0.01,0.06,0]
# hapList = [[[[[mut0,mut1],[inv0]]],6],[[[[mut0],[inv0]]],4]]


# numGens = 500

# # Now run the simulator session in length and recording intervals as desired using above parameters

# # (genomes,sexes,record) = sim.genGenomesSexes(size,initMutList,initInvList,genomes=genomes,sexes=sexes)
# (genomes,sexes,record) = sim.genGenoSexFromWholeGenHap(size,mutList,invList,hapList,genomes=genomes,sexes=sexes)

# pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,\
# 	encounterNum, choiceNoiseSD, invRecBuffer, willMutate = False, willMutInv = False,\
# 	noMaleCost = True, genomes=genomes,sexes=sexes,record=record)
# # pop = sim.SAIpop(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate,encounterNum, choiceNoiseSD, invRecBuffer)
# # pop.stepNGens(numGens)
# # pop.recordEveryNGens(10,200)
# pop.recordNGens(numGens)
# # pop.printGenomes()
# # pop.printRecord()
# pop.writeRecordTables('testOutput/Test1')
