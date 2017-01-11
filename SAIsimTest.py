
# Simulator script for testing SAI forward simulations with particular parameters

import OOSAIpop as sim

# Whoops, a possibility with randomly assigned sexes is population death
size = 10

# given 10^-8 SAmut/base/gen and 
mutRate = 1

# Approximated by 
mutRateInv = 1

# Deviation of normally distributed offset of mutation effects from [x,1-x], with x = uniform[0,1)
mutEffectDiffSD = .2

# Reasonable minimum length given that 1 distance unit on the chromosome = 1 Morgan
minInvLen = .1

# The expected # of crossovers per chromosome arm
# parameter = 1 gives 1 Morgan chromosome arm
recombRate = 1 

# Must be <= pop size
encounterNum = 4

# SD for the normally distributed additive noise to male quality in the female choice
choiceNoiseSD = .5

# The probability per instance of heterozygosity that a gene undergoes conversion in a heterozygous gamete,
#  which one is converted is selected at random, treated independently from number of crossover events
conversionRate = 10.0**-2.0

# Inversion record buffer represents how far from inversion edges to examine mutation effects
invRecBuffer = .1


numGens = 10

# Now run the simulator session in length and recording intervals as desired using above parameters
pop = sim.simSAIpopulation(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate, encounterNum, choiceNoiseSD, invRecBuffer)
# pop.stepNGens(numGens)
# pop.recordEveryNGens(10,500)
pop.recordNGens(numGens)
pop.printRecord()
# pop.testSharedData()
