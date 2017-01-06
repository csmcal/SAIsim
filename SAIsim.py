
# Simulator script for running SAI forward simulations with particular parameters

import OOSAIpop as sim

# Whoops, a possibility with randomly assigned sexes is population death
size = 1000

# given 10^-8 SAmut/base/gen and 
mutRate = 10^-3

# Approximated by 
mutRateInv = 10^-3

# Deviation of normally distributed offset of mutation effects from [x,1-x], with x = uniform[0,1)
mutEffectDiffSD = .2

# Reasonable minimum length given that 1 distance unit on the chromosome = 1 Morgan
minInvLen = .1

# The expected # of crossovers per chromosome arm
# parameter = 1 gives 1 Morgan chromosome arm
recombRate = 1 

# Must be <= pop size
encounterNum = 100

# SD for the normally distributed additive noise to male quality in the female choice
choiceNoiseSD = .5

# The probability that a gene undergoes conversion in a heterozygous gamete,
#  which one is converted is selected at random
conversionRate = 10.0**-4.0

# Inversion record buffer represents how far from inversion edges to examine mutation effects
invRecBuffer = .1


numGens = 50000

pop = sim.simSAIpopulation(size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate, recombRate, encounterNum, choiceNoiseSD, invRecBuffer)
pop.recordNGens(numGens)
