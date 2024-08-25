
# The Object Oriented implementation of Sexually Antagonistic mutation accumulation
#     in Inversions in a forward Population simulation
# Ideally, models large inversion polymorphisms in populations with high reproductive skew


# REWRITE WITH NUMPY ARRAYS/VECTORS?

# TODO:
#  fix chrom length by collapsing to just recomb. rate or length param. (prob. recomb. rate)
#  add full record of offspring/values for individuals/gen. instead of means/var (also use coeff. of var. for surv.)
#  combine record keeping into the step function to speed up/make simpler (how much benefit?)
#  

import numpy as np
# import pickle as p
import itertools
from collections import Counter


# For counting chromosome proportions to compare to expected proportions
class HashableList(object): 
 
 def __init__(self, val):
    self.val = val 
 
 def __hash__(self):
    return hash(str(self.val)) 
 
 def __repr__(self):
    # define this method to get clean output
    return str(self.val) 
 
 def __eq__(self, other):
    return str(self.val) == str(other.val)


# The object representing an individual in the population under simulation,
# contains methods for generating mutation position and effects, and new recombinant gametes
# recombRate is the expected number of recombination events per chromosome per reproductive event
# meanMutEffect is the expected effect size of a new mutation
# mutEffectDiffSD is the SD of the difference in survival/rep effect size of a new mutation (normal dist)
# genome defaults to 2 copies of a single chormosome,
# ideally: genome = [[[chr1hom1],[chr1hom2]],[[chr2hom1],[chr2hom2]],..]
class individual(object):
	"""individuals represent members of simSAIpopulations"""
	def __init__(self, sex, mutEffectDiffSD, recombRate, conversionRate, minInvLen, lenChrom,
			maleAchiasmy, willConvert, willRecombine, oneInvPerChrom, genome = [[[[],[]],[[],[]]]]):
		# super(individual, self).__init__()
		self.sex = sex
		self.mutEffectDiffSD = mutEffectDiffSD
		self.recombRate = recombRate
		self.conversionRate = conversionRate
		self.minInvLen = minInvLen
		self.lenChrom = lenChrom
		# CONSIDER JUST PASSING VARIABLES WHEN NEEDED
		self.maleAchiasmy = maleAchiasmy
		self.willConvert = willConvert
		self.willRecombine = willRecombine
		self.oneInvPerChrom = oneInvPerChrom
		# genome = [[[chr1hom1],[chr1hom2]],[[chr2hom1],[chr2hom2]],..]
		# where chomosome homologs = [mutist, InvList]
		self.genome = genome


	# For getting the insertion index of a position from a list of mutations 
	#   or inversions, represented individually as [position, ..] lists
	def __getInsInd(self,mutInvList,mutInvPos):
		i = 0
		while (i < len(mutInvList)) and (mutInvList[i][0] <= mutInvPos):
			i += 1
		return i

	# For inserting a mutation or inversion into a corresponding list
	# May rely on passing the list as a reference (doesn't currently)
	def __insert(self,mutInvList,mutInv):
		i = 0
		while (i < len(mutInvList)) and (mutInvList[i][0] <= mutInv[0]):
			i += 1
		mutInvList[i:i] = [mutInv]
		return mutInvList


	# Generates mutation survival and reproductive quality effect values,
	#   from the population class parameters of mean and SD
	# Consider alternate distributions? Lognormal?
	def __genEffectSizes(self):
		base = np.random.ranf()
		offset = np.random.normal(scale=self.mutEffectDiffSD)
		mutEffects = []
		while (min(1-base,base) < offset) or (offset < max(base-1,-base)):
			offset = np.random.normal(scale=self.mutEffectDiffSD)
		# Survival (multiplicative), Quality (additive)
		return [1-base+offset,base+offset]

	# Generates mutation survival and reproductive quality effect values,
	#   with ranges (0, 1] and [0, 1) 
	def __genEffectSizesRand(self):
		return [1-np.random.ranf(),np.random.ranf()]

	# Generates a mutation of the form [position, survival effect, rep effect, ID]
	#   and inserts it into a random position in the genome
	def mutate(self,ID):
		# print(self.genome)
		mutPos = np.random.ranf()*self.lenChrom
		mutEffects = self.__genEffectSizesRand()
		mutation = [mutPos]+mutEffects+[ID]
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Pick a homolog
		homIndex = np.random.randint(0,2)
		# Insert the mutation into the list of mutations for the homolog
		# print(self.genome[chromIndex][homIndex][0])
		self.genome[chromIndex][homIndex][0] = \
			self.__insert(self.genome[chromIndex][homIndex][0],mutation)
		# print(self.genome[chromIndex][homIndex][0])
		# Return the mutation data for recording
		return mutation[0:3]+[chromIndex]

	# Generates data on uninverted segments of the 
	# In order to pick where the inversion mutant may be placed
	def __genOpenUninvRegData(self,chromHomInv):
		openRegLengths = []
		openRegStarts = []
		openRegIndexes = []
		potentialStart = 0
		index = 0
		for inversion in chromHomInv:
			length = inversion[0] - potentialStart
			if length > self.minInvLen:
				openRegLengths += [length]
				openRegStarts += [potentialStart]
				openRegIndexes += [index]
			potentialStart = inversion[1]
			index += 1
		length = self.lenChrom - potentialStart
		# length = 1 - potentialStart
		if length > self.minInvLen:
			openRegLengths += [length]
			openRegStarts += [potentialStart]
			openRegIndexes += [index]
		return (openRegLengths,openRegStarts,openRegIndexes)

	# Generates and inserts an inversion into the genome,
	#   unless there is no open region >= minInvLen
	# Inversions cover the region [pos1,pos2)
	# Does not use resampling 
	# ADD A VERSION THAT ONLY CHECKS ONE RANDOM POSITION AND GIVES UP, OR HANDLING OF OVERLAPS
	def mutateInvOpenReg(self,ID):
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		chromHomInv = self.genome[chromIndex][homIndex][1]
		# return None if one one inversion per chromosome is allowed,
		#   and there is already an inversion present
		if self.oneInvPerChrom and len(chromHomInv) > 0:
			return
		# Pick where the insertion may be placed
		(openRegLengths,openRegStarts,openRegIndexes) = self.__genOpenUninvRegData(chromHomInv)
		# If there is no space for a new inversion >= minInvLen, don't add one, return None
		if len(openRegStarts) == 0:
			return
		# Now sample which open region to add the inversion to
		regChoice = np.random.choice(len(openRegStarts),p=[l/sum(openRegLengths) for l in openRegLengths])
		# Generate the inversion
		posA = openRegStarts[regChoice] + np.random.ranf()*(openRegLengths[regChoice]-self.minInvLen)
		posB = openRegStarts[regChoice] + np.random.ranf()*(openRegLengths[regChoice]-self.minInvLen)
		inversion = [min(posA,posB),max(posA,posB)+self.minInvLen,ID]
		index = openRegIndexes[regChoice]
		chromHomInv[index:index] = [inversion]
		# Return the inversion data to record
		return inversion[0:2]+[chromIndex]

	# Generates a random inversion position and length >= minInvLen,
	#   and inserts an it into the genome,
	#   if that position isn't already occupied in the selected chromosome homolog
	# Inversions cover the region [pos1,pos2)
	# Does not use resampling 
	def mutateInv(self,ID):
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		chromHomInv = self.genome[chromIndex][homIndex][1]
		# Randomly generate an inversion
		posA = np.random.ranf()*(self.lenChrom-self.minInvLen)
		posB = np.random.ranf()*(self.lenChrom-self.minInvLen)
		inversion = [min(posA,posB),max(posA,posB)+self.minInvLen,ID]
		# Check if this inversion may be added, and do so if it can
		if len(chromHomInv) > 0:
			# return None if only one inversion per chromosome is allowed,
			#   and there is already an inversion present 
			if self.oneInvPerChrom:
				return
			else:
				# Check if the inversion may be placed
				invIndex = 0
				while invIndex < len(chromHomInv) and chromHomInv[invIndex][1]<=inversion[0]:
					invIndex += 1
				if invIndex < len(chromHomInv):
					if chromHomInv[invIndex][0]<inversion[1]:
						# Found a conflicting inversion: return None, don't mutate
						return
					# Otherwise, add the inversion and return the inversion itself
					else:
						chromHomInv[invIndex:invIndex] = [inversion]
						return inversion[0:2]+[chromIndex]
						# return inversion
				else:
					chromHomInv += [inversion]
					return inversion[0:2]+[chromIndex]
					# return inversion
		else:
			chromHomInv += [inversion]
			return inversion[0:2]+[chromIndex]
			# return inversion


	# # For updating the repQuality and survival if pre-calculated
	# # Not currently used (obviously)
	# def updatePhenotypes(self):
	# 	return

	# Negative survival effects are independent and multiplicative
	def survival(self):
		# listSurvival = 1.0
		# for chrom in self.genome:
		# 	for mut in chrom[0][0]:
		# 		listSurvival *= mut[1]
		# 	for mut in chrom[1][0]:
		# 		listSurvival *= mut[1]
		muts = []
		for chrom in self.genome:
			muts += chrom[0][0]+chrom[1][0]
		if muts != []:
			# print(np.array(muts))
			survival = np.prod(np.array(muts)[...,1])
			# print(survival)
		else:
			survival = 1
		# print("Survival: "+str(survival)+" by vect; "+str(listSurvival)+" by for loop")
		# assert(survival == listSurvival), "Survival calculation differs between np and list calcs: "+\
		# 	str(survival)+" and "+str(listSurvival)+" respectively"
		return survival

	# FIX THIS - what quality accounting process is biological? (This doesn't seem tooo bad, tbh)
	# ALSO - maybe calculate it once to save time? probably recalculated several times in choosing fathers
	def repQuality(self):
		# listQuality = 0.0
		# for chrom in self.genome:
		# 	for mut in chrom[0][0]:
		# 		listQuality += mut[2]
		# 		# For modeling effects multiplicatively? - requires changing effect generation
		# 		# listQuality *= mut[2]
		# 	for mut in chrom[1][0]:
		# 		listQuality += mut[2]
		# 		# For modeling effects multiplicatively? - requires changing effect generation
		# 		# listQuality *= mut[2]
		muts = []
		for chrom in self.genome:
			muts += chrom[0][0]+chrom[1][0]
		if muts != []:
			# print(np.array(muts))
			quality = np.sum(np.array(muts)[...,2])
			# print(quality)
		else:
			quality = 0
		# print("Rep Quality: "+str(quality)+" by vect; "+str(listQuality)+" by for loop")
		# assert(quality == listQuality), "Reproductive quality calculation differs between np and list calcs: "+\
		# 	str(quality)+" and "+str(listQuality)+" respectively"
		return quality


	def __convCheck(self, chrom, parChrom):
		mutIDs0 = [x[3] for x in chrom[0][0]]
		invIDs0 = [x[2] for x in chrom[0][1]]
		mutIDs1 = [x[3] for x in chrom[1][0]]
		invIDs1 = [x[2] for x in chrom[1][1]]
		# print(mutIDs0+mutIDs1)
		mutCounts = Counter(mutIDs0+mutIDs1)
		# print(mutCounts)
		invCounts = Counter(invIDs0+invIDs1)
		if len(mutCounts) > 0:
			highestMutCount = mutCounts.most_common(1)[0][1]
			if highestMutCount > 2:
				print("Failed Conversion: "+str(highestMutCount)+" copies of mut "+str(mutCounts.most_common(1)[0][0]))
				print("ParHap0: "+str([x[3] for x in parChrom[0][0]]))
				print("ParHap1: "+str([x[3] for x in parChrom[0][1]]))
				print("ConHap0: "+str(mutIDs0))
				print("ConHap1: "+str(mutIDs1))
				print(parChrom)
				print(chrom)
		if len(invCounts) > 0:
			highestInvCount = invCounts.most_common(1)[0][1]
			if highestInvCount > 2:
				print("Failed Conversion from crossover? / too many inversions? Inv data shouldn't be changed")


	# Takes a chromosome and returns the chromosome with converted mutations
	# Setting the boolean willConvert in the population object will turn conversion off in the simulation
	def __convertChrom(self,chrom):
		# parChrom = [[list(chrom[0][0]),list(chrom[0][1])],[list(chrom[1][0]),list(chrom[1][1])]]
		homInd1 = 0
		homInd2 = 0
		# Check for heterozygosity and add/remove mutations following conversionRate probability
		# !!!! Relies on the invariant that earlier indexed mutations have earlier position
		while homInd1 < len(chrom[0][0]) and homInd2 < len(chrom[1][0]):
			mut1 = chrom[0][0][homInd1]
			mut2 = chrom[1][0][homInd2]
			if mut1[0] == mut2[0]:
				homInd1 += 1
				homInd2 += 1
			if mut1[0] < mut2[0]:
				if np.random.ranf() < self.conversionRate:
					# print("Attempting Conversion")
					if np.random.randint(2):
						chrom[0][0][homInd1:homInd1+1] = []
					else:
						chrom[1][0][homInd2:homInd2] = [mut1]
						homInd2 += 1
						homInd1 += 1
				else:
					homInd1 += 1
			if mut1[0] > mut2[0]:
				if np.random.ranf() < self.conversionRate:
					if np.random.randint(2):
						chrom[1][0][homInd2:homInd2+1] = []
					else:
						chrom[0][0][homInd1:homInd1] = [mut2]
						homInd1 += 1
						homInd2 += 1
				else:
					homInd2 += 1
		# Deal with the remaining mutations
		# when one chromosome has no mutations it is dealt with here immediately
		while homInd1 < len(chrom[0][0]):
			if np.random.ranf() < self.conversionRate:
				if np.random.randint(2):
					chrom[0][0][homInd1:homInd1+1] = []
				else:
					chrom[1][0] += [chrom[0][0][homInd1]]
					homInd1 += 1
					homInd2 += 1 # IMPORTANT - so homInd2 is still == len
			else:
				homInd1 += 1
		while homInd2 < len(chrom[1][0]):
			if np.random.ranf() < self.conversionRate:
				if np.random.randint(2):
					chrom[1][0][homInd2:homInd2+1] = []
				else:
					chrom[0][0] += [chrom[1][0][homInd2]]
					homInd2 += 1
			else:
				homInd2 += 1
		# self.__convCheck(chrom, parChrom)
		return chrom


	# Takes two inversion lists for homologous chromosomes,
	# indexes of the closest inversions to start <= the position of interest,
	# and the position of interest
	# Returns the position following the position of interest such that an odd # of crossovers
	# between the two would generate aneuploid gametes, returns -1 if it isn't in such a region
	def __getAllAneupRegion(self,invHom1,invHom2):
		unbalancedRegions = []
		ind1 = 0
		ind2 = 0
		startForward = False # Check unbalanced region starts inside one of the inversions
		unbalStart = 0 # The start of the unbalanced region ^ (when )
		retainFirst = True # For knowing which inversion it is inside of
		while (ind1 < len(invHom1)) and (ind2 < len(invHom2)):
			(start1,end1,ID1) = invHom1[ind1]
			(start2,end2,ID2) = invHom2[ind2]
			if ID1 == ID2: # CHECK POSITIONS DIRECTLY?
				ind1 += 1
				ind2 += 1
			elif startForward:
				if retainFirst:
					if start2 < end1:
						if end2 < end1:
							unbalancedRegions += [(unbalStart,start2),(start2,end2)]
							ind2 += 1
							unbalStart = end2
						if end2 == end1:
							unbalancedRegions += [(unbalStart,start2),(start2,end2)]
							ind1 += 1
							ind2 += 1
							startForward = False
						else:
							unbalancedRegions += [(unbalStart,start2),(start2,end1)]
							ind1 += 1
					else:
						unbalancedRegions += [(unbalStart,end1)]
						ind1 += 1

			elif start1 <= start2:
				if end1 < start2:
					unbalancedRegions += [(start1,end1)]
					ind1 += 1
			elif start2 <= start1:
				if end2 < start1:
					unbalancedRegions += [(start2,end2)]
					ind2 += 1
				elif end2 < end1:
					unbalancedRegions += []

		# Deal with no inversions present on the chromosome
		if len(invHom1) == 0:
			if (len(invHom2) == 0) or (recPos >= invHom2[ind2][1]):
				return unbalancedRegions
			else:
				return invHom2[ind2][1]
		elif len(invHom2) == 0:
			if recPos >= invHom1[ind1][1]:
				return -1
			else:
				return invHom1[ind1][1]
		else:
			# print(invHom1[ind1][:2])
			# print(invHom2[ind2][:2])
			# print(invHom1[ind1] == invHom2[ind2])
			if invHom1[ind1][:2] == invHom2[ind2][:2]:
				return -1
			end1 = invHom1[ind1][1]
			end2 = invHom2[ind2][1]
			# Currently, inversions are counted as [1,2), so crossovers are inside or not following such
			if recPos < end1:
				if recPos < end2:
					return min(end1,end2)
				if (ind2 < len(invHom2)-1) and invHom2[ind2+1][0] < end1:
					return invHom2[ind2+1][0]
				return end1
			elif recPos < end2:
				if (ind1 < len(invHom1)-1) and invHom1[ind1+1][0] < end2:
					return invHom1[ind1+1][0]
				return end2
			else:
				return -1

	# # Takes two inversion lists for homologous chromosomes,
	# # indexes of the closest inversions to start <= the position of interest,
	# # and the position of interest
	# # Returns the position following the position of interest such that an odd # of crossovers
	# # between the two would generate aneuploid gametes, returns -1 if it isn't in such a region
	# def __getAllAneupRegion(self,invHom1,invHom2):
	# 	unbalancedRegions = []
	# 	ind1 = 0
	# 	ind2 = 0
	# 	breaks1 = []
	# 	breaks2 = []
	# 	(start1,end1,ID1) = invHom1[ind1]
	# 	(start2,end2,ID2) = invHom2[ind2]
	# 	while (ind1 < len(invHom1)) and (ind2 < len(invHom2)):
	# 		if ID1 == ID2: # CHECK POSITIONS DIRECTLY?
	# 			ind1 += 1
	# 			ind2 += 1
	# 		while end1 < start2
	# 		elif start1 <= start2:
	# 			if end1 < start2:
	# 				unbalancedRegions += [(start1,end1)]
	# 				ind1 += 1
	# 		elif start2 <= start1:
	# 			if end2 < start1:
	# 				unbalancedRegions += [(start2,end2)]
	# 				ind2 += 1
	# 			elif end2 < end1:
	# 				unbalancedRegions += []

	# 	isEnd1 = False
	# 	isEnd2 = False
	# 	while (ind1 < len(invHom1)) and (ind2 < len(invHom2)):
	# 		(start1,end1,ID1) = invHom1[ind1]
	# 		(start2,end2,ID2) = invHom2[ind2]
	# 		if ID1 == ID2: # CHECK POSITIONS DIRECTLY?
	# 			ind1 += 1
	# 			ind2 += 1
	# 		elif start1 <= start2:
	# 			if end1 < start2:
	# 				unbalancedRegions += [(start1,end1)]
	# 				ind1 += 1
	# 		elif start2 <= start1:
	# 			if end2 < start1:
	# 				unbalancedRegions += [(start2,end2)]
	# 				ind2 += 1
	# 			elif end2 < end1:
	# 				unbalancedRegions += []

	# 	# Deal with no inversions present on the chromosome
	# 	if len(invHom1) == 0:
	# 		if (len(invHom2) == 0) or (recPos >= invHom2[ind2][1]):
	# 			return unbalancedRegions
	# 		else:
	# 			return invHom2[ind2][1]
	# 	elif len(invHom2) == 0:
	# 		if recPos >= invHom1[ind1][1]:
	# 			return -1
	# 		else:
	# 			return invHom1[ind1][1]
	# 	else:
	# 		# print(invHom1[ind1][:2])
	# 		# print(invHom2[ind2][:2])
	# 		# print(invHom1[ind1] == invHom2[ind2])
	# 		if invHom1[ind1][:2] == invHom2[ind2][:2]:
	# 			return -1
	# 		end1 = invHom1[ind1][1]
	# 		end2 = invHom2[ind2][1]
	# 		# Currently, inversions are counted as [1,2), so crossovers are inside or not following such
	# 		if recPos < end1:
	# 			if recPos < end2:
	# 				return min(end1,end2)
	# 			if (ind2 < len(invHom2)-1) and invHom2[ind2+1][0] < end1:
	# 				return invHom2[ind2+1][0]
	# 			return end1
	# 		elif recPos < end2:
	# 			if (ind1 < len(invHom1)-1) and invHom1[ind1+1][0] < end2:
	# 				return invHom1[ind1+1][0]
	# 			return end2
	# 		else:
	# 			return -1

	# For returning recombination points that are balanced, with correct resampling
	def __getBalancedRecomb(self):
		return

	# For efficiently resampling recombination breakpoints
	# when they would otherwise generate unbalanced gametes
	def genGameteResampleCrossovers(self):
		gamGenome = []
		for parChrom in self.genome:
			(numRecomb, recombPositions) = self.__getBalancedRecomb()
			# print('Recombinations: '+str(numRecomb))
			currHom1 = bool(np.random.randint(2))
			# This first copy step is necessary independent of conversion/recombination activity, 
			# to prevent changes to shared-referent chromosomes
			chrom = [[list(parChrom[0][0]),list(parChrom[0][1])],[list(parChrom[1][0]),list(parChrom[1][1])]]
			# Perform fly-specific non-recombination non-conversion check for males
			if self.willConvert and (self.sex == 'F' or (not self.maleAchiasmy)):
				chrom = self.__convertChrom(chrom)
			if (not self.willRecombine) or numRecomb == 0 or (self.maleAchiasmy and self.sex == 'M'): # WARNING - poorly defined numRecomb
				if currHom1:
					gamGenome += [chrom[0]]
				else:
					gamGenome += [chrom[1]]
			else:
				# NEED to allow for the span of positions in the chromosome if given positions outside [0,1)
				recombPositions = np.random.ranf(numRecomb)*self.lenChrom
				# recombPositions = np.random.ranf(numRecomb)
				recombPositions.sort()
				# print('Positions: '+str(recombPositions))
				
				# Construct a recombinant gamete from the recombination positions
				gamMutations = []
				gamInversions = []

				
				recIndex = 0

				invInd1 = 0
				invInd2 = 0
				invHom1 = chrom[0][1]
				invHom2 = chrom[1][1]
				mutInd1 = 0
				mutInd2 = 0

				prevInvInd1 = 0
				prevInvInd2 = 0
				prevMutInd1 = 0
				prevMutInd2 = 0

				while recIndex < numRecomb:
					recPos = recombPositions[recIndex]
					while (mutInd1 < len(chrom[0][0])) and (chrom[0][0][mutInd1][0] <= recPos):
						mutInd1 += 1
					while (mutInd2 < len(chrom[1][0])) and (chrom[1][0][mutInd2][0] <= recPos):
						mutInd2 += 1
					while (invInd1 < len(chrom[0][1])) and (chrom[0][1][invInd1][0] <= recPos):
						invInd1 += 1
					while (invInd2 < len(chrom[1][1])) and (chrom[1][1][invInd2][0] <= recPos):
						invInd2 += 1
					# Add the mutations and inversions to the gamete chromosome
					if currHom1:
						gamMutations += chrom[0][0][prevMutInd1:mutInd1]
						gamInversions += chrom[0][1][prevInvInd1:invInd1]
					else:
						gamMutations += chrom[1][0][prevMutInd2:mutInd2]
						gamInversions += chrom[1][1][prevInvInd2:invInd2]
					# Update the mutation/inversion indexes of the previous breakpoint
					prevInvInd1 = invInd1
					prevInvInd2 = invInd2
					prevMutInd1 = mutInd1
					prevMutInd2 = mutInd2
					# Deal with crossovers in regions potentially generating aneuploidy
					aneupRegEnd = self.__getAneupRegion(invHom1,invHom2,invInd1-1,invInd2-1,recPos)
					if aneupRegEnd == -1:
						currHom1 = not currHom1
						recIndex += 1
					else:
						otherRecInReg = []
						recIndex += 1
						while (recIndex < numRecomb) and (recombPositions[recIndex] < aneupRegEnd):
							otherRecInReg += [recombPositions[recIndex]]
							recIndex += 1
						# If there's an even number of recombinations in the region in total:
						if len(otherRecInReg)%2:
							currHom1 = not currHom1
							for recPos in otherRecInReg:
								while (mutInd1 < len(chrom[0][0])) and (chrom[0][0][mutInd1][0] <= recPos):
									mutInd1 += 1
								while (mutInd2 < len(chrom[1][0])) and (chrom[1][0][mutInd2][0] <= recPos):
									mutInd2 += 1
								# Add the mutations to the gamete chromosome (no inversion change in region)
								if currHom1:
									gamMutations += chrom[0][0][prevMutInd1:mutInd1]
								else:
									gamMutations += chrom[1][0][prevMutInd2:mutInd2]
								# Update the mutation indexes of the previous breakpoint
								prevMutInd1 = mutInd1
								prevMutInd2 = mutInd2
								# Switch the current chromosome from the crossover event
								currHom1 = not currHom1
				# Add the last mutations and inversions to the gamete chromosome
				if currHom1:
					gamMutations += chrom[0][0][prevMutInd1:len(chrom[0][0])]
					gamInversions += chrom[0][1][prevInvInd1:len(chrom[0][1])]
				else:
					gamMutations += chrom[1][0][prevMutInd2:len(chrom[1][0])]
					gamInversions += chrom[1][1][prevInvInd2:len(chrom[1][1])]
				gamGenome += [[gamMutations,gamInversions]]
		return gamGenome


	# Takes two inversion lists for homologous chromosomes,
	# indexes of the closest inversions to start > the position of interest,
	# and the position of interest
	# Returns the position following the position of interest such that an odd # of crossovers
	# between the two would generate aneuploid gametes, returns -1 if it isn't in such a region
	def __getAneupRegion(self,invHom1,invHom2,nextInd1,nextInd2,recPos):
		# switch to the closest inversions to start <= the position of interest
		precIndex1 = nextInd1 - 1
		precIndex2 = nextInd2 - 1
		# Deal with non-aneuploid region cases, overlapping inversion on only one chromosome:
		#   no inversions present on the chromosome, (should be covered by indexes == -1)
		#   inversions end or aren't encountered before the crossover
		if precIndex1 == -1 or recPos >= invHom1[precIndex1][1]:
			if precIndex2 == -1 or recPos >= invHom2[precIndex2][1]:
				return -1
			else:
				end = invHom2[precIndex2][1]
				if nextInd1 < len(invHom1):
					end = min(end, invHom1[nextInd1][0])
				return end
		elif precIndex2 == -1 or recPos >= invHom2[precIndex2][1]:
			# if recPos >= invHom1[precIndex1][1]:
			# 	return -1
			# else:
			# 	return invHom1[precIndex1][1]
			end = invHom1[precIndex1][1]
			if nextInd2 < len(invHom2):
				end = min(end, invHom2[nextInd2][0])
			return end
		elif invHom1[precIndex1][:2] == invHom2[precIndex2][:2]:
			return -1
		else:
			# Expect to KNOW here that recPos < both inversion ends
			assert(recPos < invHom1[precIndex1][1] and recPos < invHom2[precIndex2][1]), "Conditional tree should guarantee this"
			# print(invHom1[precIndex1][:2])
			# print(invHom2[precIndex2][:2])
			# print(invHom1[precIndex1] == invHom2[precIndex2])
			# end1 = invHom1[precIndex1][1]
			# end2 = invHom2[precIndex2][1]
			return min(invHom1[precIndex1][1],invHom2[precIndex2][1])

			# 	if (precIndex2 < len(invHom2)-1) and invHom2[precIndex2+1][0] < end1:
			# 		return invHom2[precIndex2+1][0]
			# 	return end1
			# elif recPos < end2:
			# 	if (precIndex1 < len(invHom1)-1) and invHom1[precIndex1+1][0] < end2:
			# 		return invHom1[precIndex1+1][0]
			# 	return end2
			# else:
			# 	return -1

	# Same as __getAneupRegion, but considers regions between inversions of opposite strands as aneuploidy generating,
	#   to remove recombinants that have multiple inversions
	# ASSUMES only 1 inv/chrom, don't need to account for others on each homolog
	def __getAneupRegionOneInvPerChrom(self,invHom1,invHom2,nextInd1,nextInd2,recPos):
		if len(invHom1) > 1 or len(invHom2) > 1:
			raise self.SimulationStateError((invHom1,invHom2),"When oneInvPerChrom is set, "+\
				"the simulation should never have chromosomes with more than one inversion.")
		# switch to the closest inversions to start <= the position of interest
		precIndex1 = nextInd1 - 1
		precIndex2 = nextInd2 - 1
		if precIndex1 == -1:
			if precIndex2 == -1:
				return -1
			elif recPos >= invHom2[precIndex2][1]:
				if nextInd1 < len(invHom1):
					return invHom1[nextInd1][0]
				else:
					return -1
			else:
				end = invHom2[precIndex2][1]
				if nextInd1 < len(invHom1):
					end = min(end, invHom1[nextInd1][0])
				return end
		elif recPos >= invHom1[precIndex1][1]:
			if precIndex2 == -1:
				if nextInd2 < len(invHom2):
					return invHom2[nextInd2][0]
				else:
					return -1
			elif recPos >= invHom2[precIndex2][1]:
				# ASSUMES only 1 inv/chrom, don't need to account for others on either
				return -1
			else:
				# ASSUMES only 1 inv/chrom, don't need to account for others on hom1
				return invHom2[precIndex2][1]
		else:
			if precIndex2 == -1:
				end = invHom1[precIndex1][1]
				if nextInd2 < len(invHom2):
					end = min(end, invHom2[nextInd2][0])
				return end
			elif recPos >= invHom2[precIndex2][1]:
				# ASSUMES only 1 inv/chrom, don't need to account for others on hom2
				return invHom1[precIndex1][1]
			else:
				# ASSUMES only 1 inv/chrom, don't need to account for others on either
				if invHom1[precIndex1][:2] == invHom2[precIndex2][:2]:
					return -1
				else:
					return min(invHom1[precIndex1][1],invHom2[precIndex2][1])

	# Model independently per chromosome, currently doesn't redraw if recombination events fail,
	# just throws out all recombination events in an aneuploidy-causing region if they are odd in number
	def genGameteDropUnbCrossovers(self):
		gamGenome = []
		for parChrom in self.genome:
			numRecomb = np.random.poisson(self.recombRate)
			# print('Recombinations: '+str(numRecomb))
			currHom1 = bool(np.random.randint(2))
			# This first copy step is necessary independent of conversion/recombination activity, 
			# to prevent changes to shared-referent chromosomes
			chrom = [[list(parChrom[0][0]),list(parChrom[0][1])],[list(parChrom[1][0]),list(parChrom[1][1])]]
			# Perform fly-specific non-recombination non-conversion check for males
			if self.willConvert and (self.sex == 'F' or (not self.maleAchiasmy)):
				chrom = self.__convertChrom(chrom)
			if (not self.willRecombine) or numRecomb == 0 or (self.maleAchiasmy and self.sex == 'M'):
				if currHom1:
					gamGenome += [chrom[0]]
				else:
					gamGenome += [chrom[1]]
			else:
				# NEED to allow for the span of positions in the chromosome if given positions outside [0,1)
				recombPositions = np.random.ranf(numRecomb)*self.lenChrom
				# recombPositions = np.random.ranf(numRecomb)
				recombPositions = np.sort(recombPositions)
				# print('Positions: '+str(recombPositions))
				
				# Construct a recombinant gamete from the recombination positions
				gamMutations = []
				gamInversions = []
				
				recIndex = 0

				invInd1 = 0
				invInd2 = 0
				invHom1 = chrom[0][1]
				invHom2 = chrom[1][1]
				mutInd1 = 0
				mutInd2 = 0

				prevInvInd1 = 0
				prevInvInd2 = 0
				prevMutInd1 = 0
				prevMutInd2 = 0

				while recIndex < numRecomb:
					recPos = recombPositions[recIndex]
					while (mutInd1 < len(chrom[0][0])) and (chrom[0][0][mutInd1][0] <= recPos):
						mutInd1 += 1
					while (mutInd2 < len(chrom[1][0])) and (chrom[1][0][mutInd2][0] <= recPos):
						mutInd2 += 1
					while (invInd1 < len(chrom[0][1])) and (chrom[0][1][invInd1][0] <= recPos):
						invInd1 += 1
					while (invInd2 < len(chrom[1][1])) and (chrom[1][1][invInd2][0] <= recPos):
						invInd2 += 1
					# Add the mutations and inversions to the gamete chromosome
					if currHom1:
						gamMutations += chrom[0][0][prevMutInd1:mutInd1]
						gamInversions += chrom[0][1][prevInvInd1:invInd1]
					else:
						gamMutations += chrom[1][0][prevMutInd2:mutInd2]
						gamInversions += chrom[1][1][prevInvInd2:invInd2]
					# Update the mutation/inversion indexes of the previous crossover
					prevInvInd1 = invInd1
					prevInvInd2 = invInd2
					prevMutInd1 = mutInd1
					prevMutInd2 = mutInd2
					# Deal with crossovers in regions potentially generating aneuploidy
					if self.oneInvPerChrom:
						aneupRegEnd = self.__getAneupRegionOneInvPerChrom(invHom1,invHom2,invInd1,invInd2,recPos)
					else:
						aneupRegEnd = self.__getAneupRegion(invHom1,invHom2,invInd1,invInd2,recPos)
					if aneupRegEnd == -1:
						currHom1 = not currHom1
						recIndex += 1
					else:
						otherRecInReg = []
						recIndex += 1
						while (recIndex < numRecomb) and (recombPositions[recIndex] < aneupRegEnd):
							otherRecInReg += [recombPositions[recIndex]]
							recIndex += 1
						# If there's an even number of recombinations in the region in total:
						if len(otherRecInReg)%2:
							currHom1 = not currHom1
							for recPos in otherRecInReg:
								while (mutInd1 < len(chrom[0][0])) and (chrom[0][0][mutInd1][0] <= recPos):
									mutInd1 += 1
								while (mutInd2 < len(chrom[1][0])) and (chrom[1][0][mutInd2][0] <= recPos):
									mutInd2 += 1
								# Add the mutations to the gamete chromosome (no inversion change in region)
								if currHom1:
									gamMutations += chrom[0][0][prevMutInd1:mutInd1]
								else:
									gamMutations += chrom[1][0][prevMutInd2:mutInd2]
								# Update the mutation indexes of the previous breakpoint
								prevMutInd1 = mutInd1
								prevMutInd2 = mutInd2
								# Switch the current chromosome from the crossover event
								currHom1 = not currHom1
				# Add the last mutations and inversions to the gamete chromosome
				if currHom1:
					gamMutations += chrom[0][0][prevMutInd1:len(chrom[0][0])]
					gamInversions += chrom[0][1][prevInvInd1:len(chrom[0][1])]
				else:
					gamMutations += chrom[1][0][prevMutInd2:len(chrom[1][0])]
					gamInversions += chrom[1][1][prevInvInd2:len(chrom[1][1])]
				gamGenome += [[gamMutations,gamInversions]]
		return gamGenome

	# For deciding which recombination breakpoint sampling to use
	def genGamete(self):
		resampleUnbalanced = False
		if resampleUnbalanced:
			return self.genGameteResampleCrossovers()
		else:
			return self.genGameteDropUnbCrossovers()

	# For printing a (potentially) command-line friendly version of the genome
	def printGenome(self, i, printMut = True, printInv = True):
		indStr = 'Ind {:6d} '.format(i)
		for c in range(len(self.genome)):
			for h in range(2):
				homStr = indStr + 'Chrom '+str(c)+' Hom '+str(h)
				if printMut:
					mutStr = homStr + ' Mut: ' + str(self.genome[c][h][0])
					print(mutStr)
				if printInv:
					invStr = homStr + ' Inv: ' + str(self.genome[c][h][1])
					print(invStr)
		return



# Module methods for preparing populations (genomes/sexes) from mutation/haplotype data
# Could be defined as Static methods inside the population class

# For testing the inputs for correct size parameters
# MAYY WANT TO TEST lenChrom
def __inputSizeCheck(size, numChrom, genomes, sexes):
	numGenomes = len(genomes)
	if numGenomes > 0:
		# Ignore numChrom if genomes pre-specified, peg it to the number of chormosomes in genomes[0]
		numChrom = len(genomes[0])
		if numGenomes > size:
			raise SAIpop.InputError(genomes,'More genomes provided than population size')
	# Either require equal length sexes and genomes or just that sexes is smaller than population size
	#   and allows for one of each sex?
	if (len(sexes) != 0) and (len(sexes) != numGenomes):
		raise SAIpop.InputError(sexes,'Sexes must be of equal length to genomes')
	return (numGenomes,numChrom)

# For generating the sexes list from a current sexes list 
def __genSexes(size, sexes, randomSex, minNumPerSex):
	# Account the provided sexes
	numFemales = 0
	numMales = 0
	for sex in sexes:
		if sex == 'M':
			numMales += 1
		elif sex == 'F':
			numFemales += 1
		else:
			raise SAIpop.InputError(sexes,'Sexes must only contain \'F\' and \'M\' elements')
	newSexes = []
	# Needs at least one of each sex to be a viable population
	while numMales < minNumPerSex:
		newSexes += ['M']
		numMales += 1
	while numFemales < minNumPerSex:
		newSexes += ['F']
		numFemales += 1
	numRemaining = size - numFemales - numMales
	if numRemaining < 0:
		raise SAIpop.InputError(sexes,'Sexes must allow at least one \'F\' and \'M\' in the population')
	# Reuse numMales for the number of males to be added
	if randomSex:
		numMales = np.random.binomial(numRemaining,0.5)
	else:
		numMales = int(numRemaining/2)
	newSexes += ['M' for i in range(numMales)]
	newSexes += ['F' for i in range(numRemaining-numMales)]
	np.random.shuffle(newSexes)
	popSexes = sexes + newSexes
	# Mutates the sexes list but will be returned anyway
	return popSexes

# For generating a set of genomes, sexes, and a record from specified mutations and inversions, 
# along with partial genomes and sexes lists to specify specific genomes,
# input inversion and mutation lists formatted as described in self.record but with count instead of initial generation:
# mutList: a list of all mutations, indexed by ID, with each entry as 
#      [position,survival effect,reproductive effect,chromosome,initial count]
# invList: a list of all inversions, indexed by ID, with each entry as 
#      [start position,end position,chromosome,initial count]
# populates the remaining genome slots with genomes to which mutations are randomly assigned from the count pools
# Currently doesn't allow description of prior record
#   How to model starting with mutations specific to sexes?
#   Allow limiting the number of further genomes to generate? How will this interact with having at least 1 'M'/'F'?
# @staticmethod
def genGenomesSexes(size, mutList, invList, randomSex = True, numChrom = 1, 
		lenChrom = 1.0, genomes = [], sexes = [], minNumPerSex = 1, verbose = False):
	# Handling input size checking
	if verbose: print('Num Chromosomes provided: '+str(numChrom))
	# print(genomes)
	# print(sexes)
	(numGenomes,numChrom) = __inputSizeCheck(size, numChrom, genomes, sexes) # UPDATE FOR NUMCHROM FROM MUT/INV LIST
	if verbose: print('Num Genomes already provided: '+str(numGenomes))
	if verbose: print('Num Chromosomes after checking if present in genomes provided: '+str(numChrom))
	# Generate the remaining sexes list
	popSexes = __genSexes(size, sexes, randomSex, minNumPerSex)

	# Generate the remaining genomes
	if verbose: print('Pop Size to be filled: '+str(size))
	numNewGenomes = size-numGenomes
	# print('Num Genomes not already provided: '+str(numNewGenomes))
	newGenomes = [[[[[],[]],[[],[]]] for i in range(numChrom)] for j in range(numNewGenomes)]
	# Will need to have separate lists per chromosome within these
	posOrderedMut = [[] for i in range(numChrom)]
	posOrderedInv = [[] for i in range(numChrom)]
	# Record should be returned with only [1,3] actually populated, [2,4...] with empty lists 
	#   for 0th generation update
	# record = [[],[],[],[],[],[],[],[],[],[],[]]
	# record = [[] for x in range(10)]
	record = [[],[],{},[],{},[],[],[],[],[],[],[],{},[],[],[],[]]

	numMut = len(mutList)
	for m in range(numMut):
		if verbose: print('Accounting mutation: '+str(mutList[m]))
		count = mutList[m][4]
		if count > 2*numNewGenomes:
			raise SAIpop.InputError(mutList,'Mutation counts must be less than the remaining population size')
		mutation = mutList[m][0:3]+[m]
		chrom = mutList[m][3]
		record[1] += [mutList[m][0:4]+[0,-1]]
		# record[2] += [[]]
		# record[2][m] = []
		# Insert the mutation and count to the ordered list for addition to the genomes
		i = 0
		while (i < len(posOrderedMut[chrom])) and (posOrderedMut[chrom][i][0][0] < mutation[0]):
			i += 1
		posOrderedMut[chrom][i:i] = [(mutation,count)]
	if verbose: print('Position-ordered mutations: '+str(posOrderedMut))
	numInv = len(invList)
	for i in range(numInv):
		if verbose: print('Accounting inversion: '+str(invList[i]))
		count = invList[i][3]
		if count > 2*numNewGenomes:
			raise SAIpop.InputError(invList,'Inversion counts must be less than the remaining population size')
		inversion = invList[i][0:2]+[i]
		chrom = invList[i][2]
		record[3] += [invList[i][0:3]+[0,-1]]
		# record[4] += [[]]
		# record[5] += [[]]
		# record[6] += [[]]
		# record[7] += [[]]
		# record[8] += [[]]
		# record[9] += [[]]
		# record[10] += [[]]
		# record[4][i] = [[],[],[],[],[],[],[],[],[],[]]
		# Insert the inversion and count to the ordered list for addition to the genomes
		i = 0
		while (i < len(posOrderedInv[chrom])) and (posOrderedInv[chrom][i][0][1] <= inversion[0]):
			i += 1
		if i < len(posOrderedInv[chrom]) and posOrderedInv[chrom][i+1][0][0] <= inversion[1]:
			raise SAIpop.InputError(invList,'Inversions cannot overlap')
		posOrderedInv[chrom][i:i] = [(inversion,count)]
	if verbose: print('Position-ordered inversions: '+str(posOrderedInv))

	# Preserve order by adding the mutations/inversions from position-ordered lists
	# Need to ensure no double-placement at any position
	possibleSpots = []
	newGenomes = []
	for i in range(numNewGenomes):
		for h in [0,1]:
			possibleSpots += [(i,h)]
		for chrom in range(numChrom):
			newGenomes += [[[[[],[]],[[],[]]]]]
	for chrom in range(numChrom):
		for (mut,count) in posOrderedMut[chrom]:
			np.random.shuffle(possibleSpots)
			for (i,hom) in possibleSpots[0:count]:
				newGenomes[i][chrom][hom][0] += [mut]
		for (inv,count) in posOrderedInv[chrom]:
			np.random.shuffle(possibleSpots)
			for (i,hom) in possibleSpots[0:count]:
				newGenomes[i][chrom][hom][1] += [inv]
	popGenomes = genomes + newGenomes
	# print(popGenomes)

	return (popGenomes,popSexes,record)

# For generating a set of genomes, sexes, and a record from specified mutations, inversions,
# and haplotypes with frequency data
# along with partial genomes and sexes lists to specify specific genomes,
# input lists formatted as described in self.record but with no initial generation:
# mutList: a list of all mutations, indexed by ID, with each entry as 
#      [position,survival effect,reproductive effect,chromosome]
# invList: a list of all inversions, indexed by ID, with each entry as 
#      [start position,end position,chromosome]
# hapList: a list of haplotype sets for the whole genome and counts, with each entry as 
#      [[hap chrom 0, hap chrom 1, ...], count]
# populates the remaining genome slots with genomes to which chromosomes are assigned from haplotypes randomly
# Currently doesn't allow description of prior record
#   How to model starting with mutations specific to sexes?
# IMPORTANT - FOR CORRECT RECORD KEEPING, ALL MUT OF THE SAME ID IN HAPS MUST BE SAME OBJECT/REFERENCE
# AND THE HAPLOTYPES MUST HAVE MUTATIONS AND INVERSIONS IN INCREASING POSITION (write a check function?)
# @staticmethod
def genGenoSexFromWholeGenHap(size, mutList, invList, hapList, randomSex = True, numChrom = 1, 
		lenChrom = 1.0, genomes = [], sexes = [], minNumPerSex = 1):
	# Handling input size checking
	(numGenomes,numChrom) = __inputSizeCheck(size, numChrom, genomes, sexes)
	# Check that the haplotype chromosome size matches the pre-specified genomes
	if (len(genomes) > 0) and (len(genomes[0]) != len(hapList[0][0])):
		# Ignore numChrom if genomes pre-specified, peg it to the number of chormosomes in genomes[0]
		raise SAIpop.InputError(hapList,'Haplotype lists must have chromosome number matching prespecified genomes')

	# Generate the remaining sexes list
	sexes = __genSexes(size, sexes, randomSex, minNumPerSex)

	# Generate the remaining genomes
	numNewGenomes = size-numGenomes

	#Build a haplotype pool to sample from
	totHapCount = 0
	hapPool = []
	for (haplotype,count) in hapList:
		totHapCount += count
		for i in range(count):
			hapCopy = []
			for chromHap in haplotype:
				hapCopy += [[list(chromHap[0]),list(chromHap[1])]]
			hapPool += [hapCopy]
	totHapSlots = 2*numNewGenomes
	# Check to ensure reasonable number of haplotypes provided
	if totHapCount > totHapSlots:
		raise SAIpop.InputError(hapList,'Haplotype counts larger than available population')
	numHapRemaining = totHapSlots - totHapCount
	# Fill the remaining slots with blank haplotypes
	for h in range(numHapRemaining):
		blankHap = []
		for c in range(numChrom):
			blankHap += [[[],[]]]
		hapPool += [blankHap]

	np.random.shuffle(hapPool)

	# Generate new genomes
	newGenomes = []
	for h in range(0,totHapSlots,2):
		genome = []
		for chrom in range(numChrom):
			genome += [[hapPool[h][chrom],hapPool[h+1][chrom]]]
		newGenomes += [genome]
	genomes += newGenomes

	# Record should be returned with only [1,3] actually populated, [2,4..7] with empty lists 
	#   for 0th generation update
	# record = [[],[],[],[],[],[],[],[],[],[],[]]
	# record = [[] for x in range(10)]
	record = [[],[],{},[],{},[],[],[],[],[],[],[],{},[],[],[],[]]
	for m in range(len(mutList)):
		record[1] += [mutList[m][0:4]+[0,-1]]
		# record[2] += [[]]
		# record[2][m] = []
	for i in range(len(invList)):
		record[3] += [invList[i][0:3]+[0,-1]]
		# record[4] += [[]]
		# record[5] += [[]]
		# record[6] += [[]]
		# record[7] += [[]]
		# record[8] += [[]]
		# record[9] += [[]]
		# record[10] += [[]]
		# record[4][i] = [[],[],[],[],[],[],[],[],[],[]]

	return (genomes,sexes,record)


# For loading pickled population simulations
def readPopPickle(filename):
	import pickle as p
	with open(filename,'rb') as infile:
		pop = p.load(infile)
	return pop


# The object representing the population under simulation, contains methods for stepping through generations
# size is the population size
# mutRate is the expected number of new mutations per chromosome per reproductive event
# mutRateInv is the expected number of new inversions per chromosome per reproductive event
# mutEffectDiffSD is the SD of normally distributed offset of survival/rep mutation effects
#    offsetting from [x,1-x], with x = uniform[0,1)
# minInvLen is the minimum length of a new inversion
# conversionRate is the probability heterozygous locus that a conversion occurs, direction is random
# recombRate is the expected number of recombination events per chromosome per reproductive event
# encounterNum is the number of males a female chooses from in a reproductive event
# choiceNoiseSD is the SD for the normally distributed additive noise to male quality in the female choice
# record is a general variable for storing generation statistics about the population, see __updateRecord
# invRecBuffer is the distance outside of inversion boundaries in which to record mutations for an inversion
class SAIpop(object):
	"""simSAIpopulation represents populations for sexually antagonistic inversion simulation"""
	def __init__(self, size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, conversionRate,
			recombRate, encounterNum, choiceNoiseSD, invRecBuffer,
			maleAchiasmy = True, randomSex = True, numChrom = 1, lenChrom = 1.0, willMutate = True,
			willMutInv = True, willConvert = True, willRecombine = True, 
			oneInvPerChrom = True, resampleInvMut = False, 
			removeFixed = False, recordReproductiveAdults = True, recordHaplotypes = None,
			minNumGensKeepRecord = 200, minFreqKeepMutRecord = 0.05, minFreqKeepInvRecord = 0.05, 
			noMaleCost = False, noFemaleCost = False, 
			sampleSurvivorsPerOffspring = True, sampleEncountersWithReplacement = True,
			verbose = False, testing = False, storePickledPop = False,
			genomes = [], sexes = [], record = [[],[],{},[],{},[],[],[],[],[],[],[],{},[],[],[],[]]):
		if verbose: print(record)
		# super(simSAIpopulation, self).__init__()
		self.size = size
		self.mutRate = mutRate
		self.mutRateInv = mutRateInv
		self.expectedNumMut = mutRate * size
		self.expectedNumInvMut = mutRateInv * size
		# Keep track of next mutation (or inversion) ID
		# self.__mutIDcount = 0
		# self.__invIDcount = 0
		self.__mutIDcount = len(record[1])
		self.__invIDcount = len(record[3])
		self.mutEffectDiffSD = mutEffectDiffSD
		self.encounterNum = encounterNum
		self.choiceNoiseSD = choiceNoiseSD
		self.minInvLen = minInvLen
		self.oneInvPerChrom = oneInvPerChrom
		self.resampleInvMut = resampleInvMut
		self.conversionRate = conversionRate
		self.recombRate = recombRate
		self.males = []
		self.females = []
		# For modeling D.mel specific recombination/other biology
		self.maleAchiasmy = maleAchiasmy
		self.randomSex = randomSex
		self.numChrom = numChrom
		self.lenChrom = lenChrom #FINISH IMPLEMENTATION
		# For modeling equlibrium dynamics without further mutation or other mechanics
		self.willMutate = willMutate
		self.willMutInv = willMutInv
		self.willConvert = willConvert
		self.willRecombine = willRecombine
		# For modeling survival and reproductive effects differently
		self.noMaleCost = noMaleCost
		self.noFemaleCost = noFemaleCost
		# To determine if mothers and fathers are sampled from the survival probabilities
		#   for each offspring, or if the survivors are sampled once and then parents
		#   are chosen uniformly from survivors
		self.sampleSurvivorsPerOffspring = sampleSurvivorsPerOffspring
		self.sampleEncountersWithReplacement = sampleEncountersWithReplacement
		# Record keeping variables, record data structure defined in __updateRecord
		self.recordReproductiveAdults = recordReproductiveAdults
		self.recordHaps = recordHaplotypes is not None
		self.removeFixed = removeFixed
		self.minNumGensKeepRecord = minNumGensKeepRecord
		self.minFreqKeepMutRecord = minFreqKeepMutRecord
		self.minFreqKeepInvRecord = minFreqKeepInvRecord
		self.invRecBuffer = invRecBuffer
		# State printing and testing variables
		self.verbose = verbose
		self.testing = testing
		self.storePickledPop = storePickledPop
		# if self.storePickledPop:
		# 	import pickle as p
		# Internal accounting variables
		self.__invFixed = [False for i in range(self.__invIDcount)]
		self.__locFixed = [False for m in range(self.__mutIDcount)]
		self.__invLost = [False for i in range(self.__invIDcount)]
		self.__mutLost = [False for m in range(self.__mutIDcount)]
		self.record = record
		self.arrangementRecord = []
		if self.recordHaps:
			for hapDef in recordHaplotypes:
				self.record[11] += [hapDef+[0,0]]
		self.age = 0
		# Handling input errors
		if self.encounterNum > int(self.size/2):
			raise SAIpop.InputError((self.size,self.encounterNum),
				'The number of males encountered cannot be more than half the pop size.'+\
				'\nPop size\t'+str(self.size)+'\nEncounter #\t'+str(self.encounterNum)+\
				'\nEnc/Pop Ratio\t\t'+str(self.encounterNum/self.size))
		if self.encounterNum > int(self.size/10):
			print('WARNING - The number of males encountered should be much smaller '+\
				'than half of the pop size, as otherwise the simulation may force males '+\
				'to survive in order to fulfill the needed encounter number.\n'+\
				'Pop size\t'+str(self.size)+'\nEncounter #\t'+str(self.encounterNum)+\
				'\nEnc/Pop Ratio\t\t'+str(self.encounterNum/self.size)+\
				'\nSample Survivors per Offspring:\t\t'+str(self.sampleSurvivorsPerOffspring)+\
				'\nSample Encounters with Replacement:\t'+str(self.sampleEncountersWithReplacement))
		# Handling default populations and related input errors
		numGenomes = len(genomes)
		if numGenomes > size:
			raise SAIpop.InputError(genomes,'More genomes provided than population size')
		elif len(sexes) != 0 and len(sexes) != numGenomes:
			raise SAIpop.InputError(sexes,'Sexes must be of equal length to genomes')
		else:
			if numGenomes > 0:
				# Ignore numChrom if genomes pre-specified, peg it to the number of chormosomes in genomes[0]
				self.numChrom = len(genomes[0])
				for i in range(numGenomes):
					if sexes[i] == 'F':
						self.females += [individual('F',mutEffectDiffSD,recombRate,conversionRate,\
							minInvLen,lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,genomes[i])]
					else:
						self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,\
							minInvLen,lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,genomes[i])]
			numLeft = size-numGenomes
			if len(self.females) < 1:
				numLeft -= 1
				self.females += [individual('F',mutEffectDiffSD,recombRate,conversionRate,\
					minInvLen,lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,\
					[[[[],[]],[[],[]]] for i in range(self.numChrom)])]
			if len(self.males) < 1:
				numLeft -= 1
				self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,\
					minInvLen,lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,\
					[[[[],[]],[[],[]]] for i in range(self.numChrom)])]
			# Generate the remaining set of 'size' individuals
			#   with random or equal probability sex and no mutations or inversions
			if randomSex:
				numMales = np.random.binomial(numLeft,0.5)
			else:
				numMales = int(numLeft/2)
			self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
				lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,\
				[[[[],[]],[[],[]]] for i in range(self.numChrom)]) for i in range(numMales)]
			self.females += [individual('F',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
				lenChrom,maleAchiasmy,willConvert,willRecombine,oneInvPerChrom,\
				[[[[],[]],[[],[]]] for i in range(self.numChrom)]) for j in range(numLeft-numMales)]
			# Generate the generation 0 record
			# print("Gen 0 Record: "+str(record))
			self.__updateRecord()
			# print("Update: "+str(self.record))


	class InputError(Exception):
		"""Exception raised for errors in the input.

		Attributes:
			expression -- input expression in which the error occurred
			message -- explanation of the error
		"""

		def __init__(self, expression, message):
			# super(self).__init__(message)
			Exception.__init__(self,message)
			self.expression = expression
			self.message = message

	class SimulationStateError(Exception):
	    """Exception raised for simulation states that should be impossible.
	    Only some such states are checked.

	    Attributes:
	        expression -- input expression in which the error occurred
	        message -- explanation of the error
	    """

	    def __init__(self, expression, message):
	        # super(self).__init__(message)
	        Exception.__init__(self,message)
	        self.expression = expression
	        self.message = message

	# For changing the population size during simulation
	def setSize(self,newSize):
		self.size = newSize
		return

	# For extracting the population sizes during recorded generations
	#  that overlap with the period specified
	def __getPopSizesDuringPeriod(self,initGen,finalGen):
		numRecords = len(self.record[0])
		i = 0
		while i < numRecords and self.record[0][i] < initGen:
			i += 1
		j = i
		while j < numRecords and self.record[0][j] < finalGen:
			j += 1
		popSizes = self.record[5][i:j]
		return popSizes


	# For extracting the generations at which recording occurred
	#  that overlap with the period specified
	def __getRecordGensDuringPeriod(self,initGen,finalGen):
		numRecords = len(self.record[0])
		i = 0
		while i < numRecords and self.record[0][i] < initGen:
			i += 1
		j = i
		while j < numRecords and self.record[0][j] < finalGen:
			j += 1
		recGens = self.record[0][i:j]
		return recGens

	# For extracting the life history stages at which recording occurred
	#  that overlap with the period specified
	def __getRecordLHSDuringPeriod(self,initGen,finalGen):
		numRecords = len(self.record[0])
		i = 0
		while i < numRecords and self.record[0][i] < initGen:
			i += 1
		j = i
		while j < numRecords and self.record[0][j] < finalGen:
			j += 1
		recLifeHistStages = self.record[10][i:j]
		return recLifeHistStages

	# For extracting the population sizes during recorded generations
	#  that overlap with the period specified
	def __getPopSizesAndCheckLHS(self,initGen,finalGen,counts,ID,isMut=True):
		numRecords = len(self.record[0])
		i = 0
		while i < numRecords and self.record[0][i] < initGen:
			i += 1
		j = i
		numPopSizes = 0
		while j < numRecords and self.record[0][j] < finalGen:
			j += 1
			numPopSizes += 1
		popSizes = self.record[5][i:j]
		numCounts = len(counts)

		if numCounts == numPopSizes:
			return popSizes
		# Accounting for when the mutation is lost in the previous generation,
		#    after zygote fixation check but before adult/later LHS recording
		# ONLY 1 LHS CURRENTLY, will need to be changed if more added
		elif numCounts == numPopSizes - 1 and self.record[10][j-1] != 'Zygote':
			if self.verbose:
				if isMut:
					print("WARNING - mutation appears to have been lost after fixation check but before adult recording")
					print("Removing Mut "+str(ID))
					print(self.record[1][ID])
					print("During generation "+str(self.age))
					print("Mut Count Record:")
					print(counts)
					print("Recovered pop sizes:")
					print(popSizes)
					print("Recovered life history stages:")
					print(self.__getRecordLHSDuringPeriod(initGen,finalGen))
				else:
					print("WARNING - inversion appears to have been lost after fixation check but before adult recording")
					print("Removing Inv "+str(ID))
					print(self.record[3][ID])
					print("During generation "+str(self.age))
					print("Inv Count Record:")
					print(counts)
					print("Recovered pop sizes:")
					print(popSizes)
					print("Recovered life history stages:")
					print(self.__getRecordLHSDuringPeriod(initGen,finalGen))
			return popSizes[:-1]
		else:
			if self.verbose:
				if isMut:
					print("Removing Mut "+str(ID))
					print(self.record[1][ID])
					print("During generation "+str(self.age))
					print("Mut Count Record:")
					print(counts)
					print("Recovered pop sizes:")
					print(popSizes)
					print("Recovered life history stages:")
					print(self.__getRecordLHSDuringPeriod(initGen,finalGen))
				else:
					print("Removing Inv "+str(ID))
					print(self.record[3][ID])
					print("During generation "+str(self.age))
					print("Inv Count Record:")
					print(counts)
					print("Recovered pop sizes:")
					print(popSizes)
					print("Recovered life history stages:")
					print(self.__getRecordLHSDuringPeriod(initGen,finalGen))
			return popSizes


	# Checks for mutation age and frequency thresholds, 
	#  removing the mutation records from the record dict if these aren't met
	def __mutRecordDictRemoval(self,ID):
		if self.testing:
			popSizes = self.__getPopSizesDuringPeriod(self.record[1][ID][4],self.record[1][ID][5])
			recGens = self.__getRecordGensDuringPeriod(self.record[1][ID][4],self.record[1][ID][5])
			if ID not in self.record[2]:
				assert len(popSizes) == len(recGens), "Number of recorded generations and population sizes "+\
					"at those generations do not have the same number of entries"
				if len(popSizes) > 0:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes), "Error: attempting to remove a mutation record "+\
						"when no record is present for mutation "+str(ID)+", despite recording population sizes "+\
						repr(popSizes)+" at generations "+repr(recGens))
				else:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes), "Error: attempting to remove a mutation record "+\
						"when no record is present from instatiation of the population or mutation, or the "+\
						"record was already removed."+\
						"\nPop sizes: "+str(popSizes))
			else:
				if len(popSizes) > len(self.record[2][ID][0]) + 1:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes), "Error: attempting to remove the record " +\
						"of mutation "+str(ID)+"; "+str(self.record[1][ID])+" and encountered a significant length difference in "+\
						"the number of expected records and the number of counts recorded."+\
						"\nPop sizes: "+str(popSizes)+"\nCounts: "+str(self.record[2][ID][0]))
		if ID in self.record[2]:
			initGen = self.record[1][ID][4]
			finalGen = self.record[1][ID][5]
			numGensPolymorphic = finalGen - initGen
			mutCounts = self.record[2][ID][0]
			maxFreq = 0
			if len(mutCounts) > 0:
				popSizes = self.__getPopSizesAndCheckLHS(initGen,finalGen,mutCounts,ID)
				maxFreq = np.max(np.divide(mutCounts,popSizes))
			if numGensPolymorphic < self.minNumGensKeepRecord and maxFreq < self.minFreqKeepMutRecord:
				del self.record[2][ID]
		# else:
		# 	raise self.SimulationStateError((ID), "Error: attempting to remove a mutation record "+\
		# 		"when no record is present for mutation "+str(ID))
		return

	# Checks for inversion age and frequency thresholds, 
	#  removing the inversion records from the record dict if these aren't met
	def __invRecordDictRemoval(self,ID):
		if self.testing:
			popSizes = self.__getPopSizesDuringPeriod(self.record[3][ID][3],self.record[3][ID][4])
			recGens = self.__getRecordGensDuringPeriod(self.record[3][ID][3],self.record[3][ID][4])
			if ID not in self.record[4]:
				assert len(popSizes) == len(recGens), "Number of recorded generations and population sizes "+\
					"at those generations do not have the same number of entries"
				if len(popSizes) > 0:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes), "Error: attempting to remove a inversion record "+\
						"when no record is present for inversion "+str(ID)+", despite recording population sizes "+\
						repr(popSizes)+" at generations "+repr(recGens))
				else:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes), "Error: attempting to remove a inversion record "+\
						"when no record is present from instantiation of the population or inversion, or the "+\
						"record was already removed."+\
						"\nPop sizes: "+str(popSizes))
			else:
				if len(popSizes) > len(self.record[4][ID][0]) +1:
					if self.storePickledPop:
						self.__storePickleBytes(p.dumps(self),'errorPop.pickle')
					raise self.SimulationStateError((ID,popSizes,self.record[4][ID][0]), "Error: attempting to remove the record " +\
						"of inversion "+str(ID)+"; "+str(self.record[3][ID])+" and encountered a significant length difference in "+\
						"the number of expected records and the number of counts recorded."+\
						"\nPop sizes: "+str(popSizes)+"\nCounts: "+str(self.record[4][ID][0]))

		if ID in self.record[4]:
			initGen = self.record[3][ID][3]
			finalGen = self.record[3][ID][4]
			numGensPolymorphic = finalGen - initGen
			invCounts = self.record[4][ID][0]
			maxFreq = 0
			if len(invCounts) > 0:
				popSizes = self.__getPopSizesAndCheckLHS(initGen,finalGen,invCounts,ID,isMut=False)
				maxFreq = np.max(np.divide(invCounts,popSizes))
			if numGensPolymorphic < self.minNumGensKeepRecord and maxFreq < self.minFreqKeepInvRecord:
				del self.record[4][ID]
		return


	# # Removes inversions that are fixed in the population,
	# # changes the position of internal mutations to reflect true position (flips around the center of the inversion)
	# # the old position of the mutations is lost
	# # WORTH IMPROVING PERFORMANCE: merge with record keeping in record generations? also reduce loop number
	# # ^ Just add a record flag to the step, you *could* mask this with step{step(false)}, recordStep{step(true)}
	# def __removeFixedInv(self):
	# 	wholePop = self.males + self.females
	# 	invCounts = [0]*self.__invIDcount
	# 	invPos = [[]]*self.__invIDcount
	# 	# All mutations are passed as references, so don't want to change more than one
	# 	firstMutEncounter = [True]*self.__mutIDcount

	# 	inversionRemoved = False
	# 	# Generate count and position data
	# 	for i in range(len(wholePop)):
	# 		for c in range(self.numChrom):
	# 			for h in range(2):
	# 				chromHomInv = wholePop[i].genome[c][h][1]
	# 				for invIndex in range(len(chromHomInv)):
	# 					inversion = chromHomInv[invIndex]
	# 					# print(inversion)
	# 					invCounts[inversion[2]] += 1
	# 					invPos[inversion[2]] += [[i,c,h,invIndex]]
	# 	for ID in range(self.__invIDcount):
	# 		if invCounts[ID] == 2*self.size:
	# 			inversionRemoved = True
	# 			# print("Inversion Removed")
	# 			self.__invFixed[ID] = True
	# 			inversion = self.record[3][ID]
	# 			length = self.record[3][ID][1] - self.record[3][ID][0]
	# 			invCenter = self.record[3][ID][0] + length/2.0
	# 			for [i,c,h,invIndex] in invPos[ID]:
	# 				# Remove the inversion
	# 				wholePop[i].genome[c][h][1][invIndex:invIndex+1] = []
	# 				# Rearrange the mutation positions to reflect the loss of inversion data
	# 				mutHom = wholePop[i].genome[c][h][0]
	# 				index = 0
	# 				mutInside = []
	# 				while index < len(mutHom) and mutHom[index][0] < inversion[0]:
	# 					index += 1
	# 				startIndex = index
	# 				while index < len(mutHom) and mutHom[index][0] < inversion[1]:
	# 					mutInside = [mutHom[index]] + mutInside
	# 					index += 1
	# 				# Replace the mutations inside with ones in inverted order
	# 				mutHom[startIndex:index] = mutInside
	# 				# Flip the mutation positions around the center of the removed inversion
	# 				for mut in mutInside:
	# 					if firstMutEncounter[mut[3]]:
	# 						firstMutEncounter[mut[3]] = False
	# 						# Flip the mutation across the center
	# 						mut[0] = 2*invCenter-mut[0]
	# 					# print(mut)
	# 		elif invCounts[ID] > 2*self.size:
	# 			print("Error: more inversions than alleles possible for ID " + str(i))
	# 	return inversionRemoved

	# Removes mutations and inversions that are fixed in the population,
	# changes the position of internal mutations to reflect true position (flips around the center of the inversion)
	# the old position of the mutations is lost
	# WORTH IMPROVING PERFORMANCE: merge with record keeping in record generations? and generating new individuals. also reduce loop number
	# ^ Just add a record flag to the step, you *could* mask this with step{step(false)}, recordStep{step(true)}
	def __removeFixed(self):  # LOVE <3
		if self.testing:
			if self.storePickledPop:
				import pickle as p
				initPopBytes = p.dumps(self)
			else:
				initPopHapStr = self.__getStrPopChromCounts()
		wholePop = self.males + self.females
		mutCounts = np.zeros(self.__mutIDcount)
		mutPos = [[] for i in range(self.__mutIDcount)]
		invCounts = np.zeros(self.__invIDcount)
		invPos = [[] for i in range(self.__invIDcount)]
		# Generate count and position data
		for i in range(len(wholePop)):
			for c in range(self.numChrom):
				for h in range(2):
					chromHomMut = wholePop[i].genome[c][h][0]
					for mutIndex in range(len(chromHomMut)):
						mutation = chromHomMut[mutIndex]
						mutCounts[mutation[3]] += 1
						mutPos[mutation[3]] += [[i,c,h,mutIndex]]
					chromHomInv = wholePop[i].genome[c][h][1]
					for invIndex in range(len(chromHomInv)):
						inversion = chromHomInv[invIndex]
						# print(inversion)
						invCounts[inversion[2]] += 1
						invPos[inversion[2]] += [[i,c,h,invIndex]]
		if self.verbose:
			if np.sum(mutCounts==(2*len(wholePop)))>1:
				print("Multiple alleles fixing in one generation:")
			print(mutCounts)
			for ID in np.arange(self.__mutIDcount):
				if mutCounts[ID] > 0:
					print("Mut\t"+repr(ID)+"\tCount\t"+repr(mutCounts[ID]))

		# Remove fixed mutations
		mutationRemoved = False
		mutsToBeRemoved = []
		for ID in range(self.__mutIDcount):
			count = mutCounts[ID]
			if self.__locFixed[ID]:
				if count != 0:
					raise self.SimulationStateError((ID,count,self.age), "Error: mutation ID "+str(ID)+\
						" fixation flag was already set, but it is still present with count "+\
						str(count)+" in generation "+str(self.age))
			else:
				# print("Mut "+str(ID)+" Count "+str(count))
				if count == 2*self.size:
					mutationRemoved = True
					# print("Mutation Removed, mutFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(mutCounts[ID]))
					# print(str(len(mutPos[ID]))+" mutPos entries")
					# Update the fix/loss record
					self.__locFixed[ID] = True
					if self.record[1][ID][5] != -1:
						# should mean that self.__mutLost[ID] = True
						raise self.SimulationStateError((ID,count,self.age), "Error: mutation ID "+str(ID)+\
							" final generation was already set, but it is being set again with count "+\
							str(count)+" in generation "+str(self.age))
					self.record[1][ID][5] = self.age
					self.__mutRecordDictRemoval(ID)
					mutsToBeRemoved += [ID]
				elif count == 0:
					# print("Mutation Lost, mutFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(mutCounts[ID]))
					# print(str(len(mutPos[ID]))+" mutPos entries")
					# Update the fix/loss record
					self.__locFixed[ID] = True
					self.__mutLost[ID] = True
					if self.record[1][ID][5] != -1:
						# should mean that self.__mutLost[ID] = True
						raise self.SimulationStateError((ID,count,self.age), "Error: mutation ID "+str(ID)+\
							" final generation was already set, but it is being set again with count "+\
							str(count)+" in generation "+str(self.age))
					self.record[1][ID][5] = self.age
					self.__mutRecordDictRemoval(ID)
				elif count > 2*self.size:
					raise self.SimulationStateError((ID,count,2*self.size), "Error: more mutations than alleles possible for ID "+\
						str(ID)+"in population size "+str(self.size)+": Count "+str(count))
					# self.printGenomes()
				elif count < 0:
					raise self.SimulationStateError((ID,count), "Error: negative count of mutations ID "+\
						str(ID)+": Count "+str(count))
					# self.printGenomes()
		# Sort the mutations that need removal by position in the chromosome,
		#   so that removing one mutation will not change the position of subsequent ones to be removed
		sortedMuts = sorted(mutsToBeRemoved,key=lambda ID:-self.record[1][ID][0])
		if self.verbose:
			print("Mutations flagged for removal:\n"+repr(mutsToBeRemoved))
			print("Sorted in removal order:\n"+repr(sortedMuts))
			print("Chromosomal position of mutations:\n"+repr([self.record[1][ID][0] for ID in mutsToBeRemoved]))

		# print(sortedMuts)
		for ID in sortedMuts:
			# Remove the mutation from all population chromosomes
			if self.verbose:
				print("Removing mutation: "+str(self.record[1][ID])+" ID: "+str(ID))
			for [i,c,h,mutIndex] in mutPos[ID]:
				# print('Attempting removal at '+str([i,c,h,mutIndex]))
				# print(wholePop[i].genome[c][h][0])
				wholePop[i].genome[c][h][0][mutIndex:mutIndex+1] = []
				# print(wholePop[i].genome[c][h][0])

		if self.testing:
			# Testing to ensure removal of mutations
			testMutCounts = np.zeros(self.__mutIDcount)
			for i in np.arange(len(wholePop)):
				for c in np.arange(self.numChrom):
					for h in np.arange(2):
						# chromHomMut = wholePop[i].genome[c][h][0]
						for mutation in wholePop[i].genome[c][h][0]:
						# for mutIndex in range(len(chromHomMut)):
						# 	mutation = chromHomMut[mutIndex]
							testMutCounts[mutation[3]] += 1
							if mutation[0] != self.record[1][mutation[3]][0] or \
								mutation[1] != self.record[1][mutation[3]][1] or \
								mutation[2] != self.record[1][mutation[3]][2]:
								if self.storePickledPop:
									self.__storePickleBytes(initPopBytes,'errorPop.pickle')
								raise self.SimulationStateError((mutation,self.record[1][mutation[3]]), \
									"Mutation data found in chromosome does not match mutation record. ID "+\
									repr(mutation[3])+"; Chromosomal copy "+repr(mutation)+"; Record copy "+\
									repr(self.record[1][mutation[3]]))
			for j in np.arange(self.__mutIDcount):
				if self.__locFixed[j] == True and testMutCounts[j] != 0:
					errorStr = "Error: fixed mutation intended to be "+\
						"removed from the simulation has nonzero count:\nID "+\
						str(j)+": Count "+str(mutCounts[j])+": Recount "+str(testMutCounts[j])+\
						": Mutation "+str(self.record[1][j])+": Record "+str(self.record[2].get(j))+'\n'
					if self.storePickledPop:
						self.__storePickleBytes(initPopBytes,'errorPop.pickle')
						raise self.SimulationStateError((j,mutCounts[j],testMutCounts[j]), errorStr)
					else:
						raise self.SimulationStateError((j,mutCounts[j],testMutCounts[j]), errorStr +\
							"Population haplotype summary before removals:\n"+initPopHapStr+\
							'\nPopulation haplotype summary after mutation removal:\n'+self.__getStrPopChromCounts()+'\n')
			# print(str(testMutCounts[ID])+" remaining after removal")
		# if self.storePickledPop and self.age > 200:
		# 	self.__storePickleBytes(initPopBytes,'errorPop.pickle')
		# 	raise self.SimulationStateError(self.age, "simulated a stop to test pickling")

		if self.verbose:
			if np.sum(invCounts==(2*len(wholePop)))>1:
				print("Multiple inversions fixing in one generation:")
			print(invCounts)
			for ID in np.arange(self.__invIDcount):
				if invCounts[ID] > 0:
					print("Inv\t"+repr(ID)+"\tCount\t"+repr(invCounts[ID]))

		# Remove fixed inversions
		inversionRemoved = False
		# For use when reordering mutations on the chromosome when a fixed inversion is removed, because
		#   all mutations are passed as references, so don't want to change the reference more than once
		# Also, no mutations will be in multiple fixed inversions, so the scope is fine :P
		firstMutEncounter = [True]*self.__mutIDcount
		for ID in range(self.__invIDcount):
			count = invCounts[ID]
			if self.__invFixed[ID]:
				if count != 0:
					raise self.SimulationStateError((ID,count,self.age), "Error: inversion ID "+str(ID)+\
						" fixation flag was already set, but it is but it is still present with count "+\
						str(count)+" in generation "+str(self.age))
			else:
				if count == 2*self.size:
					# print("Inversion fixed: ID "+str(ID))
					inversionRemoved = True
					# print("Inversion Removed, invFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(invCounts[ID]))
					# Update the fix/loss record
					self.__invFixed[ID] = True
					if self.record[3][ID][4] != -1:
						# should mean that self.__mutLost[ID] = True
						raise self.SimulationStateError((ID,count,self.age), "Error: inversion ID "+str(ID)+\
							" final generation was already set, but it is being set again with count "+\
							str(count)+" in generation "+str(self.age))
					self.record[3][ID][4] = self.age
					self.__invRecordDictRemoval(ID)
					# Remove the inversion,
					# and reorder all variants within the inversion in order to account
					#    for the new default orientation when the inversion is removed
					inversion = self.record[3][ID]
					length = inversion[1] - inversion[0]
					invCenter = inversion[0] + length/2.0
					for [i,c,h,invIndex] in invPos[ID]:
						# print("Removal Attempt")
						# print(wholePop[i].genome[c][h][1])
						wholePop[i].genome[c][h][1][invIndex:invIndex+1] = []
						# print(wholePop[i].genome[c][h][1])
						# Rearrange the mutation positions to reflect the loss of inversion orientation data
						mutHom = wholePop[i].genome[c][h][0]
						index = 0
						mutInside = []
						while index < len(mutHom) and mutHom[index][0] < inversion[0]:
							index += 1
						startIndex = index
						while index < len(mutHom) and mutHom[index][0] < inversion[1]:
							mutInside = [mutHom[index]] + mutInside
							index += 1
						# Replace the mutations inside with ones in inverted order
						mutHom[startIndex:index] = mutInside
						# Flip the mutation positions around the center of the removed inversion
						for mut in mutInside:
							if firstMutEncounter[mut[3]]:
								firstMutEncounter[mut[3]] = False
								# Flip the mutation across the center: pos + 2*(center-pos)
								mut[0] = 2*invCenter-mut[0]
								# Update the mutation record to reflect the new position
								self.record[1][mut[3]][0] = mut[0]
							# print(mut)
				elif count == 0:
					# print("Inversion Lost, invFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(invCounts[ID]))
					# Update the fix/loss record
					self.__invFixed[ID] = True
					self.__invLost[ID] = True
					if self.record[3][ID][4] != -1:
						# should mean that self.__mutLost[ID] = True
						raise self.SimulationStateError((ID,count,self.age), "Error: inversion ID "+str(ID)+\
							" final generation was already set, but it is being set again with count "+\
							str(count)+" in generation "+str(self.age))
					self.record[3][ID][4] = self.age
					self.__invRecordDictRemoval(ID)
				elif count > 2*self.size:
					raise self.SimulationStateError((ID,count,2*self.size), "Error: more inversions than alleles possible for ID "+\
						str(ID)+"in population size "+str(self.size)+": Count "+str(count))
					# self.printGenomes()
				elif count < 0:
					raise self.SimulationStateError((ID,count), "Error: negative count of inversion ID "+\
						str(ID)+": Count "+str(count))
					# self.printGenomes()

		if self.testing:
			# Testing to ensure removal of inversions
			testInvCounts = np.zeros(self.__invIDcount)
			for i in np.arange(len(wholePop)):
				for c in np.arange(self.numChrom):
					for h in np.arange(2):
						# chromHomInv = wholePop[i].genome[c][h][1]
						# for invIndex in range(len(chromHomInv)):
						# 	inversion = chromHomInv[invIndex]
						for inversion in wholePop[i].genome[c][h][1]:
							testInvCounts[inversion[2]] += 1
							if inversion[0] != self.record[3][inversion[2]][0] or \
								inversion[1] != self.record[3][inversion[2]][1]:
								if self.storePickledPop:
									self.__storePickleBytes(initPopBytes,'errorPop.pickle')
								raise self.SimulationStateError((inversion,self.record[3][inversion[2]]), \
									"Inversion data found in chromosome does not match mutation record. ID "+\
									repr(inversion[2])+"; Chromosomal copy "+repr(inversion)+"; Record copy "+\
									repr(self.record[3][inversion[2]]))
			for j in np.arange(self.__invIDcount):
				if self.__invFixed[j] == True and testInvCounts[j] != 0:
					errorStr = "Error: fixed inversion intended to be "+\
						"removed from the simulation has nonzero count:\nID "+\
						str(j)+": Count "+str(invCounts[j])+": Recount "+str(testInvCounts[j])+\
						": Inversion "+str(self.record[3][j])+": Record "+str(self.record[4].get(j))+'\n'
					if self.storePickledPop:
						self.__storePickleBytes(initPopBytes,'errorPop.pickle')
						raise self.SimulationStateError((j,invCounts[j],testInvCounts[j]), errorStr)
					else:
						raise self.SimulationStateError((j,invCounts[j],testInvCounts[j]), errorStr +\
							"Population haplotype summary before removals:\n"+initPopHapStr+\
							'\nPopulation haplotype summary after mutation removal:\n'+self.__getStrPopChromCounts()+'\n')
		return (mutationRemoved, inversionRemoved)

	# Simulating a single generational step
	def step(self,recReprodParents=False,recParentalArr=False):
		# print ("Generation " + str(self.age))
		if self.verbose:
			print ("Generation " + str(self.age))
			print("Current Haplotypes:")
			self.__printPopChromCounts()
		# Pick mothers
		motherIndexes = []
		if self.noFemaleCost:
			motherIndexes = np.random.randint(low=0, high=len(self.females), size=self.size)
		else:
			numFemales = len(self.females)
			femaleSurvivalProbs = np.ones(numFemales)
			for f in np.arange(numFemales):
				# print(female.genome)
				femaleSurvivalProbs[f] = self.females[f].survival()
			if self.verbose:
				print("Mean female survival:\t\t\t"+str(np.mean(femaleSurvivalProbs)))
			# print("Female Survivals:" + str(femaleSurvivalProbs))
			# If mothers and fathers are sampled from the survival probabilities separately for each offspring
			if self.sampleSurvivorsPerOffspring:
				scaledFemaleSurvivalProbs = np.divide(femaleSurvivalProbs,np.sum(femaleSurvivalProbs))
				motherIndexes = np.random.choice(numFemales,self.size,p=scaledFemaleSurvivalProbs)
			# If the surviving parents are sampled once and then parents are chosen uniformly from survivors
			else:
				survivorFlags = femaleSurvivalProbs >= np.random.uniform(size=numFemales)
				numSurvivors = np.sum(survivorFlags)
				while numSurvivors < 1:
					print("WARNING - No surviving females in survivor sampling. Resampling and continuing.")
					survivorFlags = femaleSurvivalProbs >= np.random.uniform(size=numFemales)
					numSurvivors = np.sum(survivorFlags)
				survivorUnifProb = np.divide(survivorFlags,numSurvivors)
				motherIndexes = np.random.choice(numFemales,self.size,p=survivorUnifProb)


		# Calculate necessary male values - saves time to calculate male reproductive quality here
		numMales = len(self.males)
		maleSurvivalProbs = np.ones(numMales)
		maleRepQualities = np.zeros(numMales)
		if self.noMaleCost:
			if self.sampleSurvivorsPerOffspring:
				scaledMaleSurvivalProbs = np.repeat(1/numMales,numMales)
			else:
				survivingMaleUnif = np.repeat(1/numMales,numMales)
			for m in np.arange(numMales):
				maleRepQualities[m] = self.males[m].repQuality()
		else:
			for m in np.arange(numMales):
				maleSurvivalProbs[m] = self.males[m].survival()
			# maleGenQuality += [male.repQuality()]
			if self.verbose:
				print("Mean male survival:\t\t\t"+str(np.mean(maleSurvivalProbs)))
			if self.sampleSurvivorsPerOffspring:
				scaledMaleSurvivalProbs = np.divide(maleSurvivalProbs,np.sum(maleSurvivalProbs))
				for m in np.arange(numMales):
					maleRepQualities[m] = self.males[m].repQuality()
				if self.verbose:
					print("Mean male repr. quality:\t\t"+str(np.mean(maleRepQualities)))
			else:
				survivorFlags = maleSurvivalProbs >= np.random.uniform(size=numMales)
				numSurvivors = np.sum(survivorFlags)
				# CHECK - allow male encounters with replacement???
				if self.sampleEncountersWithReplacement:
					requiredNumMales = 1
					warnStr = " Resampling until at least one male survives, and continuing."
				else:
					requiredNumMales = self.encounterNum
					warnStr = " Resampling until the number of survivng males is >= the encounter size ("+\
						str(self.encounterNum)+"), and continuing."
				while numSurvivors < requiredNumMales:
					print("WARNING - "+str(numSurvivors)+" surviving males in survivor sampling."+warnStr)
					survivorFlags = maleSurvivalProbs >= np.random.uniform(size=numMales)
					numSurvivors = np.sum(survivorFlags)
				survivingMaleUnif = np.divide(survivorFlags,numSurvivors)
				for m in np.arange(numMales):
					# Save a little time by only calculating for the survivors that will be sampled
					if survivorFlags[m]:
						maleRepQualities[m] = self.males[m].repQuality()
				if self.verbose:
					print("Mean (surviving) male repr. quality:\t"+str(np.mean(maleRepQualities[survivorFlags])))
		# print("Male Survival Probabilities: " + str(maleSurvivalProbs))

		# Populate the next generation by choosing fathers per mother and generating child genomes
		newMales = []
		newFemales = []

		if recReprodParents or recParentalArr:
			# Get the set of mothers/fathers and their number of offspring, for recording the reproductive population statistics
			obsSurvFemInds = set(motherIndexes)
			obsSurvMalInds = set()
			motherCounter = Counter(motherIndexes)
			fatherCounter = Counter()

		# numExpectedMales used for pegged equal-sex populations,
		#   and to guarantee at least one of each sex per generation
		if self.randomSex:
			if self.sampleEncountersWithReplacement:
				numExpectedMales = 1
				numExpectedFemales = 1
			else:
				numExpectedMales = self.encounterNum
				numExpectedFemales = self.encounterNum
		else:
			numExpectedMales = int(self.size/2)
			numExpectedFemales = self.size - numExpectedMales
		# For each mother, select a father
		for motherIndex in motherIndexes:
			mother = self.females[motherIndex]
			if self.sampleSurvivorsPerOffspring:
				encounteredMaleIndexes = np.random.choice(numMales,self.encounterNum,
					replace=self.sampleEncountersWithReplacement,p=scaledMaleSurvivalProbs)
				if recReprodParents or recParentalArr:
					obsSurvMalInds.update(set(encounteredMaleIndexes))
			else:
				encounteredMaleIndexes = np.random.choice(numMales,self.encounterNum,
					replace=self.sampleEncountersWithReplacement,p=survivingMaleUnif)
				if recReprodParents or recParentalArr:
					obsSurvMalInds.update(set(encounteredMaleIndexes))
			# fatherIndex = encounteredMaleIndexes[0]
			# # maxScore = self.males[fatherIndex].repQuality() + np.random.normal(scale=self.choiceNoiseSD)
			# maxScore = 0 # (assumes rep quality score >= 0)
			# for f in encounteredMaleIndexes:
			# 	score = self.males[f].repQuality() + np.random.normal(scale=self.choiceNoiseSD)
			# 	if score > maxScore:
			# 		fatherIndex = f
			# 		maxScore = score
			fatherIndex = encounteredMaleIndexes[np.argmax(np.add(maleRepQualities[encounteredMaleIndexes],
				np.random.normal(scale=self.choiceNoiseSD,size=self.encounterNum)))]
			father = self.males[fatherIndex]
			if recReprodParents or recParentalArr:
				fatherCounter.update([fatherIndex])

			# Generate the genome of the new member from the parents
			genome = []
			fatherGamete = father.genGamete()
			# if len(fatherGamete[0][0]) > 2:
			# 	print("father gamete" + str(fatherGamete))
			# 	print("father genome" + str(father.genome))
			motherGamete = mother.genGamete()
			# if len(motherGamete[0][0]) > 2:
			# 	print("mother gamete" + str(motherGamete))
			# 	print("mother genome" + str(mother.genome))
			# print("mother gamete" + str(motherGamete))
			for chrom in range(len(fatherGamete)):
				genome += [[fatherGamete[chrom],motherGamete[chrom]]]
			# print(genome)

			# Instantiate the new population member with an expected or random sex
			newIndSex = 'M'
			if numExpectedMales > 0:
				numExpectedMales -= 1
			elif numExpectedFemales > 0:
				newIndSex = 'F'
				numExpectedFemales -= 1
			elif np.random.choice([True,False]):
				newIndSex = 'F'
			newInd = individual(newIndSex,self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
						self.minInvLen,self.lenChrom,self.maleAchiasmy,self.willConvert,self.willRecombine,\
						self.oneInvPerChrom,genome)
			if newIndSex == 'F':
				newFemales += [newInd]
			else:
				newMales += [newInd]

			# # Check if using random sexes
			# if self.randomSex:
			# 	# Instantiate the new population member with a random sex
			# 	if np.random.randint(2):
			# 		newFemales += [individual('F',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
			# 			self.minInvLen,self.lenChrom,self.maleAchiasmy,self.willConvert,self.willRecombine,\
			# 			self.oneInvPerChrom,genome)]
			# 	else:
			# 		newMales += [individual('M',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
			# 			self.minInvLen,self.lenChrom,self.maleAchiasmy,self.willConvert,self.willRecombine,\
			# 			self.oneInvPerChrom,genome)]
			# else:
			# 	# Instantiate the new population member with determinate sex
			# 	# MAY NEED TO ENSURE THAT MOTHERS AREN'T ALWAYS GIVING SONS, ETC
			# 	if numExpectedMales > 0:
			# 		newMales += [individual('M',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
			# 			self.minInvLen,self.lenChrom,self.maleAchiasmy,self.willConvert,self.willRecombine,\
			# 			self.oneInvPerChrom,genome)]
			# 		numExpectedMales -= 1
			# 	else:
			# 		newFemales += [individual('F',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
			# 			self.minInvLen,self.lenChrom,self.maleAchiasmy,self.willConvert,self.willRecombine,\
			# 			self.oneInvPerChrom,genome)]

		if recReprodParents:
			# Record the parental population, and distribution of reproductive success
			inFemales = [self.females[i] for i in obsSurvFemInds]
			inMales = [self.males[i] for i in obsSurvMalInds]
			numObsFemales = len(obsSurvFemInds)
			numObsMales = len(obsSurvMalInds)
			numMothers = len(motherCounter)
			numFathers = len(fatherCounter)
			femRepDist = np.zeros(self.size+1)
			malRepDist = np.zeros(self.size+1)
			for m in motherCounter:
				numOffspring = motherCounter[m]
				femRepDist[numOffspring] += 1
			femRepDist[0] += numObsFemales - numMothers
			for f in fatherCounter:
				numOffspring = fatherCounter[f]
				malRepDist[numOffspring] += 1
			malRepDist[0] += numObsMales - numFathers
			self.__updateRecord(lifeHistStage="Adult",
				inFemaleInds = obsSurvFemInds,
				inMaleInds = obsSurvMalInds,
				motherCounter = motherCounter,
				fatherCounter = fatherCounter)

		# Record the arrangements in the parental population
		if recParentalArr:
			survFemales = [self.females[i] for i in obsSurvFemInds]
			survMales = [self.males[i] for i in obsSurvMalInds]
			# numSurvFemales = len(obsSurvFemInds)
			# numSurvMales = len(obsSurvMalInds)
			self.__recAllArrangementData(lifeHistStage="Adult",
				inFemales = survFemales,
				inMales = survMales)
			mothers = [self.females[i] for i in motherCounter]
			fathers = [self.males[i] for i in fatherCounter]
			# numMothers = len(motherCounter)
			# numFathers = len(fatherCounter)
			self.__recAllArrangementData(lifeHistStage="Parents",
				inFemales = mothers,
				inMales = fathers)


		self.females = newFemales
		self.males = newMales

		# Update the age of the population
		self.age = self.age + 1

		# # Remove fixed inversions and rearrange internal mutations
		# self.__removeFixedInv()
		# Remove fixed mutations and inversions, and rearrange mutations in removed inversions
		if self.removeFixed:
			self.__removeFixed()

		# Add new SA allele mutations
		if self.willMutate:
			# print ("Generation " + str(self.age) + " Mutation")
			# Sprinkle on mutations
			wholePop = self.males + self.females
			# Add new SA mutations
			numMuts = np.random.poisson(self.expectedNumMut)
			if self.verbose and numMuts > 0:
				print("Drawing "+str(numMuts)+" point mutations")
			indivWithMut = np.random.randint(self.size, size=numMuts)
			for i in indivWithMut:
				# Add the mutation to the record, with start and end generations of polymorphism
				self.record[1] += [wholePop[i].mutate(self.__mutIDcount) + [self.age,-1]]
				self.record[2][self.__mutIDcount] = [[] for x in np.arange(18)]
				self.__mutIDcount += 1
				self.__locFixed += [False]
				self.__mutLost += [False]
				# self.record[2] += [[0 for t in range(len(self.record[0]))]]

		# Add new inversion mutations
		if self.willMutInv:
			wholePop = self.males + self.females
			# Add new inversions
			numInvMuts = np.random.poisson(self.expectedNumInvMut)
			if self.verbose and numInvMuts > 0:
				print("Attempting to draw "+str(numInvMuts)+" new inversions")
			addedInvs = []
			addedInvIDs = []
			if self.resampleInvMut:
				possChrom = np.arange(2*self.numChrom*self.size)
				# if self.oneInvPerChrom:
				# 	homWithoutInv = [i for i in possChrom if \
				# 		len(wholePop[i/(2*self.numChrom)].genome[i%(2*self.size)][i%(self.numChrom*self.size)][1]) == 0]

				# arbitrarily limiting failed resamples of inversion locations to 1000 draws
				numFailedDraws = 0
				maxFailedDraws = 1000
				numSuccessfulDraws = 0
				while numFailedDraws < maxFailedDraws and numSuccessfulDraws < numInvMuts:
					randomInd = np.random.randint(self.size)
					invData = wholePop[randomInd].mutateInv(self.__invIDcount)
					if invData is None:
						numFailedDraws += 1
					else:
						numSuccessfulDraws += 1
						addedInvIDs += [self.__invIDcount]
						self.__invIDcount += 1
						addedInvs += [invData]
				if numFailedDraws == maxFailedDraws:
					print("\nWARNING - resampling inversions for new inversion mutants reached the maximum of "+\
						str(maxFailedDraws)+" failed samplings")
			# Add inversions without resampling failures
			else:
				indivWithInvMut = np.random.randint(self.size, size=numInvMuts)
				for i in indivWithInvMut:
					invData = wholePop[i].mutateInv(self.__invIDcount)
					# Must account for failure due to no space on the genome, or conflicting inversion
					if not invData is None:
						addedInvIDs += [self.__invIDcount]
						self.__invIDcount += 1
						addedInvs += [invData]
			for i in np.arange(len(addedInvs)):
				self.record[3] += [addedInvs[i] + [self.age,-1]]
				self.record[4][addedInvIDs[i]] = [[] for x in np.arange(27)]
				# self.__invIDcount += 1
				self.__invFixed += [False]
				self.__invLost += [False]

		# Could update phenotype values here, if precalculated
		return

	# Repeats the last element of the list, useful for fixed mutation and inversion accounting
	# # Has 
	# def __repeat(self,l):
	# 	l += l[len(l)-1]

	# Mode function (so that scipy isn't needed)
	# returns a single value which is either the mode (if there is only one most frequenct),
	#   or a random sample from among the equally-most-frequent classes
	def __randomizedMode(self, data):
		# from collections import Counter
		c = Counter(data)
		highFreqVals = [k for (k, v) in c.items() if v == c.most_common(1)[0][1]]
		return highFreqVals[np.random.choice(len(highFreqVals))]

	# Takes a haplotype and a haplotype definition (as lists of mutation and inversion IDs),
	#   and returns an index to store the haplotype in within the population record
	def __indexFromHap(self,hap,hapDef):
		indexSum = 0
		for mutID in hap[0]:
			indexSum += 2**(hapDef[0].index(mutID))
		numDefMuts = len(hapDef[0])
		for invID in hap[1]:
			indexSum += 2**(numDefMuts+(hapDef[1].index(invID)))
		return indexSum

	# For collecting record info on a pair of chromosomes from an individual
	def __countIndChromMutInvGenotypes(self,chrom,sex,mutCounters,invCounters,hapCounters,invDataDict,
					mutOffCounters,invOffCounters,hapOffCounters,indivOffCount):
		chromSur = 1
		chromRep = 0

		assert((mutOffCounters is not None and invOffCounters is not None and hapOffCounters is not None and indivOffCount is not None) or \
			(mutOffCounters is None and invOffCounters is None and hapOffCounters is None and indivOffCount is None)), \
			"In counting variants in an individual for record keeping, either all offspring count inputs must be given or all set to None."
		recOffCounts = False
		if mutOffCounters is not None:
			recOffCounts = True

		# Account for haplotype recording. If no haplotypes are given, it should never run
		if self.recordHaps:
			numHaps = len(self.record[11])
			hapRange = np.arange(numHaps)
			hapsHom1 = [[[],[]] for i in hapRange]
			hapsHom2 = [[[],[]] for i in hapRange]
			for h in hapRange:
				for m in self.record[11][h][0]:
					if self.__locFixed[m] and not self.__mutLost[m]:
						hapsHom1[h][0] += [m]
						hapsHom2[h][0] += [m]
				for i in self.record[11][h][1]:
					if self.__invFixed[i] and not self.__invLost[i]:
						hapsHom1[h][1] += [i]
						hapsHom2[h][1] += [i]
		else:
			hapRange = []


		# Figure out which inversion buffers to check, and account for inversions and their zygosity
		numInvHom1 = len(chrom[0][1])
		numInvHom2 = len(chrom[1][1])
		if numInvHom1 == 0:
			invCounters[0][0].update(["Std"])
			if sex == 'F':
				invCounters[1][0].update(["Std"])
			else:
				invCounters[2][0].update(["Std"])
			if recOffCounts:
				invOffCounters[0][0].update(["Std"]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][0].update(["Std"]*indivOffCount)
				else:
					invOffCounters[2][0].update(["Std"]*indivOffCount)
			invIDsHom1 = ["Std"]
			invBuffersHom1 = [[0-self.invRecBuffer,self.lenChrom+self.invRecBuffer]]
			invsHom1numMutInside = [0]
			invsHom1SurvEffect = [1]
			invsHom1ReprEffect = [0]
		else:
			invIDsHom1 = []
			invBuffersHom1 = []
			invsHom1numMutInside = [0]*numInvHom1
			invsHom1SurvEffect = [1]*numInvHom1
			invsHom1ReprEffect = [0]*numInvHom1
		if numInvHom2 == 0:
			invCounters[0][0].update(["Std"])
			if sex == 'F':
				invCounters[1][0].update(["Std"])
			else:
				invCounters[2][0].update(["Std"])
			if recOffCounts:
				invOffCounters[0][0].update(["Std"]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][0].update(["Std"]*indivOffCount)
				else:
					invOffCounters[2][0].update(["Std"]*indivOffCount)
			invIDsHom2 = ["Std"]
			invBuffersHom2 = [[0-self.invRecBuffer,self.lenChrom+self.invRecBuffer]]
			invsHom2numMutInside = [0]
			invsHom2SurvEffect = [1]
			invsHom2ReprEffect = [0]
		else:
			invIDsHom2 = []
			invBuffersHom2 = []
			invsHom2numMutInside = [0]*numInvHom2
			invsHom2SurvEffect = [1]*numInvHom2
			invsHom2ReprEffect = [0]*numInvHom2

		if numInvHom1 == 0 and numInvHom2 == 0:
			invCounters[0][1].update(["Std"])
			if sex == 'F':
				invCounters[1][1].update(["Std"])
			else:
				invCounters[2][1].update(["Std"])
		elif numInvHom1 == 0 or numInvHom2 == 0:
			invCounters[0][2].update(["Std"])
			if sex == 'F':
				invCounters[1][2].update(["Std"])
			else:
				invCounters[2][2].update(["Std"])

		if recOffCounts:
			if numInvHom1 == 0 and numInvHom2 == 0:
				invOffCounters[0][1].update(["Std"]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][1].update(["Std"]*indivOffCount)
				else:
					invOffCounters[2][1].update(["Std"]*indivOffCount)
			elif numInvHom1 == 0 or numInvHom2 == 0:
				invOffCounters[0][2].update(["Std"]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][2].update(["Std"]*indivOffCount)
				else:
					invOffCounters[2][2].update(["Std"]*indivOffCount)


		homInd1 = 0
		homInd2 = 0
		# Check for heterozygosity and account for inversions in the record
		# !!!! Relies on the invariant that earlier indexed inversions have earlier position
		while homInd1 < numInvHom1 and homInd2 < numInvHom2:
			inv1 = chrom[0][1][homInd1]
			inv2 = chrom[1][1][homInd2]
			if inv1[0] == inv2[0]:
				assert(str(inv1)==str(inv2)), "two inversions should only have the same first breakpoint if they are the same inversion"
				homInd1 += 1
				homInd2 += 1
				currID = inv1[2]
				# invCounts.update([currID,currID])
				# invHomCounts.update([currID])
				invCounters[0][0].update([currID,currID])
				invCounters[0][1].update([currID])
				if sex == 'F':
					invCounters[1][0].update([currID,currID])
					invCounters[1][1].update([currID])
				else:
					invCounters[2][0].update([currID,currID])
					invCounters[2][1].update([currID])
				if recOffCounts:
					invOffCounters[0][0].update([currID,currID]*indivOffCount)
					invOffCounters[0][1].update([currID]*indivOffCount)
					if sex == 'F':
						invOffCounters[1][0].update([currID,currID]*indivOffCount)
						invOffCounters[1][1].update([currID]*indivOffCount)
					else:
						invOffCounters[2][0].update([currID,currID]*indivOffCount)
						invOffCounters[2][1].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][1]:
						hapsHom1[hapInd][1] += [currID]
						hapsHom2[hapInd][1] += [currID]
				invIDsHom1 += [currID]
				invBuffersHom1 += [[inv1[0]-self.invRecBuffer,inv1[1]+self.invRecBuffer]]
				invIDsHom2 += [currID]
				invBuffersHom2 += [[inv1[0]-self.invRecBuffer,inv1[1]+self.invRecBuffer]]
				

			if inv1[0] < inv2[0]:
				homInd1 += 1
				currID = inv1[2]
				# invCounts.update([currID])
				# invHetCounts.update([currID])
				invCounters[0][0].update([currID])
				invCounters[0][2].update([currID])
				if sex == 'F':
					invCounters[1][0].update([currID])
					invCounters[1][2].update([currID])
				else:
					invCounters[2][0].update([currID])
					invCounters[2][2].update([currID])
				if recOffCounts:
					invOffCounters[0][0].update([currID]*indivOffCount)
					invOffCounters[0][2].update([currID]*indivOffCount)
					if sex == 'F':
						invOffCounters[1][0].update([currID]*indivOffCount)
						invOffCounters[1][2].update([currID]*indivOffCount)
					else:
						invOffCounters[2][0].update([currID]*indivOffCount)
						invOffCounters[2][2].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][1]:
						hapsHom1[hapInd][1] += [currID]
				invIDsHom1 += [currID]
				invBuffersHom1 += [[inv1[0]-self.invRecBuffer,inv1[1]+self.invRecBuffer]]

			if inv1[0] > inv2[0]:
				homInd2 += 1
				currID = inv2[2]
				# invCounts.update([currID])
				# invHetCounts.update([currID])
				invCounters[0][0].update([currID])
				invCounters[0][2].update([currID])
				if sex == 'F':
					invCounters[1][0].update([currID])
					invCounters[1][2].update([currID])
				else:
					invCounters[2][0].update([currID])
					invCounters[2][2].update([currID])
				if recOffCounts:
					invOffCounters[0][0].update([currID]*indivOffCount)
					invOffCounters[0][2].update([currID]*indivOffCount)
					if sex == 'F':
						invOffCounters[1][0].update([currID]*indivOffCount)
						invOffCounters[1][2].update([currID]*indivOffCount)
					else:
						invOffCounters[2][0].update([currID]*indivOffCount)
						invOffCounters[2][2].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][1]:
						hapsHom2[hapInd][1] += [currID]
				invIDsHom2 += [currID]
				invBuffersHom2 += [[inv2[0]-self.invRecBuffer,inv2[1]+self.invRecBuffer]]

		# Deal with the remaining inversions
		# when one chromosome has no inversions it is dealt with in one of the following two code blocks
		while homInd1 < numInvHom1:
			inv1 = chrom[0][1][homInd1]
			homInd1 += 1
			currID = inv1[2]
			# invCounts.update([currID])
			# invHetCounts.update([currID])
			invCounters[0][0].update([currID])
			invCounters[0][2].update([currID])
			if sex == 'F':
				invCounters[1][0].update([currID])
				invCounters[1][2].update([currID])
			else:
				invCounters[2][0].update([currID])
				invCounters[2][2].update([currID])
			if recOffCounts:
				invOffCounters[0][0].update([currID]*indivOffCount)
				invOffCounters[0][2].update([currID]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][0].update([currID]*indivOffCount)
					invOffCounters[1][2].update([currID]*indivOffCount)
				else:
					invOffCounters[2][0].update([currID]*indivOffCount)
					invOffCounters[2][2].update([currID]*indivOffCount)
			for hapInd in hapRange:
				if currID in self.record[11][hapInd][1]:
					hapsHom1[hapInd][1] += [currID]
			invIDsHom1 += [currID]
			invBuffersHom1 += [[inv1[0]-self.invRecBuffer,inv1[1]+self.invRecBuffer]]

		while homInd2 < numInvHom2:
			inv2 = chrom[1][1][homInd2]
			homInd2 += 1
			currID = inv2[2]
			# invCounts.update([currID])
			# invHetCounts.update([currID])
			invCounters[0][0].update([currID])
			invCounters[0][2].update([currID])
			if sex == 'F':
				invCounters[1][0].update([currID])
				invCounters[1][2].update([currID])
			else:
				invCounters[2][0].update([currID])
				invCounters[2][2].update([currID])
			if recOffCounts:
				invOffCounters[0][0].update([currID]*indivOffCount)
				invOffCounters[0][2].update([currID]*indivOffCount)
				if sex == 'F':
					invOffCounters[1][0].update([currID]*indivOffCount)
					invOffCounters[1][2].update([currID]*indivOffCount)
				else:
					invOffCounters[2][0].update([currID]*indivOffCount)
					invOffCounters[2][2].update([currID]*indivOffCount)
			for hapInd in hapRange:
				if currID in self.record[11][hapInd][1]:
					hapsHom2[hapInd][1] += [currID]
			invIDsHom2 += [currID]
			invBuffersHom2 += [[inv2[0]-self.invRecBuffer,inv2[1]+self.invRecBuffer]]


		# Account the mutations, and add inversion accounting details if they're within inversion buffers
		invHom1Range = np.arange(len(invIDsHom1))
		invHom2Range = np.arange(len(invIDsHom2))
		homInd1 = 0
		homInd2 = 0
		# Check for heterozygosity and account for mutations in the record
		# !!!! Relies on the invariant that earlier indexed mutations have earlier position
		# Account for multiplicative survival effects - average effect is average of multiplied effects
		while homInd1 < len(chrom[0][0]) and homInd2 < len(chrom[1][0]):
			mut1 = chrom[0][0][homInd1]
			mut2 = chrom[1][0][homInd2]
			if mut1[0] == mut2[0]:
				assert(str(mut1)==str(mut2)), "Two mutations should only be in the same place if they are the same"
				homInd1 += 1
				homInd2 += 1
				currID = mut1[3]
				chromSur *= mut1[1]**2
				chromRep += mut1[2]*2
				# mutCounts.update([currID,currID])
				# mutHomCounts.update([currID])
				mutCounters[0][0].update([currID,currID])
				mutCounters[0][1].update([currID])
				if sex == 'F':
					mutCounters[1][0].update([currID,currID])
					mutCounters[1][1].update([currID])
				else:
					mutCounters[2][0].update([currID,currID])
					mutCounters[2][1].update([currID])
				if recOffCounts:
					mutOffCounters[0][0].update([currID,currID]*indivOffCount)
					mutOffCounters[0][1].update([currID]*indivOffCount)
					if sex == 'F':
						mutOffCounters[1][0].update([currID,currID]*indivOffCount)
						mutOffCounters[1][1].update([currID]*indivOffCount)
					else:
						mutOffCounters[2][0].update([currID,currID]*indivOffCount)
						mutOffCounters[2][1].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][0]:
						hapsHom1[hapInd][0] += [currID]
						hapsHom2[hapInd][0] += [currID]
				for i in invHom1Range:
					if mut1[0] > invBuffersHom1[i][0] and mut1[0] < invBuffersHom1[i][1]:
						invsHom1numMutInside[i] += 1
						invsHom1SurvEffect[i] *= mut1[1]
						invsHom1ReprEffect[i] += mut1[2]
				for i in invHom2Range:
					if mut2[0] > invBuffersHom2[i][0] and mut2[0] < invBuffersHom2[i][1]:
						invsHom2numMutInside[i] += 1
						invsHom2SurvEffect[i] *= mut2[1]
						invsHom2ReprEffect[i] += mut2[2]

			if mut1[0] < mut2[0]:
				homInd1 += 1
				currID = mut1[3]
				chromSur *= mut1[1]
				chromRep += mut1[2]
				# mutCounts.update([currID])
				# mutHetCounts.update([currID])
				mutCounters[0][0].update([currID])
				mutCounters[0][2].update([currID])
				if sex == 'F':
					mutCounters[1][0].update([currID])
					mutCounters[1][2].update([currID])
				else:
					mutCounters[2][0].update([currID])
					mutCounters[2][2].update([currID])
				if recOffCounts:
					mutOffCounters[0][0].update([currID]*indivOffCount)
					mutOffCounters[0][2].update([currID]*indivOffCount)
					if sex == 'F':
						mutOffCounters[1][0].update([currID]*indivOffCount)
						mutOffCounters[1][2].update([currID]*indivOffCount)
					else:
						mutOffCounters[2][0].update([currID]*indivOffCount)
						mutOffCounters[2][2].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][0]:
						hapsHom1[hapInd][0] += [currID]
				for i in invHom1Range:
					if mut1[0] > invBuffersHom1[i][0] and mut1[0] < invBuffersHom1[i][1]:
						invsHom1numMutInside[i] += 1
						invsHom1SurvEffect[i] *= mut1[1]
						invsHom1ReprEffect[i] += mut1[2]

			if mut1[0] > mut2[0]:
				homInd2 += 1
				currID = mut2[3]
				chromSur *= mut2[1]
				chromRep += mut2[2]
				# mutCounts.update([currID])
				# mutHetCounts.update([currID])
				mutCounters[0][0].update([currID])
				mutCounters[0][2].update([currID])
				if sex == 'F':
					mutCounters[1][0].update([currID])
					mutCounters[1][2].update([currID])
				else:
					mutCounters[2][0].update([currID])
					mutCounters[2][2].update([currID])
				if recOffCounts:
					mutOffCounters[0][0].update([currID]*indivOffCount)
					mutOffCounters[0][2].update([currID]*indivOffCount)
					if sex == 'F':
						mutOffCounters[1][0].update([currID]*indivOffCount)
						mutOffCounters[1][2].update([currID]*indivOffCount)
					else:
						mutOffCounters[2][0].update([currID]*indivOffCount)
						mutOffCounters[2][2].update([currID]*indivOffCount)
				for hapInd in hapRange:
					if currID in self.record[11][hapInd][0]:
						hapsHom2[hapInd][0] += [currID]
				for i in invHom2Range:
					if mut2[0] > invBuffersHom2[i][0] and mut2[0] < invBuffersHom2[i][1]:
						invsHom2numMutInside[i] += 1
						invsHom2SurvEffect[i] *= mut2[1]
						invsHom2ReprEffect[i] += mut2[2]

		# Deal with the remaining mutations
		# when one chromosome has no mutations it is dealt with here immediately
		while homInd1 < len(chrom[0][0]):
			mut1 = chrom[0][0][homInd1]
			homInd1 += 1
			currID = mut1[3]
			chromSur *= mut1[1]
			chromRep += mut1[2]
			# mutCounts.update([currID])
			# mutHetCounts.update([currID])
			mutCounters[0][0].update([currID])
			mutCounters[0][2].update([currID])
			if sex == 'F':
				mutCounters[1][0].update([currID])
				mutCounters[1][2].update([currID])
			else:
				mutCounters[2][0].update([currID])
				mutCounters[2][2].update([currID])
			if recOffCounts:
				mutOffCounters[0][0].update([currID]*indivOffCount)
				mutOffCounters[0][2].update([currID]*indivOffCount)
				if sex == 'F':
					mutOffCounters[1][0].update([currID]*indivOffCount)
					mutOffCounters[1][2].update([currID]*indivOffCount)
				else:
					mutOffCounters[2][0].update([currID]*indivOffCount)
					mutOffCounters[2][2].update([currID]*indivOffCount)
			for hapInd in hapRange:
				if currID in self.record[11][hapInd][0]:
					hapsHom1[hapInd][0] += [currID]
			for i in invHom1Range:
				if mut1[0] > invBuffersHom1[i][0] and mut1[0] < invBuffersHom1[i][1]:
					invsHom1numMutInside[i] += 1
					invsHom1SurvEffect[i] *= mut1[1]
					invsHom1ReprEffect[i] += mut1[2]

		while homInd2 < len(chrom[1][0]):
			mut2 = chrom[1][0][homInd2]
			homInd2 += 1
			currID = mut2[3]
			chromSur *= mut2[1]
			chromRep += mut2[2]
			# mutCounts.update([currID])
			# mutHetCounts.update([currID])
			mutCounters[0][0].update([currID])
			mutCounters[0][2].update([currID])
			if sex == 'F':
				mutCounters[1][0].update([currID])
				mutCounters[1][2].update([currID])
			else:
				mutCounters[2][0].update([currID])
				mutCounters[2][2].update([currID])
			if recOffCounts:
				mutOffCounters[0][0].update([currID]*indivOffCount)
				mutOffCounters[0][2].update([currID]*indivOffCount)
				if sex == 'F':
					mutOffCounters[1][0].update([currID]*indivOffCount)
					mutOffCounters[1][2].update([currID]*indivOffCount)
				else:
					mutOffCounters[2][0].update([currID]*indivOffCount)
					mutOffCounters[2][2].update([currID]*indivOffCount)
			for hapInd in hapRange:
				if currID in self.record[11][hapInd][0]:
					hapsHom2[hapInd][0] += [currID]
			for i in invHom2Range:
				if mut2[0] > invBuffersHom2[i][0] and mut2[0] < invBuffersHom2[i][1]:
					invsHom2numMutInside[i] += 1
					invsHom2SurvEffect[i] *= mut2[1]
					invsHom2ReprEffect[i] += mut2[2]


		# Record the collected inversion-relevant data
		for i in invHom1Range:
			if invIDsHom1[i] not in invDataDict:
				invDataDict[invIDsHom1[i]] = [[],[],[]]
			invDataDict[invIDsHom1[i]][0] += [invsHom1numMutInside[i]]
			invDataDict[invIDsHom1[i]][1] += [invsHom1SurvEffect[i]]
			invDataDict[invIDsHom1[i]][2] += [invsHom1ReprEffect[i]]
		for i in invHom2Range:
			if invIDsHom2[i] not in invDataDict:
				invDataDict[invIDsHom2[i]] = [[],[],[]]
			invDataDict[invIDsHom2[i]][0] += [invsHom2numMutInside[i]]
			invDataDict[invIDsHom2[i]][1] += [invsHom2SurvEffect[i]]
			invDataDict[invIDsHom2[i]][2] += [invsHom2ReprEffect[i]]


		# Record the collected haplotype data
		for h in hapRange:
			# hap1 = HashableList(hapsHom1[h])
			# hap2 = HashableList(hapsHom2[h])
			# hapCountsList[h].update([hap1,hap2])
			# if hap1 == hap2:
			# 	hapHomCountsList[h].update([hap1])
			# else:
			# 	hapHetCountsList[h].update([hap1,hap2])
			hap1i = self.__indexFromHap(hapsHom1[h],self.record[11][h])
			hap2i = self.__indexFromHap(hapsHom2[h],self.record[11][h])
			hapCounters[h][0][0].update([hap1i,hap2i])
			if hap1i == hap2i:
				hapCounters[h][0][1].update([hap1i])
			else:
				hapCounters[h][0][2].update([hap1i,hap2i])
			if sex == 'F':
				hapCounters[h][1][0].update([hap1i,hap2i])
				if hap1i == hap2i:
					hapCounters[h][1][1].update([hap1i])
				else:
					hapCounters[h][1][2].update([hap1i,hap2i])
			else:
				hapCounters[h][2][0].update([hap1i,hap2i])
				if hap1i == hap2i:
					hapCounters[h][2][1].update([hap1i])
				else:
					hapCounters[h][2][2].update([hap1i,hap2i])

			if recOffCounts:
				hapOffCounters[h][0][0].update([hap1i,hap2i]*indivOffCount)
				if hap1i == hap2i:
					hapOffCounters[h][0][1].update([hap1i]*indivOffCount)
				else:
					hapOffCounters[h][0][2].update([hap1i,hap2i]*indivOffCount)
				if sex == 'F':
					hapOffCounters[h][1][0].update([hap1i,hap2i]*indivOffCount)
					if hap1i == hap2i:
						hapOffCounters[h][1][1].update([hap1i]*indivOffCount)
					else:
						hapOffCounters[h][1][2].update([hap1i,hap2i]*indivOffCount)
				else:
					hapOffCounters[h][2][0].update([hap1i,hap2i]*indivOffCount)
					if hap1i == hap2i:
						hapOffCounters[h][2][1].update([hap1i]*indivOffCount)
					else:
						hapOffCounters[h][2][2].update([hap1i,hap2i]*indivOffCount)

		return (chromSur,chromRep)

	# A helper function to ensure the same haplotype class range is used in each call
	def __calcHapClassRange(self,h):
		hapClassesRange = np.arange(2**(len(self.record[11][h][0])+len(self.record[11][h][1])))
		return hapClassesRange

	# For updating all recorded information on the population, record structured as:
	# [0] a list of ages at which an update was made
	# [1] a list of all mutations, indexed by ID, with each entry as 
	#      [position,survival effect,reproductive effect,chromosome,gen at which mutation occured,gen at fixation/loss]
	#      -1 means the mutation is not fixed/lost yet
	# [2] a dictionary of mutation count data, indexed by ID, of the mutation at each gen noted in [0]
	#        KEY/LIST-DATA PAIRS ARE REMOVED TO SAVE DATA 
	#         IF THE MUTATIONS SURVIVES FOR LESS THAN __ GENERATIONS AND REACHES LESS THAT __ FREQUENCY
	#      [0] counts of the mutation at each gen noted in [0]
	#      [1] count of the number of individuals homozygous for the mutation (implicitly the derived allele)
	#      [2] count of the number of individuals heterozygous for the mutation
	#      [3] count of the number of mutations in females
	#      [4] count of the number of females homozygous for the mutation (implicitly the derived allele)
	#      [5] count of the number of females heterozygous for the mutation
	#      [6] count of the number of mutations in males
	#      [7] count of the number of males homozygous for the mutation (implicitly the derived allele)
	#      [8] count of the number of males heterozygous for the mutation
	#      [9] counts of allelic descendants of the mutation across all individuals (1x homozygous offspring, 1/2 heterozygous)
	#      [10] count of offspring of individuals homozygous for the mutation
	#      [11] count of offspring of individuals heterozygous for the mutation
	#      [12] counts of allelic descendants of the mutation across females
	#      [13] count of offspring of females homozygous for the mutation
	#      [14] count of offspring of females heterozygous for the mutation
	#      [15] counts of allelic descendants of the mutation across males
	#      [16] count of offspring of males homozygous for the mutation
	#      [17] count of offspring of males heterozygous for the mutation
	# [3] a list of all inversions, indexed by ID, with each entry as 
	#      [start position,end position,chromosome,gen at which mutation occured,gen at fixation/loss]
	#      -1 means the mutation is not fixed/lost yet
	# [4] a dictionary of collected inversion data, keyed by ID:
	#        KEY/LIST-DATA PAIRS ARE REMOVED TO SAVE DATA 
	#         IF THE INVERSION SURVIVES FOR LESS THAN __ GENERATIONS AND REACHES LESS THAT __ FREQUENCY
	#      [0] counts of the inversion at each gen noted in [0]
	#      [1] average count of mutations within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [2] variance in count of mutations within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [3] mode of the count of mutations within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [4] average total survival effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [5] variance in total survival effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [6] mode total survival effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [7] average total reproductive effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [8] variance in total reproductive effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [9] mode total reproductive effect within the buffer region of the inversion 
	#        across the entire population, at each gen noted in [0]
	#      [10] count of the number of individuals homozygous for the inversion (implicitly the derived arrangement)
	#      [11] count of the number of individuals heterozygous for the inversion
	#      [12] count of the number of inversions in females
	#      [13] count of the number of females homozygous for the inversion (implicitly the derived arrangement)
	#      [14] count of the number of females heterozygous for the inversion
	#      [15] count of the number of inversions in males
	#      [16] count of the number of males homozygous for the inversion (implicitly the derived arrangement)
	#      [17] count of the number of males heterozygous for the inversion
	#      [18] counts of allelic descendants of the inversion across all individuals (1x homozygous offspring, 1/2 heterozygous)
	#      [19] count of offspring of individuals homozygous for the inversion
	#      [20] count of offspring of individuals heterozygous for the inversion
	#      [21] counts of allelic descendants of the inversion across females
	#      [22] count of offspring of females homozygous for the inversion
	#      [23] count of offspring of females heterozygous for the inversion
	#      [24] counts of allelic descendants of the inversion across males
	#      [25] count of offspring of males homozygous for the inversion
	#      [26] count of offspring of males heterozygous for the inversion
	# [5] a list of pop size at each recorded generation
	# [6] mean population survival value at each recorded generation
	# [7] variance in population survival value at each recorded generation
	# [8] mean population reproductive value at each recorded generation
	# [9] variance in population reproductive value at each recorded generation
	# [10] the life history stage of each record, at the associated generation
	# [11] a list of all haplotypes to record, indexed by ID, with each entry as 
	#       [list of mutation IDs, list of inversion IDs, generation recording started] ##, index of first record with the haplotype]
	# [12] a dictionary collecting haplotype counts, keyed by the haplotype ID
	#      [0] counts of each haplotype option defined by the haplotype definition, ordered as
	#            [0,0,...0],[0]; 
	#            [1,0,...0],[0]; 
	#            [0,1,...0],[0]; 
	#            [1,1,...0],[0]; 
	#             ...
	#            [1,1,...1],[0]; 
	#            [0,0,...0],[1]; 
	#             ...
	#            [1,1,...0],[1]; 
	#            [1,1,...1],[1]
	#           This order is accessed using the indexFromHap function 
	#      [1] counts of homozygotes for each haplotype option, ordered as above
	#      [2] counts of heterozygotes for each haplotype option, ordered as above
	#      [3] counts of each haplotype option in females, ordered as above
	#      [4] counts of homozygotes for each haplotype option in females, ordered as above
	#      [5] counts of heterozygotes for each haplotype option in females, ordered as above
	#      [6] counts of each haplotype option in males, ordered as above
	#      [7] counts of homozygotes for each haplotype option in males, ordered as above
	#      [8] counts of heterozygotes for each haplotype option in males, ordered as above
	#      [9] counts of allelic descendants of each haplotype across all individuals (1x homozygous offspring, 1/2 heterozygous)
	#      [10] count of offspring of individuals homozygous for each haplotype option, ordered as above
	#      [11] count of offspring of individuals heterozygous for each haplotype option, ordered as above
	#      [12] counts of allelic descendants of each haplotype across females
	#      [13] count of offspring of each homozygous haplotype option in females, ordered as above
	#      [14] count of offspring of each heterozygous haplotype option in females, ordered as above
	#      [15] counts of allelic descendants of each haplotype across males
	#      [16] count of offspring of each homozygous haplotype option in males, ordered as above
	#      [17] count of offspring of each heterozygous haplotype option in males, ordered as above
	# [13] a list of number of female individuals at each record
	# [14] a list of number of male individuals at each record
	# [15] a distribution of number of offspring per mother at each recording of Adults, 'None' otherwise
	# [16] a distribution of number of offspring per father at each recording of Adults, 'None' otherwise
	def __updateRecord(self,lifeHistStage="Zygote",inFemaleInds=None,inMaleInds=None,
			motherCounter=None,fatherCounter=None):
		# print ("Record Update")
		# print ("curr record: "+str(self.record))

		assert(lifeHistStage == "Zygote" or (inFemaleInds is not None and inMaleInds is not None)), "Life history stages "+\
			"other than zygote must have an input of female and male individuals"
		assert(lifeHistStage != "Zygote" or (self.size == len(self.males)+len(self.females))), "At zygote life history stage, "+\
			"expects self.males and self.females to total the set population size"
		assert((motherCounter is not None and fatherCounter is not None) or (motherCounter is None and fatherCounter is None)), \
			"If offspring counts are given, must be given for both mother and father counters."
		recOffCounts = motherCounter is not None
		# Record the age of the record update
		currRecIndex = len(self.record[0])
		self.record[0] += [self.age]
		# Record the population size at the time of the of the record update, assemble the relevant individuals
		if lifeHistStage != "Zygote":
			# Calculate the distribution of reproductive success
			numFemales = len(inFemaleInds)
			numMales = len(inMaleInds)
			currPopSize = len(inFemaleInds)+len(inMaleInds)
			inFemales = [self.females[i] for i in inFemaleInds]
			inMales = [self.males[i] for i in inMaleInds]
			currIndividuals = inFemales+inMales
			inFemaleOffCounts = [motherCounter[i] for i in inFemaleInds]
			inMaleOffCounts = [fatherCounter[i] for i in inMaleInds]
			currIndivOffCounts = inFemaleOffCounts + inMaleOffCounts
			numMothers = len(motherCounter)
			numFathers = len(fatherCounter)
			femRepDist = np.zeros(self.size+1)
			malRepDist = np.zeros(self.size+1)
			for m in inFemaleInds:
				numOffspring = motherCounter[m]
				femRepDist[numOffspring] += 1
			# femRepDist[0] += numFemales - numMothers
			for f in inMaleInds:
				numOffspring = fatherCounter[f]
				malRepDist[numOffspring] += 1
			# malRepDist[0] += numMales - numFathers
		else:
			numFemales = len(self.females)
			numMales = len(self.males)
			currPopSize = self.size
			currIndividuals = self.females+self.males
			femRepDist = None
			malRepDist = None
			currIndivOffCounts = None
		self.record[5] += [currPopSize]
		# Record the life history stage
		self.record[10] += [lifeHistStage]
		self.record[13] += [numFemales]
		self.record[14] += [numMales]
		self.record[15] += [femRepDist]
		self.record[16] += [malRepDist]

		# Count the number of each mutation/inversion by ID

		# mutCounts = Counter()
		# mutHomCounts = Counter()
		# mutHetCounts = Counter()
		# invCounts = Counter()
		# invHomCounts = Counter()
		# invHetCounts = Counter()
		# hapCountsList = []
		# hapHomCountsList = []
		# hapHetCountsList = []
		# hapRange = np.arange(len(self.record[11]))
		# for h in hapRange:
		# 	hapCountsList += [Counter()]
		# 	hapHomCountsList += [Counter()]
		# 	hapHetCountsList += [Counter()]
		countCategories=["Allele","Homozygote","Heterozygote"]
		sexCategories=["All","Female","Male"]
		countRange = np.arange(len(countCategories))
		sexRange = np.arange(len(sexCategories))
		hapRange = np.arange(len(self.record[11]))
		mutCounters = [[Counter() for x in countRange] for y in sexRange]
		invCounters = [[Counter() for x in countRange] for y in sexRange]
		hapCounters = [[[Counter() for x in countRange] for y in sexRange] for h in hapRange]
		if (motherCounter is not None and fatherCounter is not None):
			mutOffCounters = [[Counter() for x in countRange] for y in sexRange]
			invOffCounters = [[Counter() for x in countRange] for y in sexRange]
			hapOffCounters = [[[Counter() for x in countRange] for y in sexRange] for h in hapRange]
		else:
			mutOffCounters = None
			invOffCounters = None
			hapOffCounters = None
			# mutOffCounters = [[None for x in countRange] for y in sexRange]
			# invOffCounters = [[None for x in countRange] for y in sexRange]
			# hapOffCounters = [[[None for x in countRange] for y in sexRange] for h in hapRange]


		# Generate the metrics of value associations for the inversions
		# numMutInBuffer = [0]*self.__invIDcount
		# survEffectInBuffer = [0]*self.__invIDcount
		# reprEffectInBuffer = [0]*self.__invIDcount
		# numMutsInBuffer = {}
		# survEffectsInBuffer = {}
		# reprEffectsInBuffer = {}
		invDataDict = {}

		# Account whole population statistics (mean and variance of values/individual)
		popSurvVals = np.zeros(currPopSize)
		popReprVals = np.zeros(currPopSize)
		# i = 0
		# # Manage Haplotype Accounting
		# mutsInHaps = set()
		# invsInHaps = set()
		# if recordHaps:
		# 	for haplotype in self.record[11]:
		# 		mutsInHaps.update(set(haplotype[0]))
		# 		invsInHaps.update(set(haplotype[1]))

		for i in np.arange(currPopSize):
			indiv = currIndividuals[i]
			indivSur = 1
			indivRep = 0
			indivOffCount = None
			if currIndivOffCounts is not None:
				indivOffCount = currIndivOffCounts[i]
			for chrom in indiv.genome:
				(chromSur,chromRep) = self.__countIndChromMutInvGenotypes(chrom,indiv.sex,
					mutCounters,invCounters,hapCounters,invDataDict,
					mutOffCounters,invOffCounters,hapOffCounters,indivOffCount)
				indivSur *= chromSur
				indivRep += chromRep
			popSurvVals[i] = indivSur
			popReprVals[i] = indivRep
			# i += 1
		# Update the record, get averages for record[5,6,7]
		# for i in range(self.__mutIDcount):
		# 	count = mutCounts[i]
		# 	if count == 0:
		# 		if self.__locFixed[i]:  # Also check mutLost?
		# 			self.record[2][i] += [2*self.size]
		# 		else:
		# 			self.record[2][i] += [0]
		# 	else:
		# 		self.record[2][i] += [count]
		# for j in range(self.__invIDcount):
		# 	count = invCounts[j]
		# 	# If no members of the inversion ID left, record averages as -1 ? Try 'NA' for working with R?
		# 	# IMPORTANT - add a boolean list such that you can check if fixed and record the correct value?
		# 	if count == 0:
		# 		if self.__invFixed
		# 			self.record[4]+= [2*self.size]
		# 		else:
		# 			self.record[4][j] += [0]
		# 		self.record[5][j] += [-1]
		# 		self.record[6][j] += [-1]
		# 		self.record[7][j] += [-1]
		# 	else:
		# 		self.record[4][j] += [count]
		# 		self.record[5][j] += [numMutInBuffer[j]/(float(count))]
		# 		self.record[6][j] += [survEffectTotalInBufferMultiplicative[j]/(float(count))]
		# 		self.record[7][j] += [reprEffectTotalInBuffer[j]/(float(count))]
		# print(mutCounts)
		# print(invCounts)
		for i in mutCounters[0][0].keys():
			count = mutCounters[0][0][i]
			if count != 0:
				if i not in self.record[2]:
					self.record[2][i] = [[] for x in np.arange(18)]
				self.record[2][i][0] += [count]
				self.record[2][i][1] += [mutCounters[0][1][i]]
				self.record[2][i][2] += [mutCounters[0][2][i]]
				self.record[2][i][3] += [mutCounters[1][0][i]]
				self.record[2][i][4] += [mutCounters[1][1][i]]
				self.record[2][i][5] += [mutCounters[1][2][i]]
				self.record[2][i][6] += [mutCounters[2][0][i]]
				self.record[2][i][7] += [mutCounters[2][1][i]]
				self.record[2][i][8] += [mutCounters[2][2][i]]
				if recOffCounts:
					self.record[2][i][9] += [mutOffCounters[0][0][i]/2]
					self.record[2][i][10] += [mutOffCounters[0][1][i]]
					self.record[2][i][11] += [mutOffCounters[0][2][i]]
					self.record[2][i][12] += [mutOffCounters[1][0][i]/2]
					self.record[2][i][13] += [mutOffCounters[1][1][i]]
					self.record[2][i][14] += [mutOffCounters[1][2][i]]
					self.record[2][i][15] += [mutOffCounters[2][0][i]/2]
					self.record[2][i][16] += [mutOffCounters[2][1][i]]
					self.record[2][i][17] += [mutOffCounters[2][2][i]]
				else:
					self.record[2][i][9] += [None]
					self.record[2][i][10] += [None]
					self.record[2][i][11] += [None]
					self.record[2][i][12] += [None]
					self.record[2][i][13] += [None]
					self.record[2][i][14] += [None]
					self.record[2][i][15] += [None]
					self.record[2][i][16] += [None]
					self.record[2][i][17] += [None]
			else:
				raise self.SimulationStateError((i,count),"Recording had a mutation counter "+\
					"ID/key with value 0, should never happen with a counter - Mutation ID "+str(i))
		if "Std" not in invCounters[0][0].keys():
			if "Std" not in self.record[4]:
				self.record[4]["Std"] = [[] for x in np.arange(27)]
			self.record[4]["Std"][0] += [0]
			self.record[4]["Std"][1] += [-1]
			self.record[4]["Std"][2] += [-1]
			self.record[4]["Std"][3] += [-1]
			self.record[4]["Std"][4] += [-1]
			self.record[4]["Std"][5] += [-1]
			self.record[4]["Std"][6] += [-1]
			self.record[4]["Std"][7] += [-1]
			self.record[4]["Std"][8] += [-1]
			self.record[4]["Std"][9] += [-1]
			self.record[4]["Std"][10] += [0]
			self.record[4]["Std"][11] += [0]
			self.record[4]["Std"][12] += [0]
			self.record[4]["Std"][13] += [0]
			self.record[4]["Std"][14] += [0]
			self.record[4]["Std"][15] += [0]
			self.record[4]["Std"][16] += [0]
			self.record[4]["Std"][17] += [0]
			if recOffCounts:
				self.record[4]["Std"][18] += [0]
				self.record[4]["Std"][19] += [0]
				self.record[4]["Std"][20] += [0]
				self.record[4]["Std"][21] += [0]
				self.record[4]["Std"][22] += [0]
				self.record[4]["Std"][23] += [0]
				self.record[4]["Std"][24] += [0]
				self.record[4]["Std"][25] += [0]
				self.record[4]["Std"][26] += [0]
			else:
				self.record[4]["Std"][18] += [None]
				self.record[4]["Std"][19] += [None]
				self.record[4]["Std"][20] += [None]
				self.record[4]["Std"][21] += [None]
				self.record[4]["Std"][22] += [None]
				self.record[4]["Std"][23] += [None]
				self.record[4]["Std"][24] += [None]
				self.record[4]["Std"][25] += [None]
				self.record[4]["Std"][26] += [None]
		for j in invCounters[0][0].keys():
			count = invCounters[0][0][j]
			if count != 0:
				if j not in self.record[4]:
					self.record[4][j] = [[] for x in np.arange(27)]
				self.record[4][j][0] += [count]
				self.record[4][j][1] += [np.mean(invDataDict[j][0])]
				self.record[4][j][2] += [np.var(invDataDict[j][0])]
				self.record[4][j][3] += [self.__randomizedMode(invDataDict[j][0])]
				self.record[4][j][4] += [np.mean(invDataDict[j][1])]
				self.record[4][j][5] += [np.var(invDataDict[j][1])]
				self.record[4][j][6] += [self.__randomizedMode(invDataDict[j][1])]
				self.record[4][j][7] += [np.mean(invDataDict[j][2])]
				self.record[4][j][8] += [np.var(invDataDict[j][2])]
				self.record[4][j][9] += [self.__randomizedMode(invDataDict[j][2])]
				self.record[4][j][10] += [invCounters[0][1][j]]
				self.record[4][j][11] += [invCounters[0][2][j]]
				self.record[4][j][12] += [invCounters[1][0][j]]
				self.record[4][j][13] += [invCounters[1][1][j]]
				self.record[4][j][14] += [invCounters[1][2][j]]
				self.record[4][j][15] += [invCounters[2][0][j]]
				self.record[4][j][16] += [invCounters[2][1][j]]
				self.record[4][j][17] += [invCounters[2][2][j]]
				if recOffCounts:
					self.record[4][j][18] += [invOffCounters[0][0][j]/2]
					self.record[4][j][19] += [invOffCounters[0][1][j]]
					self.record[4][j][20] += [invOffCounters[0][2][j]]
					self.record[4][j][21] += [invOffCounters[1][0][j]/2]
					self.record[4][j][22] += [invOffCounters[1][1][j]]
					self.record[4][j][23] += [invOffCounters[1][2][j]]
					self.record[4][j][24] += [invOffCounters[2][0][j]/2]
					self.record[4][j][25] += [invOffCounters[2][1][j]]
					self.record[4][j][26] += [invOffCounters[2][2][j]]
				else:
					self.record[4][j][18] += [None]
					self.record[4][j][19] += [None]
					self.record[4][j][20] += [None]
					self.record[4][j][21] += [None]
					self.record[4][j][22] += [None]
					self.record[4][j][23] += [None]
					self.record[4][j][24] += [None]
					self.record[4][j][25] += [None]
					self.record[4][j][26] += [None]
			else:
				raise self.SimulationStateError((j,count),"Recording had an inversion counter "+\
					"ID/key with value 0, should never happen with a counter - Inversion ID "+str(j))
		# Record Population Statistics
		# MAYBE JUST RECORD ALL THE VALUES FOR LATER ANALYSIS?
		self.record[6] += [np.mean(popSurvVals)]
		self.record[7] += [np.var(popSurvVals)]
		self.record[8] += [np.mean(popReprVals)]
		self.record[9] += [np.var(popReprVals)]

		# Record Haplotype Statistics
		for h in hapRange:
			hapClassRange = self.__calcHapClassRange(h)
			if h not in self.record[12]:
				self.record[12][h] = [[[] for i in hapClassRange] for x in np.arange(18)]
				assert(self.record[11][h][3] == currRecIndex), "First record of the haplotype does not match the estimated index"
			for i in hapClassRange:
				self.record[12][h][0][i] += [hapCounters[h][0][0][i]]
				self.record[12][h][1][i] += [hapCounters[h][0][1][i]]
				self.record[12][h][2][i] += [hapCounters[h][0][2][i]]
				self.record[12][h][3][i] += [hapCounters[h][1][0][i]]
				self.record[12][h][4][i] += [hapCounters[h][1][1][i]]
				self.record[12][h][5][i] += [hapCounters[h][1][2][i]]
				self.record[12][h][6][i] += [hapCounters[h][2][0][i]]
				self.record[12][h][7][i] += [hapCounters[h][2][1][i]]
				self.record[12][h][8][i] += [hapCounters[h][2][2][i]]
				if recOffCounts:
					self.record[12][h][9][i] += [hapOffCounters[h][0][0][i]/2]
					self.record[12][h][10][i] += [hapOffCounters[h][0][1][i]]
					self.record[12][h][11][i] += [hapOffCounters[h][0][2][i]]
					self.record[12][h][12][i] += [hapOffCounters[h][1][0][i]/2]
					self.record[12][h][13][i] += [hapOffCounters[h][1][1][i]]
					self.record[12][h][14][i] += [hapOffCounters[h][1][2][i]]
					self.record[12][h][15][i] += [hapOffCounters[h][2][0][i]/2]
					self.record[12][h][16][i] += [hapOffCounters[h][2][1][i]]
					self.record[12][h][17][i] += [hapOffCounters[h][2][2][i]]
				else:
					self.record[12][h][9][i] += [None]
					self.record[12][h][10][i] += [None]
					self.record[12][h][11][i] += [None]
					self.record[12][h][12][i] += [None]
					self.record[12][h][13][i] += [None]
					self.record[12][h][14][i] += [None]
					self.record[12][h][15][i] += [None]
					self.record[12][h][16][i] += [None]
					self.record[12][h][17][i] += [None]


	# # Checks the record to see if all mutations/inversions are fixed
	# def checkAllRecMutFixed(self):
	# 	lastRecord = len(self.record[0])-1
	# 	for m in range(self.__mutIDcount):
	# 		mutCount = self.record[2][m][lastRecord]
	# 		if not (mutCount == 0 or mutCount == 2*self.size):
	# 			return False
	# 	for i in range(self.__invIDcount):
	# 		invCount = self.record[4][i][lastRecord]
	# 		if not (invCount == 0 or invCount == 2*self.size):
	# 			return False
	# 	return True

	# 
	def guessBalancedHaplotypeDef(self,numMut=4,numInv=2,minCount=0):
		mutRange = np.arange(self.__mutIDcount)
		# polyMuts = []
		polyMutScores = []
		for mut in mutRange:
			if not self.__locFixed[mut] and mut in self.record[2] and len(self.record[2][mut][0])>0:
				mutCount = self.record[2][mut][0][-1]
				if mutCount > minCount:
					mutAge = self.age - self.record[1][mut][4]
					fixDistScale = (self.size-abs(self.size-mutCount))**2
					score = fixDistScale*mutAge
					polyMutScores += [(mut,score)]
					# polyMuts += [mut]

		invRange = np.arange(self.__invIDcount)
		# polyInvs = []
		polyInvScores = []
		for inv in invRange:
			if not self.__invFixed[inv] and inv in self.record[4] and len(self.record[4][inv][0])>0:
				invCount = self.record[4][inv][0][-1]
				if invCount > minCount:
					invAge = self.age - self.record[3][inv][3]
					fixDistScale = (self.size-abs(self.size-invCount))**2
					score = fixDistScale*invAge
					polyInvScores += [(inv,score)]
					# polyInvs += [inv]

		hapDef = [[],[]]
		if len(polyMutScores) <= numMut:
			hapDef[0] = [m[0] for m in polyMutScores]
		else:
			sortedMuts = sorted(polyMutScores,reverse=True,key=lambda x: x[1])
			chosenMuts = [sortedMuts[x][0] for x in np.arange(numMut)]
			hapDef[0] = chosenMuts
		if len(polyInvScores) <= numInv:
			hapDef[1] = [i[0] for i in polyInvScores]
		else:
			sortedInvs = sorted(polyInvScores,reverse=True,key=lambda x: x[1])
			chosenInvs = [sortedInvs[x][0] for x in np.arange(numInv)]
			hapDef[1] = chosenInvs

		return hapDef

	def samplePolyHap(self,numMut=2,numInv=1,minCount=0):
		mutRange = np.arange(self.__mutIDcount)
		polyMuts = []
		for mut in mutRange:
			if not self.__locFixed[mut] and mut in self.record[2] and len(self.record[2][mut][0])>0:
				mutCount = self.record[2][mut][0][-1]
				if mutCount > minCount:
					polyMuts += [mut]
		invRange = np.arange(self.__invIDcount)
		polyInvs = []
		for inv in invRange:
			if not self.__invFixed[inv] and inv in self.record[4] and len(self.record[4][inv][0])>0:
				invCount = self.record[4][inv][0][-1]
				if invCount > minCount:
					polyInvs += [inv]
		hapDef = [[],[]]
		if len(polyMuts) <= numMut:
			hapDef[0] = polyMuts
		else:
			hapDef[0] = list(np.random.choice(polyMuts,size=numMut,replace=False))
		if len(polyInvs) <= numInv:
			hapDef[1] = polyInvs
		else:
			hapDef[1] = list(np.random.choice(polyInvs,size=numInv,replace=False))
		return hapDef


	# For adding a haplotype to the record at some later generation
	# IF ADDED DURING A RECORDING GENERATION, MUST BE ADDED BEFORE ZYGOTE RECORDING
	def addRecordedHaplotype(self,hapDef):
		recHapDef = hapDef + [self.age,len(self.record[0])]
		self.recordHaps = True
		self.record[11] += [recHapDef]
		return

	def addRecordedHaplotypes(self,hapDefList):
		for hapDef in hapDefList:
			self.addRecordedHaplotype(hapDef)
		return


	# Checks the record to see if all mutations/inversions are fixed (or lost)
	def checkAllRecMutFixed(self):
		if self.verbose:
			print(self.__locFixed)
			print(self.__mutLost)
			print(self.__invFixed)
			print(self.__invLost)
			print(all(self.__locFixed) and all(self.__invFixed))
		return all(self.__locFixed) and all(self.__invFixed)


	# For simulating a number of generations sequentially, without recording
	def stepNGens(self, numGenerations):
		for i in range(numGenerations):
			self.step()
		return

	# For running a step and updating the record,
	#   including recording parental values
	def recordStep(self):
		if self.recordReproductiveAdults:
			self.step(recReprodParents=True)
		else:
			self.step()
		self.__updateRecord()
		return

	# For running a step and only recording arrangement data,
	#   including recording adult/parental values
	def arrRecordStep(self):
		if self.recordReproductiveAdults:
			self.step(recParentalArr=True)
		else:
			self.step()
		self.__recAllArrangementData()
		return

	# For running a step and then updating the record.
	# CONSIDER JUST MAKING __updateRecord PUBLIC
	# def recordStep(self):
	# 	self.step()
	# 	self.__updateRecord()
	# 	return

	# For simulating a number of generations sequentially, recording after each
	def recordNGens(self, numGenerations):
		for i in range(numGenerations):
			self.recordStep()
		return

	# For simulating a number of generations sequentially, recording at the end
	def recordNthGen(self, numGenerations, doubleRecord = False):
		if doubleRecord:
			self.stepNGens(numGenerations-2)
			self.recordStep()
			self.recordStep()
		else:
			self.stepNGens(numGenerations-1)
			self.recordStep()
		return

	# For simulating a number of generations sequentially, recording at the end
	def recordArrNthGen(self, numGenerations, doubleRecord = False):
		if doubleRecord:
			self.stepNGens(numGenerations-2)
			self.arrRecordStep()
			self.arrRecordStep()
		else:
			self.stepNGens(numGenerations-1)
			self.arrRecordStep()
		return

	# For simulating setSize*setNum generations sequentially, recording every setSize'th generation
	def recordNthForNSets(self, setSize, setNum, doubleRecord = False):
		if doubleRecord:
			for i in range(setNum):
				self.stepNGens(setSize-2)
				self.recordStep()
				self.recordStep()
		else:
			for i in range(setNum):
				self.stepNGens(setSize-1)
				self.recordStep()
		return


	# # For simulating setSize*setNum generations sequentially, recording every setSize'th generation
	# def recordNthForNSets(self, setSize, setNum):
	# 	for i in range(setNum):
	# 		self.stepNGens(setSize)
	# 		self.__updateRecord()
	# 	return

	# For simulating numGens generations sequentially, recording every setSize'th generation
	def recordNthInNGens(self, setSize, numGens, doubleRecord = False):
		# Be wary of this, it may not simulate the full number of generations,
		#   but really you wouldn't necessarily record the final step anyway
		numSets = int(numGens/setSize)
		if doubleRecord:
			for i in range(numSets):
				self.stepNGens(setSize-2)
				self.recordStep()
				self.recordStep()
		else:
			for i in range(numSets):
				self.stepNGens(setSize-1)
				self.recordStep()
		return


	# Helper function to factor out code used in filling out the record after fixation of the variants
	def __padRecord(self,gen,lhs,fixedHapIDs):
		alleleCount = self.size*2
		femalePop = self.size//2
		malePop = self.size - femalePop
		self.record[10] += [lhs] # Life history stage
		if lhs == "Adult":
			# ASSUMES that this is the parental stage being recorded
			self.record[0] += [gen] # Parental generation
		elif lhs == "Zygote":
			self.record[0] += [gen] # Zygote generation
		else:
			raise self.SimulationStateError((gen,lhs), "Error: unrecognized life history stage "+str(lhs))
		self.record[5] += [self.size] # Population size
		self.record[13] += [femalePop] # Number of females
		self.record[14] += [malePop] # Number of males
		self.record[15] += [None] # Female reproductive distribution
		self.record[16] += [None] # Male reproductive distribution
		self.record[6] += [1] # Mean individual survival proportion in pop
		self.record[7] += [0] # Variance in survival proportion
		self.record[8] += [0] # Mean individual display value in pop
		self.record[9] += [0] # Variance in display value
		# Record Haplotype Statistics
		for h in np.arange(len(self.record[11])):
			hapClassRange = self.__calcHapClassRange(h)
			for i in hapClassRange:
				self.record[12][h][2][i] += [0]
				self.record[12][h][5][i] += [0]
				self.record[12][h][8][i] += [0]
				self.record[12][h][11][i] += [0]
				self.record[12][h][14][i] += [0]
				self.record[12][h][17][i] += [0]
				if i == fixedHapIDs[h]:
					self.record[12][h][0][i] += [alleleCount]
					self.record[12][h][1][i] += [self.size]
					self.record[12][h][3][i] += [femalePop*2]
					self.record[12][h][4][i] += [femalePop]
					self.record[12][h][6][i] += [malePop*2]
					self.record[12][h][7][i] += [malePop]
				else:
					self.record[12][h][0][i] += [0]
					self.record[12][h][1][i] += [0]
					self.record[12][h][3][i] += [0]
					self.record[12][h][4][i] += [0]
					self.record[12][h][6][i] += [0]
					self.record[12][h][7][i] += [0]
				if lhs == "Adult":
					if i == fixedHapIDs[h]:
						self.record[12][h][9][i] += [self.size]
						self.record[12][h][10][i] += [self.size]
						self.record[12][h][12][i] += [femalePop]
						self.record[12][h][13][i] += [femalePop]
						self.record[12][h][15][i] += [malePop]
						self.record[12][h][16][i] += [malePop]
					else:
						self.record[12][h][9][i] += [0]
						self.record[12][h][10][i] += [0]
						self.record[12][h][12][i] += [0]
						self.record[12][h][13][i] += [0]
						self.record[12][h][15][i] += [0]
						self.record[12][h][16][i] += [0]

		return

	# For simulating setSize*setNum generations sequentially, recording every setSize'th generation
	# Recording every generation after fixation as the fixation frequencies if fillRecord = True
	# Ideally only used with non-mutating populations
	def recordNthInNStopWhenFixed(self, setSize, numGens, finishRecord = True): #recReprodParents=False
		# if willMutate or willMutInv:
		# 	raise SAIpop.InputError((willMutate,willMutInv),'For use with non-mutating populations')
		# Be wary of recording every N, it may not simulate the full number of generations,
		#   but really you wouldn't necessarily record the final step anyway if it falls within a set
		print("Recording every "+str(setSize)+" generations, stopping when all variants fix")
		numSets = int(numGens/setSize) # //
		for s in range(numSets):
			self.stepNGens(setSize-1)
			self.recordStep()
			if self.checkAllRecMutFixed():
				print("All variants fixed or lost")
				if finishRecord:
					lastRecord = len(self.record[0])-1
					if self.verbose:
						print(lastRecord)
					g = self.record[0][lastRecord]
					# Figure out the fixed haplotype class for each hap definition:
					fixedHapIDs = []
					hapRange = np.arange(len(self.record[11]))
					for h in hapRange:
						hapDef = self.record[11][h]
						fixedHap = [[],[]]
						for m in hapDef[0]:
							if self.__locFixed[m] and not self.__mutLost[m]:
								fixedHap[0] += [m]
						for i in hapDef[1]:
							if self.__invFixed[i] and not self.__invLost[i]:
								fixedHap[1] += [i]
						fixedHapIDs += [self.__indexFromHap(fixedHap,hapDef)]
						# Double check that the records have been generated?
						if h not in self.record[12]:
							print("WARNING - all variants fixed so terminating simulation and padding record, "\
								"but missing record for defined haplotype "+str(h)+": "+str(hapDef))
							hapClassRange = self.__calcHapClassRange(h)
							self.record[12][h] = [[[] for i in hapClassRange] for x in np.arange(27)]
					# Finish the record
					for addRecord in range(1,numSets-s):
						if self.recordReproductiveAdults:
							self.__padRecord(g+addRecord*setSize-1,"Adult",fixedHapIDs)
						self.__padRecord(g+addRecord*setSize,"Zygote",fixedHapIDs)
				return
		return


	# Factoring out the code to add a newly observed class of arrangement
	# Add a dictionary element keyed to the arrangement for storing accounted data:
	#      [0] allele count of the arrangement
	#      [1] count of the number of individuals homozygous for the arrangement (implicitly the derived arrangement)
	#      [2] count of the number of individuals heterozygous for the arrangement
	#      [3] count of the number of arrangements in females
	#      [4] count of the number of females homozygous for the arrangement (implicitly the derived arrangement)
	#      [5] count of the number of females heterozygous for the arrangement
	#      [6] count of the number of arrangements in males
	#      [7] count of the number of males homozygous for the arrangement (implicitly the derived arrangement)
	#      [8] count of the number of males heterozygous for the arrangement
	#      [9] list of the count of mutations on the chromosomal arrangement for each chromosomal copy
	#      [10] list of the total survival effect of the chromosomal arrangement for each chromosomal copy
	#      [11] list of the total male competitive effect of the chromosomal arrangement for each chromosomal copy

	def __checkAddArrDat(self,arrClass,arrDataDict):
		if arrClass not in arrDataDict:
			arrDataDict[arrClass] = [0,0,0,0,0,0,0,0,0,[],[],[]]
		return

	def __updateArrCount(self,arrClass,zygosity,sex,arrDataDict):
		if zygosity == "Hom":
			arrDataDict[arrClass][0] += 2
			arrDataDict[arrClass][1] += 1
			if sex == 'F':
				arrDataDict[arrClass][3] += 2
				arrDataDict[arrClass][4] += 1
			else:
				arrDataDict[arrClass][6] += 2
				arrDataDict[arrClass][7] += 1
		elif zygosity == "Het":
			arrDataDict[arrClass][0] += 1
			arrDataDict[arrClass][2] += 1
			if sex == 'F':
				arrDataDict[arrClass][3] += 1
				arrDataDict[arrClass][5] += 1
			else:
				arrDataDict[arrClass][6] += 1
				arrDataDict[arrClass][8] += 1
		else:
			raise self.SimulationStateError(zygosity,"Given zygosity value unexpected, "+\
				"expecting 'Hom' or 'Het'.")
		return

	def __updateArrValues(self,arrClass,chromMuts,sex,arrDataDict):
		numMuts = 0
		chromSur = 1
		chromRep = 0
		for mut in chromMuts:
			numMuts += 1
			chromSur *= mut[1]
			chromRep += mut[2]
		arrDataDict[arrClass][9] += [numMuts]
		arrDataDict[arrClass][10] += [chromSur]
		arrDataDict[arrClass][11] += [chromRep]
		return

	# For collecting record info on a pair of chromosomes from an individual
	def __countArrangementGenotypes(self,chrom,sex,arrDataDict):
		# Make inversion ID tuples as arrangement dictionary keys
		arrHom1 = ()
		invHom1 = chrom[0][1]
		if len(invHom1) > 0:
			invIDs = [inv[2] for inv in invHom1]
			arrHom1 = tuple(invIDs)
		arrHom2 = ()
		invHom2 = chrom[1][1]
		if len(invHom2) > 0:
			invIDs = [inv[2] for inv in invHom2]
			arrHom2 = tuple(invIDs)

		self.__checkAddArrDat(arrHom1,arrDataDict)
		if arrHom1 == arrHom2:
			self.__updateArrCount(arrHom1,'Hom',sex,arrDataDict)
		else:
			self.__checkAddArrDat(arrHom2,arrDataDict)
			self.__updateArrCount(arrHom1,'Het',sex,arrDataDict)
			self.__updateArrCount(arrHom2,'Het',sex,arrDataDict)

		mutHom1 = chrom[0][0]
		self.__updateArrValues(arrHom1,mutHom1,sex,arrDataDict)
		mutHom2 = chrom[1][0]
		self.__updateArrValues(arrHom2,mutHom2,sex,arrDataDict)

		return

	# For profiling all arrangements in a population
	# Returns a dictionary of collected arrangement data, keyed by the tuple of associated inversion IDs:
	#      [0] counts of the arrangement
	#      [1] count of the number of individuals homozygous for the arrangement (implicitly the derived arrangement)
	#      [2] count of the number of individuals heterozygous for the arrangement
	#      [3] count of the number of arrangements in females
	#      [4] count of the number of females homozygous for the arrangement (implicitly the derived arrangement)
	#      [5] count of the number of females heterozygous for the arrangement
	#      [6] count of the number of arrangements in males
	#      [7] count of the number of males homozygous for the arrangement (implicitly the derived arrangement)
	#      [8] count of the number of males heterozygous for the arrangement
	#      [9] average count of mutations within the buffer region of the arrangement 
	#        across the entire population
	#      [10] variance in count of mutations within the buffer region of the arrangement 
	#        across the entire population
	#      [11] mode of the count of mutations within the buffer region of the arrangement 
	#        across the entire population
	#      [12] average total survival effect within the buffer region of the arrangement 
	#        across the entire population
	#      [13] variance in total survival effect within the buffer region of the arrangement 
	#        across the entire population
	#      [14] mode total survival effect within the buffer region of the arrangement 
	#        across the entire population
	#      [15] average total reproductive effect within the buffer region of the arrangement 
	#        across the entire population
	#      [16] variance in total reproductive effect within the buffer region of the arrangement 
	#        across the entire population
	#      [17] mode total reproductive effect within the buffer region of the arrangement 
	#        across the entire population
	def genAllArrangementData(self,lifeHistStage="Zygote",inFemales=None,inMales=None):
		print ("All Arrangement Check")

		assert(lifeHistStage == "Zygote" or (inFemales is not None and inMales is not None)), "Life history stages "+\
			"other than zygote must have an input of female and male individuals"
		assert(lifeHistStage != "Zygote" or (self.size == len(self.males)+len(self.females))), "At zygote life history stage, "+\
			"expect self.males and self.females to total the set population size"

		arrRec = {}

		# Assemble the relevant individuals/stats for the LHS
		if lifeHistStage != "Zygote":
			currIndividuals = inFemales+inMales
		else:
			currIndividuals = self.females+self.males

		# Count arrangements and their genetic values
		arrDataDict = {}
		for i in np.arange(len(currIndividuals)):
			indiv = currIndividuals[i]
			indivSur = 1
			indivRep = 0
			for chrom in indiv.genome:
				self.__countArrangementGenotypes(chrom,indiv.sex,arrDataDict)

		# Rearrange data, calculate the mean, variance, and mode of values
		for arrClass in arrDataDict:
			arrRec[arrClass] = arrDataDict[arrClass][0:9]+[
								np.mean(arrDataDict[arrClass][9]),
								np.var(arrDataDict[arrClass][9]),
								self.__randomizedMode(arrDataDict[arrClass][9]),
								np.mean(arrDataDict[arrClass][10]),
								np.var(arrDataDict[arrClass][10]),
								self.__randomizedMode(arrDataDict[arrClass][10]),
								np.mean(arrDataDict[arrClass][11]),
								np.var(arrDataDict[arrClass][11]),
								self.__randomizedMode(arrDataDict[arrClass][11])
								]

		return arrRec

	def __recAllArrangementData(self,lifeHistStage="Zygote",inFemales=None,inMales=None):
		arrRec = self.genAllArrangementData(lifeHistStage,inFemales,inMales)
		self.arrangementRecord += [[self.age,lifeHistStage,arrRec]]
		return

	# Writing the arrangement details for a specific arrangement
	def __writeArrStep(self,outfile,recordGen,recordStage,arrRec):
		for arrClass in arrRec:
			outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(arrClass) + '\t'\
				+ str(arrRec[arrClass][0]) + '\t' + str(arrRec[arrClass][1]) + '\t'\
				+ str(arrRec[arrClass][2]) + '\t' + str(arrRec[arrClass][3]) + '\t'\
				+ str(arrRec[arrClass][4]) + '\t' + str(arrRec[arrClass][5]) + '\t'\
				+ str(arrRec[arrClass][6]) + '\t' + str(arrRec[arrClass][7]) + '\t'\
				+ str(arrRec[arrClass][8]) + '\t' + str(arrRec[arrClass][9]) + '\t'\
				+ str(arrRec[arrClass][10]) + '\t' + str(arrRec[arrClass][11]) + '\t'\
				+ str(arrRec[arrClass][12]) + '\t' + str(arrRec[arrClass][13]) + '\t'\
				+ str(arrRec[arrClass][14]) + '\t' + str(arrRec[arrClass][15]) + '\t'\
				+ str(arrRec[arrClass][16]) + '\t' + str(arrRec[arrClass][17]) + '\n')
		return

	# Takes a filename and arrangement record dictionary
	# Writes a tab delineated file of the arrangment inversion ID tuple,
	#  allele count, average/var/mode number of mutations, average/var/mode cumulative survival,
	#  and average/var/mode cumulative male competitive effect
	def writeArrangements(self,filename):
		outfile = open(filename, 'w')
		outfile.write('Generation\tLifeHistoryStage\tArrangement\tCount\tHomozygote\tHeterozygote\t'+\
			'FemaleCount\tFemaleHomozygote\tFemaleHeterozygote\tMaleCount\tMaleHomozygote\tMaleHeterozygote\t'+\
			'AvgNumMut\tVarNumMut\tModeNumMut\tAvgSurEff\tVarSurEff\tModeSurEff\tAvgRepEff\tVarRepEff\tModeRepEff\n')
		for recordGen,recordStage,arrRec in self.arrangementRecord:
			self.__writeArrStep(outfile,recordGen,recordStage,arrRec)
		outfile.close()
		return True



	def __chromHashableIter(self,individuals,chrom):
		for i in individuals:
			for hom in [0,1]:
				yield HashableList(i.genome[chrom][hom])

	def __countHaps(self,individuals,chrom):
		return Counter(self.__chromHashableIter(individuals,chrom))

	# For generating and printing chromosome counts and frequencies in a readable manner
	def __getStrChromCounts(self,individuals):
		string = ''
		for chrom in np.arange(self.numChrom):
			string += "Haplotypes of Chromosome "+str(chrom)+":\n"
			c = self.__countHaps(individuals,chrom)
			haps = np.array(list(c.keys()))
			counts = np.array(list(c.values()))
			total = np.sum(counts)
			freqs = np.divide(counts,total)
			for i in np.arange(freqs.size):
				line = "{:7.4f}".format(freqs[i])+'\t{:6d}'.format(counts[i])+\
					'\t'+str(haps[i])+'\n'
				string += line
		return string

	def __getStrPopChromCounts(self):
		string = ''
		string += str(len(self.females))+' Females:\n'
		string += self.__getStrChromCounts(self.females)
		string += str(len(self.males))+' Males:\n'
		string += self.__getStrChromCounts(self.males)
		return string

	# For generating and printing the chromosome counts and frequencies in a readable manner
	def __printPopChromCounts(self):
		print(self.__getStrPopChromCounts())
		return


	# Functions dealing with writing output to file or the terminal, mainly from self.record

	def printRecord(self):
		print('\n\nRecorded Population Data at Generation '+str(self.age)+\
			'\n-----------------------------------------------------')
		print('\nRecorded Generations:\n'+str(self.record[0]))
		print('\nRecorded Life History Stages:\n'+str(self.record[10]))
		print('\nRecorded Mutations:\n'+str(self.record[1]))
		print('\nIndividual Mutation Count Records:\n'+str(self.record[2]))
		print('\nRecorded Inversions:\n'+str(self.record[3]))
		print('\nIndividual Inversion Records:\n'+str(self.record[4]))
		print('\nPopulation Sizes:\n'+str(self.record[5]))
		print('\nPopulation Mean Survival Proportion:\n'+str(self.record[6]))
		print('\nPopulation Variance in Survival Proportion:\n'+str(self.record[7]))
		print('\nPopulation Mean Display Value:\n'+str(self.record[8]))
		print('\nPopulation Variance in Display Value:\n'+str(self.record[9]))
		print('\nRecorded Haplotypes:\n'+str(self.record[11]))
		print('\nIndividual Haplotype Class Count Records:\n'+str(self.record[12]))
		print('\nNumber of females:\n'+str(self.record[13]))
		print('\nNumber of males:\n'+str(self.record[14]))
		print('\nNumber of offspring per mother at each recording of Adults:\n'+str(self.record[15]))
		print('\nNumber of offspring per father at each recording of Adults:\n'+str(self.record[16]))
		return


	def printGenomes(self):
		print('\nIndividual Genomes at Generation '+str(self.age)+\
			'\n-----------------------------------------------------')
		print(str(len(self.females))+' Females:')
		i = 0
		for indiv in self.females:
			indiv.printGenome(i)
			# indStr = 'Ind {:6d} '.format(i)
			# # print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		homStr = indStr + 'Chrom '+str(c)+' Hom '+str(h)
			# 		mutStr = homStr + ' Mut: ' + str(indiv.genome[c][h][0])
			# 		invStr = homStr + ' Inv: ' + str(indiv.genome[c][h][1])
			# 		print(mutStr)
			# 		print(invStr)
			# 		# print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		# print('Mutations: '+str(indiv.genome[c][h][0]))
			# 		# print('Inversions: '+str(indiv.genome[c][h][1]))
			i += 1
		print('\n'+str(len(self.males))+' Males:')
		for indiv in self.males:
			indiv.printGenome(i)
			# indStr = 'Ind {:6d} '.format(i)
			# # print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		homStr = indStr + 'Chrom '+str(c)+' Hom '+str(h)
			# 		mutStr = homStr + ' Mut: ' + str(indiv.genome[c][h][0])
			# 		invStr = homStr + ' Inv: ' + str(indiv.genome[c][h][1])
			# 		print(mutStr)
			# 		print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		print('Mutations: '+str(indiv.genome[c][h][0]))
			# 		print('Inversions: '+str(indiv.genome[c][h][1]))
			i += 1
		return

	def printInversionHaps(self):
		print('\nIndividual Inversion Haplotypes at Generation '+str(self.age)+\
			'\n-----------------------------------------------------')
		print(str(len(self.females))+' Females:')
		i = 0
		for indiv in self.females:
			indiv.printGenome(i,printMut=False)
			# indStr = 'Ind {:6d} '.format(i)
			# # print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		print(str(indiv.genome[c][h][1]))
			i += 1
		print('\n'+str(len(self.males))+' Males:')
		for indiv in self.males:
			indiv.printGenome(i,printMut=False)
			# print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		print(str(indiv.genome[c][h][1]))
			i += 1
		return

	def printMutationHaps(self):
		print('\nIndividual Mutation Haplotypes at Generation '+str(self.age)+\
			'\n-----------------------------------------------------')
		print(str(len(self.females))+' Females:')
		i = 0
		for indiv in self.females:
			indiv.printGenome(i,printInv=False)
			# print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		print(str(indiv.genome[c][h][0]))
			i += 1
		print('\n'+str(len(self.males))+' Males:')
		for indiv in self.males:
			indiv.printGenome(i,printInv=False)
			# print('Individual '+str(i)+':')
			# for c in range(self.numChrom):
			# 	for h in range(2):
			# 		print('Chromosome '+str(c)+' Homolog '+str(h)+':')
			# 		print(str(indiv.genome[c][h][0]))
			i += 1
		return




	# General helper functions for writing Population objects to a pickle file
	def __storePickleBytes(self,popBytes,out_file_name):
		with open(out_file_name,"wb") as outfile:
			outfile.write(popBytes)
		return

	def writePopPickle(self,out_file_name):
		import pickle as p
		with open(out_file_name,"wb") as outfile:
			popBytes = p.dumps(self)
			outfile.write(popBytes)
		return



	# # Writing the mutation details from a single generation, for writeMutation()
	# # Formatted as Generateion	
	# def __writeMutStep(self,outfile,recordGen,ID,i):
	# 	outfile.write(str(recordGen) + '\t' + str(self.record[4][ID][i]) + '\t'\
	# 		+ str(self.record[5][ID][i]) + '\t' + str(self.record[6][ID][i]) + '\t'\
	# 		+ str(self.record[7][ID][i]) + '\t' + str(self.record[8][ID][i]) + '\t'\
	# 		+ str(self.record[9][ID][i]) + '\t' + str(self.record[10][ID][i]) + '\n')

	# Takes a filename and mutation ID
	# Writes a tab delineated file of the generation and count data for that mutation, returns True
	# If no record of the mutation was maintained, writes  no file, returns False
	# ASSUMES 'Adult' lhs means offspring counts are recorded
	def writeMutation(self,filename,ID,fillRecord = False):
		# Exit if there is no record
		if ID not in self.record[2] or len(self.record[2][ID][0]) == 0:
			return False
		outfile = open(filename, 'w')
		outfile.write('Generation\tLifeHistoryStage\tSex\tCountType\tCount\tHomozygote\tHeterozygote\n')
		# outfile.write(str(self.record[0][0]) + '\t' + str(self.record[2][ID][0]) + '\n')
		# outfile.write(str(self.record[1][ID][4]) + '\t' + str(1) + '\n')
		# for t in range(1,len(self.record[0])):
		countRange = np.arange(len(["Allele","Homozygote","Heterozygote"]))
		recSexes = ["All","Female","Male"]
		sexRange = np.arange(len(recSexes))
		(firstGen, finalGen) = self.record[1][ID][4:6]
		t = 0
		recordGen = self.record[0][t]
		recordStage = self.record[10][t]
		while recordGen < firstGen:
			if fillRecord:
				for s in sexRange:
					outfile.write(str(recordGen)+'\t'+str(recordStage)+'\t'+recSexes[s]+'\t'+'Allele'+'\t0\t0\t0\n')
				if recordStage == "Adult":
					outfile.write(str(recordGen)+'\t'+str(recordStage)+'\t'+recSexes[s]+'\t'+'Offspring'+'\t0\t0\t0\n')
			t += 1
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
		offset = t
		m = t - offset
		recLen = len(self.record[2][ID][0])
		while m < recLen:
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
			for s in sexRange:
				outfile.write(str(recordGen) + '\t' + str(recordStage) +'\t'+recSexes[s]+'\t'+'Allele' + \
					'\t' + str(self.record[2][ID][s*3+0][m]) + '\t' + str(self.record[2][ID][s*3+1][m]) + \
					'\t' + str(self.record[2][ID][s*3+2][m]) + '\n')
			if recordStage == "Adult":
				for s in sexRange:
					outfile.write(str(recordGen) + '\t' + str(recordStage) +'\t'+recSexes[s]+'\t'+'Offspring' + \
						'\t' + str(self.record[2][ID][s*3+9][m]) + '\t' + str(self.record[2][ID][s*3+10][m]) + \
						'\t' + str(self.record[2][ID][s*3+11][m]) + '\n')
			t += 1
			m = t - offset
		# print(t)
		# print(self.record[0][-1])
		if fillRecord:
			for t in range(t,len(self.record[0])):
				recordGen = self.record[0][t]
				recordStage = self.record[10][t]
				if self.__locFixed[ID] and not self.__mutLost[ID]:
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[0]+'\t'+'Allele' + \
						'\t' + str(2*self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0\n')
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[1]+'\t'+'Allele' + \
						'\t' + str(2*self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0\n')
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[2]+'\t'+'Allele' + \
						'\t' + str(2*self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0\n')
					if recordStage == "Adult":
						outfile.write(str(recordGen) + '\t' + str(recordStage) +'\t'+recSexes[0]+'\t'+'Offspring' + \
							'\t' + str(self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0\n')
						outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[1]+'\t'+'Offspring' + \
							'\t' + str(self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0\n')
						outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[2]+'\t'+'Offspring' + \
							'\t' + str(self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0\n')
				else:
					for s in sexRange:
						outfile.write(str(recordGen)+'\t'+str(recordStage)+'\t'+recSexes[s]+'\t'+'Allele'+'\t0\t0\t0\n')
					if recordStage == "Adult":
						for s in sexRange:
							outfile.write(str(recordGen)+'\t'+str(recordStage)+'\t'+recSexes[s]+'\t'+'Offspring'+'\t0\t0\t0\n')
		# offset = t
		# if finalGen < firstGen:
		# 	for t in range(t,len(self.record[0])):
		# 		m = t - offset
		# 		recordGen = self.record[0][t]
		# 		recordStage = self.record[10][t]
		# 		for s in sexRange:
		# 			outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[s] + \
		# 				'\t' + str(self.record[2][ID][s*3+0][m]) + '\t' + str(self.record[2][ID][s*3+1][m]) + \
		# 				'\t' + str(self.record[2][ID][s*3+2][m]) + '\n')
		# else:
		# 	while recordGen < finalGen:
		# 		m = t - offset
		# 		for s in sexRange:
		# 			outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[s] + \
		# 				'\t' + str(self.record[2][ID][s*3+0][m]) + '\t' + str(self.record[2][ID][s*3+1][m]) + \
		# 				'\t' + str(self.record[2][ID][s*3+2][m]) + '\n')
		# 		t += 1
		# 		recordGen = self.record[0][t]
		# 		recordStage = self.record[10][t]
		# 	if fillRecord:
		# 		for t in range(t,len(self.record[0])):
		# 			recordGen = self.record[0][t]
		# 			recordStage = self.record[10][t]
		# 			if self.__locFixed[ID] and not self.__mutLost[ID]:
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[0] + \
		# 					'\t' + str(2*self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0\n')
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[1] + \
		# 					'\t' + str(2*self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0\n')
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[2] + \
		# 					'\t' + str(2*self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0\n')
		# 			else:
		# 				for s in sexRange:
		# 					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + recSexes[s] + '\t0\t0\t0\n')
		outfile.close()
		return True

	# Writing the inversion details from a single generation, for writeInversion()
	def __writeInvStep(self,outfile,recordGen,recordStage,ID,i):
		outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(self.record[4][ID][0][i]) + '\t'\
			+ str(self.record[4][ID][10][i]) + '\t' + str(self.record[4][ID][11][i]) + '\t'\
			+ str(self.record[4][ID][12][i]) + '\t' + str(self.record[4][ID][13][i]) + '\t'\
			+ str(self.record[4][ID][14][i]) + '\t' + str(self.record[4][ID][15][i]) + '\t'\
			+ str(self.record[4][ID][16][i]) + '\t' + str(self.record[4][ID][17][i]) + '\t'\
			+ str(self.record[4][ID][1][i]) + '\t' + str(self.record[4][ID][2][i]) + '\t'\
			+ str(self.record[4][ID][3][i]) + '\t' + str(self.record[4][ID][4][i]) + '\t'\
			+ str(self.record[4][ID][5][i]) + '\t' + str(self.record[4][ID][6][i]) + '\t'\
			+ str(self.record[4][ID][7][i]) + '\t' + str(self.record[4][ID][8][i]) + '\t'\
			+ str(self.record[4][ID][9][i]) + '\t' + str(self.record[4][ID][18][i]) + '\t'\
			+ str(self.record[4][ID][19][i]) + '\t' + str(self.record[4][ID][20][i]) + '\t'\
			+ str(self.record[4][ID][21][i]) + '\t' + str(self.record[4][ID][22][i]) + '\t'\
			+ str(self.record[4][ID][23][i]) + '\t' + str(self.record[4][ID][24][i]) + '\t'\
			+ str(self.record[4][ID][25][i]) + '\t' + str(self.record[4][ID][26][i]) + '\n')

	# Takes a filename and inversion ID
	# Writes a tab delineated file of the generation, count, 
	#  average number of mutations in buffer, average cumulative survival,
	#  and average cumulative effect data for that inversion
	def writeInversion(self,filename,ID,fillRecord = False):
		# Exit if there is no record
		if ID not in self.record[4] or len(self.record[4][ID][0]) == 0:
			return False
		outfile = open(filename, 'w')
		# outfile.write('Inversion '+str(ID)+'\n')
		outfile.write('Generation\tLifeHistoryStage\tCount\tHomozygote\tHeterozygote\t'+\
			'FemaleCount\tFemaleHomozygote\tFemaleHeterozygote\tMaleCount\tMaleHomozygote\tMaleHeterozygote\t'+\
			'AvgNumMut\tVarNumMut\tModeNumMut\tAvgSurEff\tVarSurEff\tModeSurEff\tAvgRepEff\tVarRepEff\tModeRepEff\t'+\
			'OffspringCount\tOffspringCountInHomozygote\tOffspringCountInHeterozygote\t'+\
			'OffspringCountInFemale\tOffspringCountInFemaleHomozygote\tOffspringCountInFemaleHeterozygote\t'+\
			'OffspringCountInMale\tOffspringCountInMaleHomozygote\tOffspringCountInMaleHeterozygote\n')
		# outfile.write(str(self.record[0][0]) + '\t' + str(self.record[4][ID][0]) + '\t'\
		# 	+ str(self.record[5][ID][0]) + '\t' + str(self.record[6][ID][0]) + '\t'\
		# 	+ str(self.record[7][ID][0]) + '\n')
		# # May want to find a way to acount for the number of mutations in the inversion
		# #   and their effects upon generation
		# outfile.write(str(self.record[3][ID][3]) + '\t' + str(1) + '\t' + 'NA' + '\t'\
		# 	+ 'NA' + '\t' + 'NA' + '\n')
		# for t in range(1,len(self.record[0])):
		(firstGen, finalGen) = self.record[3][ID][3:5]
		t = 0
		recordGen = self.record[0][t]
		recordStage = self.record[10][t]
		while recordGen < firstGen:
			if fillRecord:
				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t0'*9+'\t-1'*9+'\t0'*9+'\n')
			t += 1
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
		offset = t
		i = t - offset
		recLen = len(self.record[4][ID][0])
		while i < recLen:
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
			self.__writeInvStep(outfile,recordGen,recordStage,ID,i)
			t += 1
			i = t - offset
		if fillRecord:
			for t in range(t,len(self.record[0])):
				recordGen = self.record[0][t]
				recordStage = self.record[10][t]
				if self.__invFixed[ID] and not self.__invLost[ID]:
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(2*self.record[5][t]) + \
						'\t' + str(self.record[5][t]) + '\t0' + '\t' + str(2*self.record[13][t]) + \
						'\t' + str(self.record[13][t]) + '\t0' + '\t' + str(2*self.record[14][t]) + \
						'\t' + str(self.record[14][t]) + '\t0' + '\t-1'*9+
						'\t' + str(self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0' + \
						'\t' + str(self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0' + \
						'\t' + str(self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0' + '\n')
				else:
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t0'*9+'\t-1'*9+'\t0'*9+'\n')
		# offset = t
		# if finalGen < firstGen:
		# 	for t in range(t,len(self.record[0])):
		# 		i = t - offset
		# 		recordGen = self.record[0][t]
		# 		recordStage = self.record[10][t]
		# 		self.__writeInvStep(outfile,recordGen,recordStage,ID,i)
		# 	# while t < len(self.record[0]):
		# 	# 	i = t - offset
		# 	# 	self.__writeInvStep(outfile,recordGen,ID,i)
		# 	# 	t += 1
		# 	# 	recordGen = self.record[0][t]
		# 	# i = t - offset
		# 	# self.__writeInvStep(outfile,recordGen,ID,i)
		# else:
		# 	while recordGen < finalGen:
		# 		i = t - offset
		# 		self.__writeInvStep(outfile,recordGen,recordStage,ID,i)
		# 		t += 1
		# 		recordGen = self.record[0][t]
		# 		recordStage = self.record[10][t]
		# 	if fillRecord:
		# 		for t in range(t,len(self.record[0])):
		# 			recordGen = self.record[0][t]
		# 			recordStage = self.record[10][t]
		# 			if self.__invFixed[ID] and not self.__invLost[ID]:
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(2*self.record[5][t]) + \
		# 					'\t' + str(self.record[5][t]) + '\t0' + '\t' + str(2*self.record[13][t]) + \
		# 					'\t' + str(self.record[13][t]) + '\t0' + '\t' + str(2*self.record[14][t]) + \
		# 					'\t' + str(self.record[14][t]) + '\t0' + '\t-1'*9+'\n')
		# 			else:
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t0'*9 + '\t-1'*9+'\n')
		# 		# 	if self.__invFixed[ID]:
		# 		# 		outfile.write(str(recordGen) + '\t' + str(self.record[11][t]) + '\t-1\t-1\t-1\n')
		# 		# 	else:
		# 		# 		outfile.write(str(recordGen) + '\t0\t-1\t-1\t-1\n')
		# 		# 	t += 1
		# 		# 	recordGen = self.record[0][t]
		# 		# if self.__invFixed[ID]:
		# 		# 	outfile.write(str(recordGen) + '\t' + str(self.record[11][t]) + '\t-1\t-1\t-1\n')
		# 		# else:
		# 		# 	outfile.write(str(recordGen) + '\t0\t-1\t-1\t-1\n')
		# # for t in range(len(self.record[0])):
		# # 	outfile.write(str(self.record[0][t]) + '\t' + str(self.record[4][ID][t]) + '\t'\
		# # 		+ str(self.record[5][ID][t]) + '\t' + str(self.record[6][ID][t]) + '\t'\
		# # 		+ str(self.record[7][ID][t]) + '\n')
		outfile.close()
		return True

	# Takes a filename,
	# Writes a tab delineated file of the generation, count, 
	#  average number of mutations, average cumulative survival,
	#  and average cumulative effect data
	# for standard arrangement chromosomes in the population
	def writeStdArrangement(self,filename):
		ID = "Std"
		# Exit if there is no record
		if ID not in self.record[4] or len(self.record[4][ID][0]) == 0:
			return False
		outfile = open(filename, 'w')
		outfile.write('Generation\tLifeHistoryStage\tCount\tHomozygote\tHeterozygote\t'+\
			'FemaleCount\tFemaleHomozygote\tFemaleHeterozygote\tMaleCount\tMaleHomozygote\tMaleHeterozygote\t'+\
			'AvgNumMut\tVarNumMut\tModeNumMut\tAvgSurEff\tVarSurEff\tModeSurEff\tAvgRepEff\tVarRepEff\tModeRepEff\t'+\
			'OffspringCount\tOffspringCountInHomozygote\tOffspringCountInHeterozygote\t'+\
			'OffspringCountInFemale\tOffspringCountInFemaleHomozygote\tOffspringCountInFemaleHeterozygote\t'+\
			'OffspringCountInMale\tOffspringCountInMaleHomozygote\tOffspringCountInMaleHeterozygote\n')
		for t in range(len(self.record[0])):
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
			self.__writeInvStep(outfile,recordGen,recordStage,ID,t)
		outfile.close()
		return True


	# For generating all possible haplotypes from the mutation and inversion lists,
	#   NOT GUARANTEED to be in the index order used in indexFromHap
	def __buildHapClasses(self,hapDef):
		muts = hapDef[0]
		invs = hapDef[1]
		potMutAlleles = []
		for m in muts:
			potMutAlleles += [([],[m])]
		potMutHaps = [list(itertools.chain(*x)) for x in itertools.product(*potMutAlleles)]
		potInvAlleles = []
		for i in invs:
			potInvAlleles += [([],[i])]
		potInvHaps = [list(itertools.chain(*x)) for x in itertools.product(*potInvAlleles)]
		hapClasses = list(itertools.product(potMutHaps,potInvHaps))
		orderedHapClasses = [[] for x in hapClasses]
		for hap in hapClasses:
			orderedHapClasses[self.__indexFromHap(hap,hapDef)] = hap
		return orderedHapClasses


	# Takes a filename and haplotype definition ID
	# Writes a tab delineated file of the generation, count, 
	#  average number of mutations in buffer, average cumulative survival, 
	#  and average cumulative effect data for that inversion
	def writeHaplotypes(self,filename,ID,fillRecord=False):
		# Exit if there is no record
		if ID not in self.record[12] or len(self.record[12][ID][0]) == 0:
			return False
		if self.verbose:
			print(self.record[12])
		firstGen = self.record[11][ID][2]
		firstRec = self.record[11][ID][3]
		hapDef = self.record[11][ID][0:2]
		hapClasses = self.__buildHapClasses(hapDef)
		hapClassRange = np.arange(len(hapClasses))
		hapClassNames = []
		for hap in hapClasses:
			hapClassNames += [':'.join([str(m) for m in hap[0]]) + ';' + ':'.join([str(i) for i in hap[1]])]
		# Write the allele, homozygote, heterozygote counts, with haplotype class
		outfile = open(filename, 'w')
		header = 'Generation\tLifeHistoryStage\tHaplotype\tCount\tHomozygote\tHeterozygote\t'+\
			'FemaleCount\tFemaleHomozygote\tFemaleHeterozygote\tMaleCount\tMaleHomozygote\tMaleHeterozygote\t'+\
			'OffspringCount\tOffspringCountInHomozygote\tOffspringCountInHeterozygote\t'+\
			'OffspringCountInFemale\tOffspringCountInFemaleHomozygote\tOffspringCountInFemaleHeterozygote\t'+\
			'OffspringCountInMale\tOffspringCountInMaleHomozygote\tOffspringCountInMaleHeterozygote\n'
		outfile.write(header)
		t = 0
		recordGen = self.record[0][t]
		recordStage = self.record[10][t]
		# Padding with 0's ?
		while t < firstRec:
			if fillRecord:
				for hapName in hapClassNames:
					outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(hapName) + '\t0'*18+'\n')
			t += 1
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
		offset = t
		for t in range(t,len(self.record[0])):
			m = t - offset
			recordGen = self.record[0][t]
			recordStage = self.record[10][t]
			for h in hapClassRange:
				# relies on the hapClassNames and record[12] being in the same order
				# hapIndex = self.__indexFromHap(hapClasses[h],hapDef)
				# hapIndex = h
				hapName = hapClassNames[h]
				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(hapName) + \
					'\t' + str(self.record[12][ID][0][h][m]) + '\t' + str(self.record[12][ID][1][h][m]) + \
					'\t' + str(self.record[12][ID][2][h][m]) + '\t' + str(self.record[12][ID][3][h][m]) + \
					'\t' + str(self.record[12][ID][4][h][m]) + '\t' + str(self.record[12][ID][5][h][m]) + \
					'\t' + str(self.record[12][ID][6][h][m]) + '\t' + str(self.record[12][ID][7][h][m]) + \
					'\t' + str(self.record[12][ID][8][h][m]) + '\t' + str(self.record[12][ID][8][h][m]) + \
					'\t' + str(self.record[12][ID][10][h][m]) + '\t' + str(self.record[12][ID][11][h][m]) + \
					'\t' + str(self.record[12][ID][12][h][m]) + '\t' + str(self.record[12][ID][13][h][m]) + \
					'\t' + str(self.record[12][ID][14][h][m]) + '\t' + str(self.record[12][ID][15][h][m]) + \
					'\t' + str(self.record[12][ID][16][h][m]) + '\t' + str(self.record[12][ID][17][h][m]) + '\n')

		# # FIXATION OF HAPLOTYPES NOT IMPLEMENTED HERE, PADDED IN RECORD
		# # - consider storing the fixed haplotype in record[11]
		# if fillRecord:
		# 	for t in range(t,len(self.recoaasdrd[0])):
		# 		recordGen = self.record[0][t]
		# 		recordStage = self.record[10][t]
		# 		if self.__locFixed[ID] and not self.__mutLost[ID]:
		# 			for hap in hapClassNames:
		# 				outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(hapName) +\
		# 					'\t' + str(2*self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0' + \
		# 					'\t' + str(2*self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0' + \
		# 					'\t' + str(2*self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0' + \
		# 					'\t' + str(self.record[5][t]) + '\t' + str(self.record[5][t]) + '\t0' + \
		# 					'\t' + str(self.record[13][t]) + '\t' + str(self.record[13][t]) + '\t0' + \
		# 					'\t' + str(self.record[14][t]) + '\t' + str(self.record[14][t]) + '\t0' + '\n')
		# 		else:
		# 			outfile.write(str(recordGen) + '\t' + str(recordStage) + '\t' + str(hapName) + '\t0'*18+'\n')
		outfile.close()
		return True



	# Takes a filename and writes a tab delineated file of the mutation frequency 
	#   by generation and ID with generation in the first column
	def writeAllMutFreqTable(self,filename):
		outfile = open(filename, 'w')
		header = 'Generation\tLifeHistoryStage\tPopSize'
		mutID = []
		mutInit = []
		mutFinal = []
		offsets = []
		recLens = []
		for m in range(self.__mutIDcount):
			if m in self.record[2]:
				recLen = len(self.record[2][m][0])
				if recLen > 0:
					mutID += [m]
					# header += '\t'+str(m)
					header += '\tM'+str(m)
					initGen = self.record[1][m][4]
					finalGen = self.record[1][m][5]
					mutInit += [initGen]
					if initGen > finalGen:
						mutFinal += [self.age+1]
					else:
						mutFinal += [finalGen]
					offsets += [-1]
					recLens += [recLen]
		# print(mutID)
		# print(mutInit)
		# print(mutFinal)
		# print(self.record[1])
		outfile.write(header + '\n')
		for g in range(len(self.record[0])):
			gen = self.record[0][g]
			lhs = self.record[10][g]
			size = self.record[5][g]
			genLine = str(gen) + '\t' + str(lhs) + '\t' + str(size)
			for m in range(len(mutID)):
				if gen < mutInit[m]:
					genLine += '\t0'
				elif offsets[m] < 0: # for gen < mutFinal[i], The last generation isn't recorded as it is removed first
					offsets[m] = g
					genLine += '\t' + str(self.record[2][mutID[m]][0][0])
				elif g-offsets[m] < recLens[m]:
					genLine += '\t' + str(self.record[2][mutID[m]][0][g-offsets[m]])
				elif self.__locFixed[mutID[m]] and not self.__mutLost[mutID[m]]:
					genLine += '\t' + str(2*self.record[5][g])
				else:
					genLine += '\t0'
			outfile.write(genLine + '\n')
		outfile.close()
		return

	# Takes a filename and writes a tab delineated file of the inversion frequency
	#   by generation and ID with generation in the first column
	def writeAllInvFreqTable(self,filename):
		outfile = open(filename, 'w')
		# header = 'Generation'
		header = 'Generation\tLifeHistoryStage\tPopSize\tStd'
		invID = []
		invInit = []
		invFinal = []
		offsets = []
		recLens = []
		for i in range(self.__invIDcount):
			if i in self.record[4]:
				recLen = len(self.record[4][i][0])
				if recLen > 0:
					invID += [i]
					# header += '\t'+str(i)
					header += '\tI'+str(i)
					initGen = self.record[3][i][3]
					finalGen = self.record[3][i][4]
					invInit += [initGen]
					if initGen > finalGen:
						invFinal += [self.age+1]
					else:
						invFinal += [finalGen]
					offsets += [-1]
					recLens += [recLen]
		# print(invInit)
		# print(invFinal)
		# print(self.record[3])
		outfile.write(header + '\n')
		for g in range(len(self.record[0])):
			gen = self.record[0][g]
			lhs = self.record[10][g]
			size = self.record[5][g]
			genLine = str(gen) + '\t' + str(lhs) + '\t' + str(size)
			stdFreq = self.record[4]["Std"][0][g]
			genLine += '\t' + str(stdFreq)
			for i in range(len(invID)):
				if gen < invInit[i]:
					genLine += '\t0'
				elif offsets[i] < 0:
					offsets[i] = g
					genLine += '\t' + str(self.record[4][invID[i]][0][0])
				elif g-offsets[i] < recLens[i]:
					# print("g "+str(g)+" gen "+str(gen)+" Offset "+str(offsets[i])+" Index "+str(g-offsets[i])+" Init "+str(invInit[i])+" Final "+str(invFinal[i]))
					# print(len(self.record[4][i]))
					# print(self.record[4][i][g-offsets[i]-2:])
					genLine += '\t' + str(self.record[4][invID[i]][0][g-offsets[i]])
				elif self.__invFixed[invID[i]] and not self.__invLost[invID[i]]:
					genLine += '\t' + str(2*self.record[5][g])
				else:
					genLine += '\t0'
			outfile.write(genLine + '\n')
		outfile.close()
		return

	# Takes a filename and writes a tab delineated file of the position, effect, chromosome,
	#   and initial generation data for all mutations
	def writeMutCharTable(self,filename):
		outfile = open(filename, 'w')
		# outfile.write('Position\tSurEffect\tRepEffect\tChromosome\tInitGen\tFinalGen\n')
		outfile.write('ID\tPos\tSurEf\tRepEf\tChrom\tInitGen\tFinalGen\n')
		for m in range(len(self.record[1])):
			line = repr(m)+'\t'
			for datum in self.record[1][m]:
				line += repr(datum) + '\t'
			outfile.write(line[:-1] + '\n')
		outfile.close()
		return

	# Takes a filename and writes a tab delineated file of the positions, chromosome,
	#   and initial generation data for all inversions
	def writeInvCharTable(self,filename):
		outfile = open(filename, 'w')
		# outfile.write('Position1\tPosition2\tChromosome\tInitGen\tFinalGen\n')
		outfile.write('ID\tPos1\tPos2\tChrom\tInitGen\tFinalGen\n')
		for i in range(len(self.record[3])):
			line = repr(i)+'\t'
			for datum in self.record[3][i]:
				line += repr(datum) + '\t'
			outfile.write(line[:-1] + '\n')
		outfile.close()
		return

	# Takes a filename and writes a tab delineated file of the haplotype definition
	#   and initial generation data for all recorded haplotype classifications
	def writeHapCharTable(self,filename):
		outfile = open(filename, 'w')
		outfile.write('ID\tDerivedHaplotype\tInitialGen\n')
		for i in range(len(self.record[11])):
			muts = self.record[11][i][0]
			invs = self.record[11][i][1]
			derivedHapClassName = ':'.join([str(m) for m in muts]) + ';' + ':'.join([str(i) for i in invs])
			line = repr(i)+'\t'+derivedHapClassName+'\t'+repr(self.record[11][i][2])+'\n'
			outfile.write(line)
		outfile.close()
		return



	# Takes a filename and writes a tab delineated file of the population statistics
	#   by generation in the first column
	def writePopStatTable(self,filename):
		with open(filename, 'w') as outfile:
			outfile.write('Generation\tLifeHistoryStage\tSurMean\tSurVar\tRepMean\tRepVar\tPopSize\tNumFemales\tNumMales\n')
			for g in range(len(self.record[0])):
				outfile.write(repr(self.record[0][g]) + '\t' + str(self.record[10][g]) + '\t' +\
					repr(self.record[6][g]) + '\t' + repr(self.record[7][g]) + '\t' + \
					repr(self.record[8][g]) + '\t' + repr(self.record[9][g]) + '\t' + \
					repr(self.record[5][g]) + '\t' + repr(self.record[13][g]) + '\t' + \
					repr(self.record[14][g]) + '\n')
		return

	# Takes a filename and writes a tab delineated file of the distribution of offspring
	#   by generation in the first column
	def writeRepDistTable(self,filename):
		with open(filename, 'w') as outfile:
			header = 'Generation\tSex'
			# CURRENTLY UNCLEAR HOW TO REPRESENT THIS UNDER CHANGING POPULATION SIZES
			maxPopSize = self.size
			recordedDistIndexes = []
			for g in range(len(self.record[0])):
				if self.record[15][g] is not None:
					maxPopSize = max(maxPopSize,self.record[5][g])
					recordedDistIndexes += [g]
			for n in np.arange(maxPopSize+1):
				header += '\t'+str(n)
			header += '\n'
			outfile.write(header)
			for g in recordedDistIndexes:
				femLine = repr(self.record[0][g]) + '\tFemale'
				for n in self.record[15][g]:
					femLine += '\t'+str(int(n))
				femLine += '\n'
				outfile.write(femLine)
				malLine = repr(self.record[0][g]) + '\tMale'
				for n in self.record[16][g]:
					malLine += '\t'+str(int(n))
				malLine += '\n'
				outfile.write(malLine)
		return


	# Writes a summary file for the parameters of the simulation run
	# UPDATE - changeable size and additional parameters
	def writeParamSummary(self,filename):
		outfile = open(filename, 'w')
		outfile.write('Parameter\tValue\n')
		outfile.write('Size\t'+str(self.size)+'\n')
		outfile.write('NumChromosomes\t'+str(self.numChrom)+'\n')
		outfile.write('MutationRate\t'+str(self.mutRate)+'\n')
		outfile.write('InversionMutationRate\t'+str(self.mutRateInv)+'\n')
		outfile.write('MutationEffectOffsetSD\t'+str(self.mutEffectDiffSD)+'\n')
		outfile.write('MinimumInversionLength\t'+str(self.minInvLen)+'\n')
		outfile.write('ConversionRate\t'+str(self.conversionRate)+'\n')
		outfile.write('RecombinationRate\t'+str(self.recombRate)+'\n')
		outfile.write('EncounterNumber\t'+str(self.encounterNum)+'\n')
		outfile.write('FemaleChoiceNoiseSD\t'+str(self.choiceNoiseSD)+'\n')
		outfile.write('InversionRecordBuffer\t'+str(self.invRecBuffer)+'\n')
		outfile.write('RandomSex\t'+str(self.randomSex)+'\n')
		outfile.write('WillMutate\t'+str(self.willMutate)+'\n')
		outfile.write('WillMutInv\t'+str(self.willMutInv)+'\n')
		outfile.write('WillConvert\t'+str(self.willConvert)+'\n')
		outfile.write('WillRecombine\t'+str(self.willRecombine)+'\n')
		outfile.write('NoMaleCost\t'+str(self.noMaleCost)+'\n')
		outfile.write('NoFemaleCost\t'+str(self.noFemaleCost)+'\n')
		outfile.write('RemoveFixedVariants\t'+str(self.removeFixed)+'\n')
		outfile.write('SampleSurvivingParentsPerDrawnOffspring\t'+str(self.sampleSurvivorsPerOffspring)+'\n')
		outfile.write('SampleMaleEncountersWithReplacement\t'+str(self.sampleEncountersWithReplacement)+'\n')
		outfile.write('RecordPoolOfReproductiveAdults\t'+str(self.recordReproductiveAdults)+'\n')
		outfile.write('RecordUserDefinedHaplotypes\t'+str(self.recordHaps)+'\n')
		outfile.write('RecordUserDefinedHaplotypes\t'+str(self.recordHaps)+'\n')
		# self.minNumGensKeepRecord
		# self.minFreqKeepMutRecord
		# self.minFreqKeepInvRecord
		# self.verbose
		# self.testing
		# self.storePickledPop
		outfile.close()
		return


	# Writes a summary for the parameters of the simulation run, to the console
	# UPDATE - changeable size and additional parameters
	def printParamSummary(self):
		print('Parameter\t\tValue')
		print('Size\t\t'+str(self.size))
		print('NumChromosomes\t'+str(self.numChrom))
		print('MutationRate\t'+str(self.mutRate))
		print('InversionMutationRate\t'+str(self.mutRateInv))
		print('MutationEffectOffsetSD\t'+str(self.mutEffectDiffSD))
		print('MinimumInversionLength\t'+str(self.minInvLen))
		print('ConversionRate\t'+str(self.conversionRate))
		print('RecombinationRate\t'+str(self.recombRate))
		print('EncounterNumber\t'+str(self.encounterNum))
		print('FemaleChoiceNoiseSD\t'+str(self.choiceNoiseSD))
		print('InversionRecordBuffer\t'+str(self.invRecBuffer))
		print('RandomSex\t'+str(self.randomSex))
		print('WillMutate\t'+str(self.willMutate))
		print('WillMutInv\t'+str(self.willMutInv))
		print('WillConvert\t'+str(self.willConvert))
		print('WillRecombine\t'+str(self.willRecombine))
		print('NoMaleCost\t'+str(self.noMaleCost))
		print('NoFemaleCost\t'+str(self.noFemaleCost))
		print('RemoveFixedVariants\t'+str(self.removeFixed))
		print('SampleSurvivingParentsPerDrawnOffspring\t'+str(self.sampleSurvivorsPerOffspring))
		print('SampleMaleEncountersWithReplacement\t'+str(self.sampleEncountersWithReplacement))
		print('RecordPoolOfReproductiveAdults\t'+str(self.recordReproductiveAdults))
		print('RecordUserDefinedHaplotypes\t'+str(self.recordHaps))
		print('RecordUserDefinedHaplotypes\t'+str(self.recordHaps))
		# self.minNumGensKeepRecord
		# self.minFreqKeepMutRecord
		# self.minFreqKeepInvRecord
		# self.verbose
		# self.testing
		# self.storePickledPop
		return

	# Writes record to outfiles as such:
	# One summary document
	# Tab delineated table of mutation characteristics for all mutations
	# Tab delineated table of inversion characteristics for all inversions
	# Tab delineated table with generation and count data for every mutation
	# Tab delineated table with generation, count, etc. additional data for every inversion
	def writeRecordTables(self,outFilePrefix,fillRecord=False,
			writeIndMuts=True,writeIndInvs=True,
			writeSummMuts=True,writeSummInvs=True,
			writeRepDist=True):
		self.writeParamSummary(outFilePrefix+'ParamSumm.txt')
		if writeSummMuts:
			self.writeMutCharTable(outFilePrefix+'MutSumm.txt')
			self.writeAllMutFreqTable(outFilePrefix+'MutFreqs.txt')
		if writeSummInvs:
			self.writeInvCharTable(outFilePrefix+'InvSumm.txt')
			self.writeAllInvFreqTable(outFilePrefix+'InvFreqs.txt')
		if len(self.record[11])>0:
			self.writeHapCharTable(outFilePrefix+'HapSumm.txt')
		self.writePopStatTable(outFilePrefix+'PopStats.txt')
		if self.recordReproductiveAdults and writeRepDist:
			filename = outFilePrefix+'RepDist.txt'
			self.writeRepDistTable(filename)

		if writeIndMuts:
			for mutID in range(self.__mutIDcount):
				filename = outFilePrefix+'Mut'+str(mutID)+'.txt'
				self.writeMutation(filename,mutID,fillRecord)
		if writeIndInvs:
			for invID in range(self.__invIDcount):
				filename = outFilePrefix+'Inv'+str(invID)+'.txt'
				self.writeInversion(filename,invID,fillRecord)
			filename = outFilePrefix+'Inv'+'Std'+'.txt'
			self.writeStdArrangement(filename)
		for hapID in range(len(self.record[11])):
			filename = outFilePrefix+'Hap'+str(hapID)+'.txt'
			self.writeHaplotypes(filename,hapID,fillRecord)

		if self.storePickledPop:
			filename = outFilePrefix+'Pop.pickle'
			self.writePopPickle(filename)
		return



# For running a simulation from a model file
def runModelFile(param_file,output_file_path):
    # Verify that the input and output filenames are provided
    assert os.path.isfile(param_file),\
        'Model specification file argument (--model {}) \
        is not a file path'.format(param_file)
    assert os.path.isfile(output_file_path),\
        'Output file path argument (--out {}) \
        is not a file path or is not given'.format(param_file)


# # For handling direct calls to SAIsim.py, for easier definition of models from input files.
# def main(args):
# 	param_file = args.model
# 	output_file_path = args.out
# 	if param_file != '':
# 		runModelFile(param_file,output_file_path)
# 	else:
# 	    # Parse input arguments
# 	    size = args.size
# 		mutRate = mutRate
# 		mutRateInv = mutRateInv
# 		#
# 		mutEffectDiffSD = mutEffectDiffSD
# 		minInvLen = minInvLen
# 		conversionRate = conversionRate
# 		recombRate = recombRate
# 		encounterNum = encounterNum
# 		choiceNoiseSD = choiceNoiseSD
# 		# For modeling D.mel specific recombination/other biology
# 		maleAchiasmy = maleAchiasmy
# 		randomSex = randomSex
# 		numChrom = numChrom
# 		# For modeling equlibrium dynamics without further mutation or other mechanics
# 		willMutate = willMutate
# 		willMutInv = willMutInv
# 		willConvert = willConvert
# 		willRecombine = willRecombine
# 		# For modeling survival and reproductive effects differently
# 		noMaleCost = noMaleCost
# 		# Record keeping variables, record data structure defined in __updateRecord
# 		invRecBuffer = invRecBuffer = args.dataset
# 	    bin_num = args.bin_num
# 	    output_file_path = args.out
# 	    bin_str = args.bin_str.lower()

# 	    assert (size>1), 'size argument (--size) {} \
# 	        must be larger than 1'.format(size)
# 	    assert (output_file_path!=''),'MI output file path (--out) {} \
# 	        not given'.format(output_file_path)

# 	    # Run the simulation


# # Runs if SAIsim.py is the main script
# if __name__ == "__main__":
#     # Collect input arguments
#     # CHECK argparse documentation for other input ideas (positional arguments? other?)
#     parser = argparse.ArgumentParser(
#         description=__doc__,
#         formatter_class=argparse.RawDescriptionHelpFormatter)
#     parser.add_argument('--model',
#                         help='parameter and epoch file, see documentation',
#                         type=str,
#                         default='')
#     parser.add_argument('--out',
#                         help='output path and file prefix',
#                         type=str,
#                         default='')

#     args = parser.parse_args()

#     main(args)




