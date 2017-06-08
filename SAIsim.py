
# The Object Oriented implementation of Sexually Antagonistic mutation accumulation
#     in Inversions in a forward Population simulation
# Ideally, models large inversion polymorphisms in populations with high reproductive skew


# REWRITE WITH NUMPY ARRAYS/VECTORS?

import numpy as np

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
			isFly, willConvert, willRecombine, genome = [[[[],[]],[[],[]]]]):
		# super(individual, self).__init__()
		self.sex = sex
		self.mutEffectDiffSD = mutEffectDiffSD
		self.recombRate = recombRate
		self.conversionRate = conversionRate
		self.minInvLen = minInvLen
		self.lenChrom = lenChrom
		# CONSIDER JUST PASSING VARIABLES WHEN NEEDED
		self.isFly = isFly
		self.willConvert = willConvert
		self.willRecombine = willRecombine
		# genome = [[[chr1hom1],[chr1hom2]],[[chr2hom1],[chr2hom2]],..]
		# where chomosome homologs = [mutist, InvList]
		self.genome = genome


	# For getting the insertion index of a position from a list of [position,..] lists
	def __getInsInd(self,mutInvList,mutInvPos):
		i = 0
		while (i < len(mutInvList)) and (mutInvList[i][0] <= mutInvPos):
			i += 1
		return i

	# For inserting a mutation/inversion into a corresponding list
	# May rely on passing the list as a reference (doesn't currently)
	def __insert(self,mutInvList,mutInv):
		i = 0
		while (i < len(mutInvList)) and (mutInvList[i][0] <= mutInv[0]):
			i += 1
		mutInvList[i:i] = [mutInv]
		return mutInvList


	# generates mutations
	# MAKE THIS BIOLOGICALLY RELEVANT? SETTLE ON A DISTRIBUTION - lognormal?
	def __genEffectSizes(self):
		base = np.random.ranf()
		offset = np.random.normal(scale=self.mutEffectDiffSD)
		mutEffects = []
		while (min(1-base,base) < offset) or (offset < max(base-1,-base)):
			offset = np.random.normal(scale=self.mutEffectDiffSD)

		# # For modeling effects multiplicatively?
		# return [base+offset,base+offset]

		# Survival mult., Quality additive
		return [1-base+offset,base+offset]

	# Generates a mutation of the form [position, survival effect, rep effect, ID]
	def mutate(self,ID):
		# print(self.genome)
		mutPos = np.random.ranf()*self.lenChrom
		# mutPos = np.random.ranf()
		mutEffects = self.__genEffectSizes()
		mutation = [mutPos]+mutEffects+[ID]
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		# print(self.genome[chromIndex][homIndex][0])
		self.genome[chromIndex][homIndex][0] = self.__insert(self.genome[chromIndex][homIndex][0],mutation)
		# print(self.genome[chromIndex][homIndex][0])
		# Return the mutation data for recording
		return mutation[0:3]+[chromIndex]

	# Generates and inserts an inversion into the genome, unless there is no open region >= minInvLen
	# Inversions cover the region [pos1,pos2)
	# Does not use resampling 
	def mutateInv(self,ID):
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		chromHomInv = self.genome[chromIndex][homIndex][1]
		# Pick where the insertion may be placed
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
		# If there is no space for a new inversion >= minInvLen, don't add one, return None
		if len(openRegStarts) == 0:
			return
		# Now sample which open region to add the inversion to
		regChoice = np.random.choice(len(openRegStarts),p=[l/sum(openRegLengths) for l in openRegLengths])
		# Generate the inversion
		posA = openRegStarts[regChoice] + np.random.ranf()*openRegLengths[regChoice]
		posB = openRegStarts[regChoice] + np.random.ranf()*openRegLengths[regChoice]
		inversion = [min(posA,posB),max(posA,posB),ID]
		index = openRegIndexes[regChoice]
		chromHomInv[index:index] = [inversion]
		# Return the inversion data to record
		return inversion[0:2]+[chromIndex]

	# # For updating the repQuality and survival if pre-calculated
	# # Not currently used (obviously)
	# def updatePhenotypes(self):
	# 	return

	# Negative survival effects are independent and multiplicative
	def survival(self):
		survival = 1.0
		for chrom in self.genome:
			for mut in chrom[0][0]:
				survival *= mut[1]
			for mut in chrom[1][0]:
				survival *= mut[1]
		return survival

	# FIX THIS - what quality accounting process is biological? (This doesn't seem tooo bad, tbh)
	# ALSO - maybe calculate it once to save time? probably recalculated several times in choosing fathers
	def repQuality(self):
		quality = 1.0
		for chrom in self.genome:
			for mut in chrom[0][0]:
				quality += mut[2]
				# For modeling effects multiplicatively? - requires changing effect generation
				# quality *= mut[2]
			for mut in chrom[1][0]:
				quality += mut[2]
				# For modeling effects multiplicatively? - requires changing effect generation
				# quality *= mut[2]
		return quality

	# Takes a chromosome and returns the chromosome with converted mutations
	# May want to allow for a boolean hasConversion to turn this off in the simulation
	def __convertChrom(self,chrom):
		homInd1 = 0
		homInd2 = 0
		# Check for heterozygosity and add/remove mutations following conversionRate probability
		# Relies on the invariant that earlier indexed mutations have earlier position
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
						chrom[1][0][homInd2-1:homInd2-1] = [mut1]
						homInd2 += 1
						homInd1 += 1
				else:
					homInd1 += 1
			if mut1[0] > mut2[0]:
				if np.random.ranf() < self.conversionRate:
					if np.random.randint(2):
						chrom[1][0][homInd2:homInd2+1] = []
					else:
						chrom[0][0][homInd1-1:homInd1-1] = [mut2]
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
		return chrom


	# Takes two inversion lists for homologous chromosomes,
	# indexes of the closest inversions to start <= the position of interest,
	# and the position of interest
	# Returns the position following the position of interest such that an odd # of crossovers
	# between the two would generate aneuploid gametes, returns -1 if it isn't in such a region
	def __getAneupRegion(self,invHom1,invHom2,ind1,ind2,recPos):
		# Deal with no inversions present on the chromosome
		if len(invHom1) == 0:
			if (len(invHom2) == 0) or (recPos >= invHom2[ind2][1]):
				return -1
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

	# Model independently per chromosome, currently doesn't redraw if recombination events fail,
	# just throws out all recombination events in an aneuploidy-causing region if they are odd in number
	def genGamete(self):
		gamGenome = []
		for parChrom in self.genome:
			numRecomb = np.random.poisson(self.recombRate)
			# print('Recombinations: '+str(numRecomb))
			currHom1 = bool(np.random.randint(2))
			# This first copy step is necessary independent of conversion/recombination activity, 
			# to prevent changes to shared-referent chromosomes
			chrom = [[list(parChrom[0][0]),list(parChrom[0][1])],[list(parChrom[1][0]),list(parChrom[1][1])]]
			# Perform fly-specific non-recombination non-conversion check for males
			if self.willConvert and (self.sex == 'F' or (not self.isFly)):
				chrom = self.__convertChrom(chrom)
			if (not self.willRecombine) or numRecomb == 0 or (self.isFly and self.sex == 'M'):
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
def __genSexes(size, sexes, randomSex):
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
	if numMales == 0:
		newSexes += ['M']
		numMales += 1
	if numFemales == 0:
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
	sexes += newSexes
	# Mutates the sexes list but will be returned anyway
	return sexes

# For generating a set of genomes, sexes, and a record from specified mutations and inversions, 
# along with partial genomes and sexes lists to specify specific genomes,
# input lists formatted as described in self.record but with count instead of initial generation:
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
		lenChrom = 1.0, genomes = [], sexes = []):
	# Handling input size checking
	(numGenomes,numChrom) = __inputSizeCheck(size, numChrom, genomes, sexes) # UPDATE FOR NUMCHROM FROM MUT/INV LIST

	# Generate the remaining sexes list
	sexes = __genSexes(size, sexes, randomSex)

	# Generate the remaining genomes
	numNewGenomes = size-numGenomes
	newGenomes = [[[[[],[]],[[],[]]] for i in range(numChrom)] for j in range(numNewGenomes)]
	# Will need to have separate lists per chromosome within these
	posOrderedMut = [[] for i in range(numChrom)]
	posOrderedInv = [[] for i in range(numChrom)]
	# Record should be returned with only [1,3] actually populated, [2,4..8] with empty lists 
	#   for 0th generation update
	record = [[],[],[],[],[],[],[],[],[]]

	numMut = len(mutList)
	for m in range(numMut):
		count = mutList[m][4]
		if count > numNewGenomes:
			raise SAIpop.InputError(mutList,'Mutation counts must be less than the remaining population size')
		mutation = mutList[m][0:3]+[m]
		chrom = mutList[m][3]
		record[1] += [mutList[m][0:4]+[0,-1]]
		record[2] += [[]]
		# Insert the mutation and count to the ordered list for addition to the genomes
		i = 0
		while (i < len(posOrderedMut[chrom])) and (posOrderedMut[chrom][i][0][0] < mutation[0]):
			i += 1
		posOrderedMut[chrom][i:i] = [(mutation,count)]
	numInv = len(invList)
	for i in range(numInv):
		count = invList[m][3]
		if count > numNewGenomes:
			raise SAIpop.InputError(invList,'Inversion counts must be less than the remaining population size')
		inversion = invList[i][0:2]+[i]
		chrom = invList[i][2]
		record[3] += [invList[i][0:3]+[0,-1]]
		record[4] += [[]]
		record[5] += [[]]
		record[6] += [[]]
		record[7] += [[]]
		# Insert the inversion and count to the ordered list for addition to the genomes
		i = 0
		while (i < len(posOrderedInv[chrom])) and (posOrderedInv[chrom][i][0][1] <= inversion[0]):
			i += 1
		if i < len(posOrderedInv[chrom]) and posOrderedInv[chrom][i+1][0][0] <= inversion[1]:
			raise SAIpop.InputError(invList,'Inversions cannot overlap')
		posOrderedInv[chrom][i:i] = [(inversion,count)]

	# Preserve order by adding the mutations/inversions from in (position) order lists
	# Need to ensure no double-placement at any position
	possibleSpots = []
	for i in range(numNewGenomes):
		for h in [0,1]:
			possibleSpots += [(i,h)]
	for chrom in range(numChrom):
		for (mut,count) in posOrderedMut[chrom]:
			np.random.shuffle(possibleSpots)
			for (i,hom) in possibleSpots[0:count]:
				newGenomes[i][chrom][hom][0] += [mut]
		for (inv,count) in posOrderedInv[chrom]:
			np.random.shuffle(possibleSpots)
			for (i,hom) in possibleSpots[0:count]:
				newGenomes[i][chrom][hom][1] += [inv]
	genomes += newGenomes
	return (genomes,sexes,record)


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
		lenChrom = 1.0, genomes = [], sexes = []):
	# Handling input size checking
	(numGenomes,numChrom) = __inputSizeCheck(size, numChrom, genomes, sexes)
	# Check that the haplotype chromosome size matches the pre-specified genomes
	if (len(genomes) > 0) and (len(genomes[0]) != len(hapList[0][0])):
		# Ignore numChrom if genomes pre-specified, peg it to the number of chormosomes in genomes[0]
		raise SAIpop.InputError(hapList,'Haplotype lists must have chromosome number matching prespecified genomes')

	# Generate the remaining sexes list
	sexes = __genSexes(size, sexes, randomSex)

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
	record = [[],[],[],[],[],[],[],[],[]]
	for m in range(len(mutList)):
		record[1] += [mutList[m][0:4]+[0,-1]]
		record[2] += [[]]
	for i in range(len(invList)):
		record[3] += [invList[i][0:3]+[0,-1]]
		record[4] += [[]]
		record[5] += [[]]
		record[6] += [[]]
		record[7] += [[]]

	return (genomes,sexes,record)

		
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
			isFly = True, randomSex = True, numChrom = 1, lenChrom = 1.0, willMutate = True,
			willMutInv = True, willConvert = True, willRecombine = True, noMaleCost = False,
			genomes = [], sexes = [], record = [[],[],[],[],[],[],[],[],[]]):
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
		self.minInvLen = minInvLen
		self.conversionRate = conversionRate
		self.recombRate = recombRate
		self.encounterNum = encounterNum
		self.choiceNoiseSD = choiceNoiseSD
		self.males = []
		self.females = []
		# For modeling D.mel specific recombination/other biology
		self.isFly = isFly
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
		# Record keeping variables, record data structure defined in __updateRecord
		self.invRecBuffer = invRecBuffer
		self.__invFixed = [False for i in range(self.__invIDcount)]
		self.__mutFixed = [False for m in range(self.__mutIDcount)]
		self.record = record
		self.age = 0
		# Handling input and default populations
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
							minInvLen,lenChrom,isFly,willConvert,willRecombine,genomes[i])]
					else:
						self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,\
							minInvLen,lenChrom,isFly,willConvert,willRecombine,genomes[i])]
			numLeft = size-numGenomes
			if len(self.females) < 1:
				numLeft -= 1
				self.females += [individual('F',mutEffectDiffSD,recombRate,conversionRate,\
					minInvLen,lenChrom,isFly,willConvert,willRecombine,\
					[[[[],[]],[[],[]]] for i in range(self.numChrom)])]
			if len(self.males) < 1:
				numLeft -= 1
				self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,\
					minInvLen,lenChrom,isFly,willConvert,willRecombine,\
					[[[[],[]],[[],[]]] for i in range(self.numChrom)])]
			# Generate the remaining set of 'size' individuals
			#   with random or equal probability sex and no mutations or inversions
			if randomSex:
				numMales = np.random.binomial(numLeft,0.5)
			else:
				numMales = int(numLeft/2)
			self.males += [individual('M',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
				lenChrom,isFly,willConvert,willRecombine,\
				[[[[],[]],[[],[]]] for i in range(self.numChrom)]) for i in range(numMales)]
			self.females += [individual('F',mutEffectDiffSD,recombRate,conversionRate,minInvLen,\
				lenChrom,isFly,willConvert,willRecombine,\
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


	# For changing the population size during simulation
	def setSize(self,newSize):
		self.size = newSize
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
	def __removeFixed(self):
		wholePop = self.males + self.females
		mutCounts = [0]*self.__mutIDcount
		mutPos = [[] for i in range(self.__mutIDcount)]
		invCounts = [0]*self.__invIDcount
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
		# Remove fixed inversions
		mutationRemoved = False
		for ID in range(self.__mutIDcount):
			count = mutCounts[ID]
			# print("Mut "+str(ID)+" Count "+str(count))
			if count == 2*self.size:
				mutationRemoved = True
				# print("Mutation Removed, mutFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(mutCounts[ID]))
				# print(str(len(mutPos[ID]))+" mutPus entries")
				# Update the fix/loss record
				self.__mutFixed[ID] = True
				self.record[1][ID][5] = self.age
				for [i,c,h,mutIndex] in mutPos[ID]:
					# Remove the mutation
					# print("Removal Attempt")
					# print(wholePop[i].genome[c][h][0])
					wholePop[i].genome[c][h][0][mutIndex:mutIndex+1] = []
					# print(wholePop[i].genome[c][h][0])
				testMutCounts = [0]*self.__mutIDcount
				for i in range(len(wholePop)):
					for c in range(self.numChrom):
						for h in range(2):
							chromHomMut = wholePop[i].genome[c][h][0]
							for mutIndex in range(len(chromHomMut)):
								mutation = chromHomMut[mutIndex]
								testMutCounts[mutation[3]] += 1
				# print(str(testMutCounts[ID])+" remaining after removal")
			elif count == 0:
				# Update the fix/loss record
				if self.record[1][ID][5] < self.record[1][ID][4]: # Change to use another flag? mutLost?
					# self.__mutLost[ID] = True
					self.record[1][ID][5] = self.age
			elif count > 2*self.size:
				print("Error: more mutations than alleles possible for ID "+str(ID)+" Count "+str(count))
				# self.printGenomes()
		# All mutations are passed as references, so don't want to change more than one
		# Also, no mutations will be in multiple fixed inversions, so the scope is fine :P
		firstMutEncounter = [True]*self.__mutIDcount
		inversionRemoved = False
		# Remove fixed inversions
		inversionRemoved = False
		for ID in range(self.__invIDcount):
			if invCounts[ID] == 2*self.size:
				inversionRemoved = True
				# print("Inversion Removed, invFixed["+str(ID)+"]=True at gen "+str(self.age)+" count "+str(invCounts[ID]))
				# Update the fix/loss record
				self.__invFixed[ID] = True
				self.record[3][ID][4] = self.age
				inversion = self.record[3][ID]
				length = self.record[3][ID][1] - self.record[3][ID][0]
				invCenter = self.record[3][ID][0] + length/2.0
				for [i,c,h,invIndex] in invPos[ID]:
					# Remove the inversion
					# print("Removal Attempt")
					# print(wholePop[i].genome[c][h][1])
					wholePop[i].genome[c][h][1][invIndex:invIndex+1] = []
					# print(wholePop[i].genome[c][h][1])
					# Rearrange the mutation positions to reflect the loss of inversion data
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
							# Flip the mutation across the center
							mut[0] = 2*invCenter-mut[0]
						# print(mut)
			elif invCounts[ID] == 0:
				# Update the fix/loss record
				if self.record[3][ID][4] < self.record[3][ID][3]: # Change to use another flag? invLost?
					# self.__invLost[ID] = True
					self.record[3][ID][4] = self.age
			elif invCounts[ID] > 2*self.size:
				print("Error: more inversions than alleles possible for ID " + str(invCounts[ID]))
				# self.printGenomes()
		return (mutationRemoved, inversionRemoved)

	# Simulating a single generational step
	def step(self):
		# print ("Generation " + str(self.age) + " Reproduction")
		# Pick mothers
		femaleSurvivals = []
		for female in self.females:
			# print(female.genome)
			femaleSurvivals += [female.survival()]
		# print("Female Survivals:" + str(femaleSurvivals))
		femaleSurvivals = [s/sum(femaleSurvivals) for s in femaleSurvivals]
		motherIndexes = np.random.choice(len(femaleSurvivals),self.size,p=femaleSurvivals)

		# Calculate necessary male values
		maleSurvivals = []
		# maleGenQuality = []
		for male in self.males:
			# print(male.genome)
			if self.noMaleCost:  #MAKE THIS LESS COSTLY TO COMPUTE (still runs .choice below)
				maleSurvivals += [1]
			else:
				maleSurvivals += [male.survival()]
			# maleGenQuality += [male.repQuality()]
		# print("Male Survivals:" + str(maleSurvivals))
		maleSurvivals = [s/sum(maleSurvivals) for s in maleSurvivals]

		# Populate the next generation by choosing fathers per mother and generating child genomes
		newMales = []
		newFemales = []
		# numMales = 0
		# numExpectedMales used for pegged equal-sex populations
		# NEED TO GUARANTEE ONE OF EACH SEX
		numExpectedMales = int(self.size/2)
		for motherIndex in motherIndexes:
			mother = self.females[motherIndex]

			encounteredMaleIndexes = np.random.choice(len(maleSurvivals),self.encounterNum,p=maleSurvivals)
			fatherIndex = encounteredMaleIndexes[0]
			# maxScore = self.males[fatherIndex].repQuality() + np.random.normal(scale=self.choiceNoiseSD)
			maxScore = 0 # (assumes rep quality score >= 0)
			for f in encounteredMaleIndexes:
				score = self.males[f].repQuality() + np.random.normal(scale=self.choiceNoiseSD)
				if score > maxScore:
					fatherIndex = f
					maxScore = score
			father = self.males[fatherIndex]

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

			# Check if using random sexes
			if self.randomSex:
				# Instantiate the new population member with a random sex
				if np.random.randint(2):
					newFemales += [individual('F',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
						self.minInvLen,self.lenChrom,self.isFly,self.willConvert,self.willRecombine,genome)]
				else:
					newMales += [individual('M',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
						self.minInvLen,self.lenChrom,self.isFly,self.willConvert,self.willRecombine,genome)]
			else:
				# Instantiate the new population member with determinate sex
				# MAY NEED TO ENSURE THAT MOTHERS AREN'T ALWAYS GIVING SONS, ETC
				if numExpectedMales > 0:
					newMales += [individual('M',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
						self.minInvLen,self.lenChrom,self.isFly,self.willConvert,self.willRecombine,genome)]
					numExpectedMales -= 1
				else:
					newFemales += [individual('F',self.mutEffectDiffSD,self.recombRate,self.conversionRate,\
						self.minInvLen,self.lenChrom,self.isFly,self.willConvert,self.willRecombine,genome)]


		self.females = newFemales
		self.males = newMales

		# Update the age of the population
		self.age = self.age + 1

		# # Remove fixed inversions and rearrange internal mutations
		# self.__removeFixedInv()
		# Remove fixed mutations and inversions, and rearrange mutations in removed inversions
		self.__removeFixed()

		if self.willMutate:
			# print ("Generation " + str(self.age) + " Mutation")
			# Sprinkle on mutations
			wholePop = self.males + self.females
			# Add new SA mutations
			numMuts = np.random.poisson(self.expectedNumMut)
			indivWithMut = np.random.randint(self.size, size=numMuts)
			for i in indivWithMut:
				# Add the mutation to the record, with start and end generations of polymorphism
				self.record[1] += [wholePop[i].mutate(self.__mutIDcount) + [self.age,-1]]
				self.record[2] += [[]]
				self.__mutIDcount += 1
				self.__mutFixed += [False]
				# self.record[2] += [[0 for t in range(len(self.record[0]))]]
		if self.willMutInv:
			wholePop = self.males + self.females
			# Add new inversions
			numInvMuts = np.random.poisson(self.expectedNumInvMut)
			indivWithInvMut = np.random.randint(self.size, size=numInvMuts)
			for i in indivWithInvMut:
				# Must account for failure due to no space on the genome
				invData = wholePop[i].mutateInv(self.__invIDcount)
				if not invData is None:
					self.record[3] += [invData + [self.age,-1]]
					self.record[4] += [[]]
					self.record[5] += [[]]
					self.record[6] += [[]]
					self.record[7] += [[]]
					self.__invIDcount += 1
					self.__invFixed += [False]
					# self.record[4] += [[0 for t in range(len(self.record[0]))]]
					# self.record[5] += [[0 for t in range(len(self.record[0]))]]
					# self.record[6] += [[0 for t in range(len(self.record[0]))]]
					# self.record[7] += [[0 for t in range(len(self.record[0]))]]
		# Could update phenotype values here
		return

	# Repeats the last element of the list, useful for fixed mutation and inversion accounting
	# # Has 
	# def __repeat(self,l):
	# 	l += l[len(l)-1]

	# For updating all recorded information on the population, record structured as:
	# [0] a list of ages at which an update was made
	# [1] a list of all mutations, indexed by ID, with each entry as 
	#      [position,survival effect,reproductive effect,chromosome,age at which mutation occured,age at fixation/loss]
	#      -1 means the mutation is not fixed/lost yet
	# [2] a list of counts, indexed by ID, of the mutation at each age noted in [0]
	# [3] a list of all inversions, indexed by ID, with each entry as 
	#      [start position,end position,chromosome,age at which mutation occured,age at fixation/loss]
	#      -1 means the mutation is not fixed/lost yet
	# [4] a list of counts, indexed by ID, of the inversion at each age noted in [0]
	# [5] a list of average count of mutations within the buffer region across the entire population,
	#      indexed by ID, of each inversion at each age noted in [0]
	# [6] a list of average survival effect within buffer region across the entire population,
	#      indexed by ID, of each inversion at each age noted in [0]
	# [7] a list of average reproductive effect within buffer region across the entire population,
	#      indexed by ID, of each inversion at each age noted in [0]
	# [8] a list of pop size at each recorded generation
	def __updateRecord(self):
		# print ("Record Update")
	 	# Record the age of the record update
		self.record[0] += [self.age]
	 	# Record the population size at the time of the of the record update
		self.record[8] += [self.size]
		# Count the number of each mutation/ID
		mutCounts = [0]*self.__mutIDcount
		invCounts = [0]*self.__invIDcount
		numMutInBuffer = [0]*self.__invIDcount
		# survEffectTotalInBuffer = [0]*self.__invIDcount
		survEffectTotalInBufferMultiplicative = [0]*self.__invIDcount
		reprEffectTotalInBuffer = [0]*self.__invIDcount
		for indiv in self.males + self.females:
			for chrom in indiv.genome:
				for hom in chrom:
					for mut in hom[0]:
						mutCounts[mut[3]] += 1
					lowerMutInd = 0
					upperMutInd = 0
					for inv in hom[1]:
						invCounts[inv[2]] += 1
						while lowerMutInd < len(hom[0]) and hom[0][lowerMutInd][0] < inv[0]-self.invRecBuffer:
							lowerMutInd += 1
						while upperMutInd < len(hom[0]) and hom[0][upperMutInd][0] < inv[1]+self.invRecBuffer:
							upperMutInd += 1
						mutInside = hom[0][lowerMutInd:upperMutInd]
						# Account for multiplicative survival effects - average effect is average of multiplied effects
						numMutInside = 0
						thisSurvEffect = 1
						thisReprEffect = 0
						for mut in mutInside:
							numMutInside += 1
							thisSurvEffect *= mut[1]
							thisReprEffect += mut[2]
						numMutInBuffer[inv[2]] += numMutInside
						survEffectTotalInBufferMultiplicative[inv[2]] += thisSurvEffect
						reprEffectTotalInBuffer[inv[2]] += thisReprEffect
		# Update the record, get averages for record[5,6,7]
		# for i in range(self.__mutIDcount):
		# 	count = mutCounts[i]
		# 	if count == 0:
		# 		if self.__mutFixed[i]:
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
		for i in range(self.__mutIDcount):
			count = mutCounts[i]
			if count != 0:
				self.record[2][i] += [count]
		for j in range(self.__invIDcount):
			count = invCounts[j]
			if count != 0:
				self.record[4][j] += [count]
				self.record[5][j] += [numMutInBuffer[j]/(float(count))]
				self.record[6][j] += [survEffectTotalInBufferMultiplicative[j]/(float(count))]
				self.record[7][j] += [reprEffectTotalInBuffer[j]/(float(count))]

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

	# Checks the record to see if all mutations/inversions are fixed
	def checkAllRecMutFixed(self):
		return all(self.__mutFixed) and all(self.__invFixed)


	# For simulating a number of generations sequentially
	def stepNGens(self, numGenerations):
		for i in range(numGenerations):
			self.step()
		return

	# For running a step and then updating the record.
	# CONSIDER JUST MAKING __updateRecord PUBLIC
	def recordStep(self):
		self.step()
		self.__updateRecord()
		return

	# For simulating a number of generations sequentially, recording after each
	def recordNGens(self, numGenerations):
		for i in range(numGenerations):
			self.step()
			self.__updateRecord()
		return

	# For simulating setSize*setNum generations sequentially, recording every setSize'th generation
	def recordNthForNSets(self, setSize, setNum):
		for i in range(setNum):
			self.stepNGens(setSize)
			self.__updateRecord()
		return

	# For simulating numGens generations sequentially, recording every setSize'th generation
	def recordNthInNGens(self, setSize, numGens):
		# Be wary of this, it may not simulate the full number of generations,
		#   but really you wouldn't necessarily record the final step anyway
		numSets = int(numGens/setSize)
		for i in range(numSets):
			self.stepNGens(setSize)
			self.__updateRecord()
		return

	# For simulating setSize*setNum generations sequentially, recording every setSize'th generation
	# Recording every generation after fixation as the fixation frequencies if fillRecord = True
	# Ideally only used with non-mutating populations
	def recordNthInNStopWhenFixed(self, setSize, numGens, fillRecord = True):
		# if willMutate or willMutInv:
		# 	raise SAIpop.InputError((willMutate,willMutInv),'For use with non-mutating populations')
		# Be wary of this, it may not simulate the full number of generations,
		#   but really you wouldn't necessarily record the final step anyway
		numSets = int(numGens/setSize)
		for i in range(numSets):
			self.stepNGens(setSize)
			self.__updateRecord()
			if self.checkAllRecMutFixed():
				if fillRecord:
					lastRecord = len(self.record[0])-1
					g = self.record[0][lastRecord]
					for addRecord in range(1,numSets-i):
						self.record[0] += [g+addRecord*setSize]
						self.record[8] += [self.size]
						# for m in range(self.__mutIDcount):
						# 	# print(self.record[2][m])
						# 	self.record[2][m] += [self.record[2][m][lastRecord]]
						# for i in range(self.__invIDcount):
						# 	self.record[4][i] += [self.record[4][i][lastRecord]]
						# 	self.record[5][i] += [self.record[5][i][lastRecord]]
						# 	self.record[6][i] += [self.record[6][i][lastRecord]]
						# 	self.record[7][i] += [self.record[7][i][lastRecord]]
				return
		return

	def printRecord(self):
		print (self.record)
		return

	def printGenomes(self):
		genomes = []
		for indiv in self.females:
			genomes += [indiv.genome]
		for indiv in self.males:
			genomes += [indiv.genome]
		print(genomes)

	# Takes a filename and mutation ID
	# Writes a tab delineated file of the generation and count data for that mutation
	def writeMutation(self,filename,ID,fillRecord = False):
		outfile = open(filename, 'w')
		# outfile.write('Mutation '+str(ID)+'\n')
		outfile.write('Generation\tCount\n')
		# outfile.write(str(self.record[0][0]) + '\t' + str(self.record[2][ID][0]) + '\n')
		# outfile.write(str(self.record[1][ID][4]) + '\t' + str(1) + '\n')
		# for t in range(1,len(self.record[0])):
		(firstGen, finalGen) = self.record[1][ID][4:6]
		if finalGen < firstGen:
			finalGen = self.age
		t = 0
		recordGen = self.record[0][t]
		while recordGen < firstGen:
			if fillRecord:
				outfile.write(str(recordGen) + '\t0\n')
			t += 1
			recordGen = self.record[0][t]
		offset = t
		while recordGen < finalGen:
			m = t - offset
			outfile.write(str(recordGen) + '\t' + str(self.record[2][ID][m]) + '\n')
			t += 1
			recordGen = self.record[0][t]
		if fillRecord:
			while t < len(self.record[0]):
				if self.__mutFixed[ID]:
					outfile.write(str(recordGen) + '\t' + str(self.record[8][t]) + '\n')
				else:
					outfile.write(str(recordGen) + '\t0\n')
				t += 1
				recordGen = self.record[0][t]
		# for t in range(len(self.record[0])):
		# 	outfile.write(str(self.record[0][t]) + '\t' + str(self.record[2][ID][t]) + '\n')
		outfile.close()

	# Takes a filename and inversion ID
	# Writes a tab delineated file of the generation, count, 
	#  average number of mutations in buffer, average cumulative survival,
	#  and average cumulative effect data for that inversion
	def writeInversion(self,filename,ID,fillRecord = False):
		outfile = open(filename, 'w')
		# outfile.write('Inversion '+str(ID)+'\n')
		outfile.write('Generation\tCount\tAvgNumMut\tAvgSurEff\tAvgRepEff\n')
		# outfile.write(str(self.record[0][0]) + '\t' + str(self.record[4][ID][0]) + '\t'\
		# 	+ str(self.record[5][ID][0]) + '\t' + str(self.record[6][ID][0]) + '\t'\
		# 	+ str(self.record[7][ID][0]) + '\n')
		# # May want to find a way to acount for the number of mutations in the inversion
		# #   and their effects upon generation
		# outfile.write(str(self.record[3][ID][3]) + '\t' + str(1) + '\t' + 'NA' + '\t'\
		# 	+ 'NA' + '\t' + 'NA' + '\n')
		# for t in range(1,len(self.record[0])):
		(firstGen, finalGen) = self.record[3][ID][3:5]
		if finalGen < firstGen:
			finalGen = self.age
		t = 0
		recordGen = self.record[0][t]
		while recordGen < firstGen:
			if fillRecord:
				outfile.write(str(recordGen) + '\t0\n')
			t += 1
			recordGen = self.record[0][t]
		offset = t
		while recordGen < finalGen:
			i = t - offset
			outfile.write(str(recordGen) + '\t' + str(self.record[4][ID][i]) + '\t'\
				+ str(self.record[5][ID][i]) + '\t' + str(self.record[6][ID][i]) + '\t'\
				+ str(self.record[7][ID][i]) + '\n')
			t += 1
			recordGen = self.record[0][t]
		if fillRecord:
			while t < len(self.record[0]):
				if self.__invFixed[ID]:
					outfile.write(str(recordGen) + '\t' + str(self.record[8][t]) + '\t-1\t-1\t-1\n')
				else:
					outfile.write(str(recordGen) + '\t0\t-1\t-1\t-1\n')
				t += 1
				recordGen = self.record[0][t]
		# for t in range(len(self.record[0])):
		# 	outfile.write(str(self.record[0][t]) + '\t' + str(self.record[4][ID][t]) + '\t'\
		# 		+ str(self.record[5][ID][t]) + '\t' + str(self.record[6][ID][t]) + '\t'\
		# 		+ str(self.record[7][ID][t]) + '\n')
		outfile.close()

	# Takes a filename and writes a tab delineated file of the mutation frequency 
	#   by generation and ID with generation in the first column
	def writeAllMutFreqTable(self,filename):
		outfile = open(filename, 'w')
		header = 'Generation'
		mutInit = []
		mutFinal = []
		offsets = []
		for m in range(self.__mutIDcount):
			# header += '\t'+str(m)
			header += '\tM'+str(m)
			initGen = self.record[1][m][4]
			finalGen = self.record[1][m][5]
			mutInit += [initGen]
			if initGen > finalGen:
				mutFinal += [self.age]
			else:
				mutFinal += [finalGen]
			offsets += [-1]
		# print(mutInit)
		# print(mutFinal)
		# print(self.record[1])
		outfile.write(header + '\n')
		for g in range(len(self.record[0])):
			gen = self.record[0][g]
			genLine = str(gen) 
			for m in range(self.__mutIDcount):
				if gen < mutInit[m]:
					genLine += '\t0'
				elif gen < mutFinal[m]: # The last generation isn't recorded as it is removed first
					if offsets[m] >= 0:
						# print("g "+str(g)+" gen "+str(gen)+" Init "+str(mutInit[m])+" Final "+str(mutFinal[m]))
						# print("Offset "+str(offsets[m]))
						# print("Index "+str(g-offsets[m]))
						# print(self.record[2][m][g-offsets[m]-2:])
						genLine += '\t' + str(self.record[2][m][g-offsets[m]])
					else:
						# print("Setting Offset:")
						# print(" g|offset "+str(g)+" gen "+str(gen)+" Init "+str(mutInit[m])+" Final "+str(mutFinal[m]))
						offsets[m] = g
						genLine += '\t' + str(self.record[2][m][0])
				elif self.__mutFixed[m]:
					genLine += '\t' + str(self.record[8][g])
				else:
					genLine += '\t0'
			outfile.write(genLine + '\n')
		outfile.close()

	# Takes a filename and writes a tab delineated file of the inversion frequency
	#   by generation and ID with generation in the first column
	def writeAllInvFreqTable(self,filename):
		outfile = open(filename, 'w')
		header = 'Generation'
		invInit = []
		invFinal = []
		offsets = []
		for i in range(self.__invIDcount):
			# header += '\t'+str(i)
			header += '\tI'+str(i)
			initGen = self.record[3][i][3]
			finalGen = self.record[3][i][4]
			invInit += [initGen]
			if initGen > finalGen:
				invFinal += [self.age]
			else:
				invFinal += [finalGen]
			offsets += [-1]
		# print(invInit)
		# print(invFinal)
		# print(self.record[3])
		outfile.write(header + '\n')
		for g in range(len(self.record[0])):
			gen = self.record[0][g]
			genLine = str(gen) 
			for i in range(self.__invIDcount):
				if gen < invInit[i]:
					genLine += '\t0'
				elif gen < invFinal[i]:
					if offsets[i] >= 0:
						# print("g "+str(g)+" gen "+str(gen)+" Offset "+str(offsets[i])+" Index "+str(g-offsets[i])+" Init "+str(invInit[i])+" Final "+str(invFinal[i]))
						# print(len(self.record[4][i]))
						# print(self.record[4][i][g-offsets[i]-2:])
						genLine += '\t' + str(self.record[4][i][g-offsets[i]])
					else:
						# print("Setting Offset:")
						# print("g|offset "+str(g)+" gen "+str(gen)+" Init "+str(invInit[i])+" Final "+str(invFinal[i]))
						# print(self.record[4][i])
						offsets[i] = g
						genLine += '\t' + str(self.record[4][i][0])
				elif self.__invFixed[i]:
					genLine += '\t' + str(self.record[8][g])
				else:
					genLine += '\t0'
			outfile.write(genLine + '\n')
		outfile.close()

	# Takes a filename and writes a tab delineated file of the position, effect, chromosome,
	#   and initial generation data for all mutations
	def writeMutCharTable(self,filename):
		outfile = open(filename, 'w')
		# outfile.write('Position\tSurEffect\tRepEffect\tChromosome\tInitGen\tFinalGen\n')
		outfile.write('Pos\tSurEf\tRepEf\tChrom\tInitGen\tFinalGen\n')
		for m in range(len(self.record[1])):
			line = ''
			for datum in self.record[1][m]:
				line += str(datum) + '\t'
			outfile.write(line[:-1] + '\n')
		outfile.close()

	# Takes a filename and writes a tab delineated file of the positions, chromosome,
	#   and initial generation data for all inversions
	def writeInvCharTable(self,filename):
		outfile = open(filename, 'w')
		# outfile.write('Position1\tPosition2\tChromosome\tInitGen\tFinalGen\n')
		outfile.write('Pos1\tPos2\tChrom\tInitGen\tFinalGen\n')
		for i in range(len(self.record[3])):
			line = ''
			for datum in self.record[3][i]:
				line += str(datum) + '\t'
			outfile.write(line[:-1] + '\n')
		outfile.close()

	# Writes a summary file for the parameters of the simulation run
	# UPDATE - changeable size and additional parameters
	def writeSummary(self,filename):
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
		outfile.close()
		return

	# Writes record to outfiles as such:
	# One summary document
	# Tab delineated table of mutation characteristics for all mutations
	# Tab delineated table of inversion characteristics for all inversions
	# Tab delineated table with generation and count data for every mutation
	# Tab delineated table with generation, count, etc. additional data for every inversion
	def writeRecordTables(self,outFilePrefix):
		self.writeSummary(outFilePrefix+'ParamSumm.txt')
		self.writeMutCharTable(outFilePrefix+'MutSumm.txt')
		self.writeInvCharTable(outFilePrefix+'InvSumm.txt')
		self.writeAllMutFreqTable(outFilePrefix+'MutFreqs.txt')
		self.writeAllInvFreqTable(outFilePrefix+'InvFreqs.txt')
		for mutID in range(self.__mutIDcount):
			filename = outFilePrefix+'Mut'+str(mutID)+'.txt'
			self.writeMutation(filename,mutID)
		for invID in range(self.__invIDcount):
			filename = outFilePrefix+'Inv'+str(invID)+'.txt'
			self.writeInversion(filename,invID)
		return

	# # Scratch for testing mutation data change (due to shared reference)
	# def testSharedData(self):
	# 	mutCounts = [0]*self.__mutIDcount
	# 	for indiv in self.males + self.females:
	# 		for chrom in indiv.genome:
	# 			for hom in chrom:
	# 				for mut in hom[0]:
	# 					mutCounts[mut[3]] += 1
	# 	for i in range(self.__mutIDcount):
	# 		if mutCounts[i] == 2*self.size:
	# 			# self.record[1][i] = mut
	# 			# mut[2] = 'test'
	# 			print('Starting test for ID '+ str(i))
	# 			# print(mut)
	# 			# print(self.record[1][i])
	# 			firstEncounter = True
	# 			for indiv in self.males + self.females:
	# 				for chrom in indiv.genome:
	# 					for hom in chrom:
	# 						for mut in hom[0]:
	# 							if mut[3] == i:
	# 								if firstEncounter:
	# 									mut[2] = 'test'
	# 									firstEncounter = False
	# 								print(mut)
	# 	return




