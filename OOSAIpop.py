
# The Object Oriented implementation of Sexually Antagonistic mutation accumulation
#     in Inversions in a forward Population simulation
# Ideally, models large inversion polymorphisms in populations with high reproductive skew

# REWRITE WITH NUMPY ARRAYS/VECTORS

import numpy as np

# The object representing an individual in the population under simulation,
# contains methods for generating mutation position and effects, and new recombinant gametes
# recombRate is the expected number of recombination events per chromosome per reproductive event
# meanMutEffect is the expected effect size of a new mutation
# mutEffectDiffSD is the SD of the difference in survival/rep effect size of a new mutation (normal dist)
# 2 copies of a single chormosome,
# CONSIDER INSTEAD: genome = [[[chr1hom1],[chr1hom2]],[[chr2hom1],[chr2hom2]],..] (implement if time)
class individual(object):
	"""individuals represent members of simSAIpopulations"""
	def __init__(self, sex, mutEffectDiffSD, recombRate, minInvLen, genome = [[[[],[]],[[],[]]]]):
		# super(individual, self).__init__()
		self.sex = sex
		self.mutEffectDiffSD = mutEffectDiffSD
		self.recombRate = recombRate
		self.minInvLen = minInvLen
		# genome = [[[chr1hom1],[chr1hom2]],[[chr2hom1],[chr2hom2]],..]
		# where chomosome homologs = [mutist, InvList]
		self.genome = genome


	# For getting the insertion index of a position from a list of [position,..] lists
	def __getInsInd(self,mutInvList,mutInvPos):
		i = 0 
		while (i < len(mutInvList)) and (mutInvList[i][0] <= mutInvPos):
			i += 1
		return i

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

	#Generates a mutation of the form [position, survival effect, rep effect, ID]
	def mutate(self,ID):
		# print self.genome
		mutPos = np.random.ranf()
		mutEffects = self.__genEffectSizes()
		mutation = [mutPos]+mutEffects+[ID]
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		print self.genome[chromIndex][homIndex][0]
		self.genome[chromIndex][homIndex][0] = self.__insert(self.genome[chromIndex][homIndex][0],mutation)
		print self.genome[chromIndex][homIndex][0]
		# Record the mutation data
		self.record[1] += [mutation[0:3]+[chromIndex]]
		return

	# Generates and inserts an inversion into the genome, umless there is no open region >= minInvLen
	# Removes resampling 
	def mutateInv(self,ID):
		# Pick a chromosome
		chromIndex = np.random.randint(0,len(self.genome))
		# Put it on one of the two homologs
		homIndex = np.random.randint(0,2)
		chromHomInv = self.genome[chromIndex][homIndex][1]
		# Pick where the in
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
		length = 1 - potentialStart
		if length > self.minInvLen:
			openRegLengths += [length]
			openRegStarts += [potentialStart]
			openRegIndexes += [index]
		# If there is no space for a new inversion >= minInvLen, don't add one
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
		# Record the inversion data
		self.record[3] += [inversion[0:2]+[chromIndex]]
		return

	# For updating the repQuality and survival if pre-calculated
	def updatePhenotypes(self):
		return

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
			end1 = invHom1[ind1][1]
			end2 = invHom2[ind2][1]
			# WHAT TO DO WITH CROSSOVERS == ONE END OF AN INVERSION
			# Currently, inversions are counted as [1,2)  (?)
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
		for chrom in self.genome:
			numRecomb = np.random.poisson(self.recombRate)
			currHom1 = bool(np.random.randint(2))
			if numRecomb == 0:
				if currHom1:
					gamGenome += [chrom[0][:]]
				else:
					gamGenome += [chrom[1][:]]
			else:
				recombPositions = np.random.ranf(numRecomb)
				
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




		
# The object representing the population under simulation, contains methods for stepping through generations
# mutRate is the expected number of new mutations per chromosome per reproductive event
# mutRateInv is the expected number of new inversions per chromosome per reproductive event
# recombRate is the expected number of recombination events per chromosome per reproductive event
# meanMutEffect is the expected effect size of a new mutation
# mutEffectDiffSD is the SD of the difference in survival/rep effect size of a new mutation (normal dist)
# record is a general variable for storing generation statistics about the population
class simSAIpopulation(object):
	"""simSAIpopulation represents sexually reproducing populations for sexually antagonistic inversion simulation"""
	def __init__(self, size, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, recombRate, encounterNum, choiceNoiseSD):
		# super(simSAIpopulation, self).__init__()
		self.size = size
		self.mutRate = mutRate
		self.mutRateInv = mutRateInv
		self.expectedNumMut = mutRate * size
		self.expectedNumInvMut = mutRateInv * size
		# Keep trck of next mutation (or inversion?) ID
		self.__mutIDcount = 0
		self.__invIDcount = 0
		self.mutEffectDiffSD = mutEffectDiffSD
		self.minInvLen = minInvLen
		self.recombRate = recombRate
		self.encounterNum = encounterNum
		self.choiceNoiseSD = choiceNoiseSD
		self.males = []
		self.females = []
		# Record keeping variables
		self.record = [[],[],[],[]]
		self.age = 0
		# Generate a set of 'size' individuals with random sex and no mutations or inversions
		numMales = np.random.binomial(self.size,0.5)
		self.males = [individual('M',mutEffectDiffSD,recombRate,minInvLen) for i in range(numMales)]
		self.females = [individual('F',mutEffectDiffSD,recombRate,minInvLen) for j in range(self.size-numMales)]


	# # Generates a new individual from a reproductive event
	# def __genChild(mother,father):
	# 	mother = self.females
	# 	return

	# Simulating a single generational step
	def step(self):
		print "Generation " + str(self.age)
		# Pick mothers
		femaleSurvivals = []
		for female in self.females:
			# print female.genome
			femaleSurvivals += [female.survival()]
		print "Female Survivals:" + str(femaleSurvivals)
		femaleSurvivals = [s/sum(femaleSurvivals) for s in femaleSurvivals]
		motherIndexes = np.random.choice(len(femaleSurvivals),self.size,p=femaleSurvivals)

		# Calculate necessary male values
		maleSurvivals = []
		# maleGenQuality = []
		for male in self.males:
			# print male.genome
			maleSurvivals += [male.survival()]
			# maleGenQuality += [male.repQuality()]
		print "Male Survivals:" + str(maleSurvivals)
		maleSurvivals = [s/sum(maleSurvivals) for s in maleSurvivals]

		# Populate the next generation by choosing fathers per mother and generating child genomes
		newMales = []
		newFemales = []
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
			# print "father gamete" + str(fatherGamete)
			motherGamete = mother.genGamete()
			# print "mother gamete" + str(motherGamete)
			for chrom in range(len(fatherGamete)):
				genome += [[fatherGamete[chrom],motherGamete[chrom]]]
			# print genome

			# Instantiate the new population member with a random sex
			if np.random.randint(0,2):
				newFemales += [individual('F',self.mutEffectDiffSD,self.recombRate,self.minInvLen,genome)]
			else:
				newMales += [individual('M',self.mutEffectDiffSD,self.recombRate,self.minInvLen,genome)]

		self.females = newFemales
		self.males = newMales

		# Sprinkle on mutations
		wholePop = self.males + self.females
		# Add new SA mutations
		numMuts = np.random.poisson(self.expectedNumMut)
		indivWithMut = np.random.randint(self.size, size=numMuts)
		for i in indivWithMut:
			wholePop[i].mutate(self.__mutIDcount)
			self.__mutIDcount += 1
		# Add new inversions
		numInvMuts = np.random.poisson(self.expectedNumInvMut)
		indivWithInvMut = np.random.randint(self.size, size=numMuts)
		for i in indivWithInvMut:
			wholePop[i].mutateInv(self.__invIDcount)
			self.__invIDcount += 1

		self.age = self.age + 1
		# Could update phenotype values here
		return

	def __updateRecord():
		#Count the number of each mutation/ID
		mutCounts = [0]*self.age
		invCounts = [0]*self.age
		for indiv in self.males + self.females:
			for chrom in indiv.genome:
				for hom in chrom:
					for mut in hom[0]:
						mutCounts[mut[3]] += 1
					for inv in hom[1]:
						invCounts[inv[2]] += 1
		for i in mutCounts:
			currentRec = 
			self.record[0][i] += [self.age,]
		for j in self.__invIDcount:
			self.record[1][i] += [self.age,]

	# For getting mutation/inversion statistics from a step of the simulation
	# returns the 
	def getStdStats(self):
		stats = []
		return stats

	# For simulating a number of generations sequentially
	def stepNGens(self, numGenerations):
		for i in range(numGenerations):
			self.step()
		return

	def recordStep(self):


	def recordNGens(self, numGenerations):
		for i in range(numGenerations):
			self.step()
			self.record += [self.getStdStats]
		return

	# Writes record to an outfile
	def writeData(self,outFileName):
		return

	def testStop(self):
		return




