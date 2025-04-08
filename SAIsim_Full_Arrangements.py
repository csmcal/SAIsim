
# Script created for running with Condor, to take in mutation rates and run the full SAIsim model

import sys
import SAIsim as sim
# import numpy as np
from os.path import exists
import time
import datetime

# print(sys.argv)

popPickleFile = str(sys.argv[1])
timeToCheckpoint = int(sys.argv[2])
mutRate = float(sys.argv[3])
mutRateInv = float(sys.argv[4])
encounterNum = int(sys.argv[5])
recombRate = float(sys.argv[6])
convRate = float(sys.argv[7])
popSize = int(sys.argv[8])
numGens = int(sys.argv[9])
wrightFisherSurvivals = bool(int(sys.argv[10]))
noMaleCost = bool(int(sys.argv[11]))
noFemaleCost = bool(int(sys.argv[12]))
replicateNum = int(sys.argv[13])


# Other set parameters as described in SAIsim
choiceNoiseSD = 1
invRecBuffer = 1
mutEffectDiffSD = .2
# minInvLen = .1
minInvLen = .05

recordEveryN = max(10,numGens//100)
recordTwice = True

if wrightFisherSurvivals:
	sampleSurvPerOff = True
	# sampleEncWReplacement = True
else:
	sampleSurvPerOff = False
	# sampleEncWReplacement = False
sampleEncWReplacement = True

# Set the script start time
startTime = time.time()
print("Current time:\t"+datetime.datetime.fromtimestamp(startTime).strftime('%Y-%m-%d %H:%M:%S'))


# Check for a checkpoint and restart from it if appropriate
if exists(popPickleFile):
	pop = sim.readPopPickle(popPickleFile)
	print("Restarting Simulation from "+popPickleFile)
else:
	# Generate a new population
	pop = sim.SAIpop(popSize, mutRate, mutRateInv, mutEffectDiffSD, minInvLen, convRate,
		recombRate, encounterNum, choiceNoiseSD, invRecBuffer,
		sampleSurvivorsPerOffspring = sampleSurvPerOff,
		sampleEncountersWithReplacement = sampleEncWReplacement, 
		noFemaleCost = noFemaleCost,
		noMaleCost = noMaleCost,
		oneInvPerChrom = False,
		removeFixed = True, 
		testing = False,
		verbose = False,
		storePickledPop = False,
		minNumGensKeepRecord = 400,
		minFreqKeepInvRecord = 0)
	print("Simulation Start")


# Start the Simulation
currGen = pop.age
# pop.printParamSummary()
numSetsTotal = numGens//recordEveryN
numSetsCurr = currGen//recordEveryN
if currGen % recordEveryN != 0:
	print("WARNING - population was restarted on a generation ("+str(currGen)+") other than the expected recording generation, "+\
		"or the end of the recording/checkpoint period (some multiple of "+str(recordEveryN)+", < "+str(numGens)+")")
currSet = numSetsCurr

# Run until checkpoint time
while currSet < numSetsTotal:
	print("Generation "+str(currSet*recordEveryN))
	if currSet == numSetsTotal-1 or currSet == numSetsTotal//2-1:
		pop.recordArrNthGen(recordEveryN,doubleRecord=recordTwice)
	else:
		pop.stepNGens(recordEveryN)
	currTime = time.time()
	if currTime-startTime > timeToCheckpoint:
		currTimeStr = datetime.datetime.fromtimestamp(currTime).strftime('%Y-%m-%d %H:%M:%S')
		elapsedTimeStr = str(datetime.timedelta(seconds=currTime-startTime))
		print("Checkpointing\nElapsed time:\t"+elapsedTimeStr+"\nCurrent time:\t"+currTimeStr)
		pop.writePopPickle(popPickleFile)
		print("Checkpointed Simulation to "+popPickleFile)
		sys.exit(85)
	currSet += 1

# Record Relevant Output
print("Generation "+str(pop.age))
print("Simulation Finished")
currTime = time.time()
currTimeStr = datetime.datetime.fromtimestamp(currTime).strftime('%Y-%m-%d %H:%M:%S')
elapsedTimeStr = str(datetime.timedelta(seconds=currTime-startTime))
print("Writing Arrangement Records\nElapsed time:\t"+elapsedTimeStr+"\nCurrent time:\t"+currTimeStr)
filePrefix = ''
outfileArr = filePrefix+'ArrangementData.txt'
pop.writeArrangements(outfileArr)
outfileInvs = filePrefix+'InvSumm.txt'
pop.writeInvCharTable(outfileInvs)
print("Arrangement Record Written")


