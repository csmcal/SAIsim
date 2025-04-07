# PLotting Per-capita or relative Paternity in SAIsim results  (with and w/o Inversions)
# For generating plots of two-mutation simulations
library(tidyverse)
library(parallel)


# Example run: Rscript PaternityPlots.R ./dmelN3/

# Collect input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  surEffects <- seq(0.0,1.0-0.025,by=0.025)
  repEffects <- seq(0.025,1.0,by=0.025)
  numReps <- 500
  popSize <- 1000
  resultsDir = './'
  print("Defaulting to using:")
  print(paste0("resultsDir = ",resultsDir))
  print(paste0("numReps = ",numReps))
  print(paste0("popSize = ",popSize))
  print("surEffects = ")
  print(surEffects)
  print("repEffects = ")
  print(repEffects)
  
  prefixList <- list("EGSE2","EGSE3","EGSE4","EGSE10","EGSE100",
                     "EGNME2","EGNME3","EGNME4","EGNME10","EGNME100",
                     "EGNFE2","EGNFE3","EGNFE4","EGNFE10","EGNFE100")
  dataDirList <- list()
  for(pre in prefixList){
    dataDirList <- append(dataDirList,paste0(pre,"/"))
  }
  print("Running on directories:")
  print(dataDirList)
} else if (length(args)==2) {
  prefix = args[1]
  dataDir = args[2]
  print(paste0("Running on dataDir = ",dataDir))
  print(paste0("Using prefix label = ",prefix))
  prefixList <- list(prefix)
  dataDirList <- list(dataDir)
  
  surEffects <- seq(0.0,1.0-0.025,by=0.025)
  repEffects <- seq(0.025,1.0,by=0.025)
  numReps <- 500
  popSize <- 1000
  resultsDir = './'
  print("For sim parameters, defaulting to using:")
  print(paste0("resultsDir = ",resultsDir))
  print(paste0("numReps = ",numReps))
  print(paste0("popSize = ",popSize))
  print("surEffects = ")
  print(surEffects)
  print("repEffects = ")
  print(repEffects)
} else if (length(args)==11) {
  prefix = args[1]
  dataDir = args[2]
  print(paste0("Running on dataDir = ",dataDir))
  print(paste0("Using prefix label = ",prefix))
  prefixList <- list(prefix)
  dataDirList <- list(dataDir)
  
  surMin <- as.numeric(args[3])
  surMax <- as.numeric(args[4])
  surStep <- as.numeric(args[5])
  surEffects <- seq(surMin,surMax,by=surStep)
  
  repMin <- as.numeric(args[6])
  repMax <- as.numeric(args[7])
  repStep <- as.numeric(args[8])
  repEffects <- seq(repMin,repMax,repStep)
  
  numReps <- as.integer(args[9])
  popSize <- as.integer(args[10])
  resultsDir <- args[11]
  
  print("For sim parameters, defaulting to using:")
  print(paste0("resultsDir = ",resultsDir))
  print(paste0("numReps = ",numReps))
  print(paste0("popSize = ",popSize))
  print("surEffects = ")
  print(surEffects)
  print("repEffects = ")
  print(repEffects)
} else {
  stop(paste0("Argument interpretation not implemented yet for ",length(args)," args"), call.=FALSE)
}

if (substr(resultsDir,nchar(resultsDir),nchar(resultsDir))!="/") {
  resultsDir <- paste0(resultsDir,'/')
}




getFinalHapData <- function(spacing,hapFileName) {
  hapFreqTable <- read.table(hapFileName,sep='\t',header=TRUE,na.strings = c("None","NA"))
  finalGen <- tail(hapFreqTable$Generation,1)
  finalGenHapData <- filter(hapFreqTable,Generation==finalGen)
  finalAdultGen <- finalGen - 1
  finalAdultGenHapData <- filter(hapFreqTable,Generation==finalAdultGen,LifeHistoryStage=="Adult")
  hapData <- rbind(finalAdultGenHapData,finalGenHapData)
  return(hapData)
}

getHapDataAtPar <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  s <- format(spacing,nsmall = 4)
  f <- format(initFreq,nsmall = 3)
  i <- as.integer(hasInv)
  fileBase <- paste(filePrefix,'pos',s,'inv',i,'freq',f,'n', sep="")
  finalHapData <- tibble()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      # print(paste0("Trying ",fileBase,n))
      hapFileName <- paste(fileBase,n,"Hap.txt",sep="")
      finalHapData <- rbind(finalHapData,getFinalHapData(spacing,hapFileName))
    })
  }
  return(finalHapData)
}

genHapTable <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  hapTable <- tibble()
  for (s in spacing) {
    hapTable <- rbind(hapTable,getHapDataAtPar(filePrefix,s,numReps,initFreq,hasInv))
  }
  # hapTable <- as.tibble(hapTable)
  return(hapTable)
}


combineHaps <- function(hapTable,hapNames) {
  comb <- tibble()
  for (name in hapNames) {
   comb <- rbind(comb,filter(hapTable,Haplotype==name)) 
  }
  # comb <- mutate_if(select(comb,!c(Generation,LifeHistoryStage,Haplotype)),is.character,as.numeric)
  comb <- select(comb,!c(Generation,LifeHistoryStage,Haplotype))
  # print(comb)
  comb <- summarise_all(comb, sum)
  # print(comb)
  # Generate the per capita parentage/descendant allele counts
  DescendantsPerAllele <- comb$OffspringCount/comb$Count
  ParentagePerHomozygote <- comb$OffspringCountInHomozygote/comb$Homozygote
  ParentagePerHeterozygote <- comb$OffspringCountInHeterozygote/comb$Heterozygote
  DescendantsPerFemaleAllele <- comb$OffspringCountInFemale/comb$FemaleCount
  MaternityPerFemaleHomozygote <- comb$OffspringCountInFemaleHomozygote/comb$FemaleHomozygote
  MaternityPerFemaleHeteroygote <- comb$OffspringCountInFemaleHeterozygote/comb$FemaleHeterozygote
  DescendantsPerMaleAllele <- comb$OffspringCountInMale/comb$MaleCount
  PaternityPerMaleHomozygote <- comb$OffspringCountInMaleHomozygote/comb$MaleHomozygote
  PaternityPerMaleHeterozygote <- comb$OffspringCountInMaleHeterozygote/comb$MaleHeterozygote
  # Attatch to the table to return
  combPerCap <- tibble(DescendantsPerAllele,ParentagePerHomozygote,ParentagePerHeterozygote,
                       DescendantsPerFemaleAllele,MaternityPerFemaleHomozygote,MaternityPerFemaleHeteroygote,
                       DescendantsPerMaleAllele,PaternityPerMaleHomozygote,PaternityPerMaleHeterozygote)
  comb <- cbind(comb,combPerCap)
  # print(comb)
  return(comb)
}

calcPerCapHapParentage <- function(hapTable) {
  adultTable <- filter(hapTable,LifeHistoryStage=="Adult")
  AB <- cbind(Haplotype=c("AB"),combineHaps(adultTable,c("0:1;","0:1;0")))
  Ab <- cbind(Haplotype=c("Ab"),combineHaps(adultTable,c("1;","1;0")))
  aB <- cbind(Haplotype=c("aB"),combineHaps(adultTable,c("0;","0;0")))
  ab <- cbind(Haplotype=c("ab"),combineHaps(adultTable,c(";",";0")))
  combPerCapHapParentage <- rbind(AB,Ab,aB,ab)
  return(combPerCapHapParentage)
}


combineHapsByLHS <- function(hapTable,hapNames) {
  selHaps <- tibble()
  for (name in hapNames) {
    selHaps <- rbind(selHaps,filter(hapTable,Haplotype==name)) 
  }
  adult <- select(filter(selHaps,LifeHistoryStage=="Adult"),!c(Generation,LifeHistoryStage,Haplotype))
  adultSums <- summarise_all(adult, sum)
  zygote <- select(filter(selHaps,LifeHistoryStage=="Zygote"),!c(Generation,LifeHistoryStage,Haplotype))
  zygoteSums <- summarise_all(zygote, sum)
  # Generate the per capita parentage/descendant allele counts
  DescendantsPerAllele <- adultSums$OffspringCount/zygoteSums$Count
  ParentagePerHomozygote <- adultSums$OffspringCountInHomozygote/zygoteSums$Homozygote
  ParentagePerHeterozygote <- adultSums$OffspringCountInHeterozygote/zygoteSums$Heterozygote
  DescendantsPerFemaleAllele <- adultSums$OffspringCountInFemale/zygoteSums$FemaleCount
  MaternityPerFemaleHomozygote <- adultSums$OffspringCountInFemaleHomozygote/zygoteSums$FemaleHomozygote
  MaternityPerFemaleHeteroygote <- adultSums$OffspringCountInFemaleHeterozygote/zygoteSums$FemaleHeterozygote
  DescendantsPerMaleAllele <- adultSums$OffspringCountInMale/zygoteSums$MaleCount
  PaternityPerMaleHomozygote <- adultSums$OffspringCountInMaleHomozygote/zygoteSums$MaleHomozygote
  PaternityPerMaleHeterozygote <- adultSums$OffspringCountInMaleHeterozygote/zygoteSums$MaleHeterozygote
  # Attatch to the table to return
  combPerCap <- tibble(DescendantsPerAllele,ParentagePerHomozygote,ParentagePerHeterozygote,
                       DescendantsPerFemaleAllele,MaternityPerFemaleHomozygote,MaternityPerFemaleHeteroygote,
                       DescendantsPerMaleAllele,PaternityPerMaleHomozygote,PaternityPerMaleHeterozygote)
  # print(combPerCap)
  return(combPerCap)
}

calcDescendantsPerAllelle <- function(hapTable) {
  AB <- cbind(Haplotype=c("AB"),combineHapsByLHS(hapTable,c("0:1;","0:1;0")))
  Ab <- cbind(Haplotype=c("Ab"),combineHapsByLHS(hapTable,c("1;","1;0")))
  aB <- cbind(Haplotype=c("aB"),combineHapsByLHS(hapTable,c("0;","0;0")))
  ab <- cbind(Haplotype=c("ab"),combineHapsByLHS(hapTable,c(";",";0")))
  combPerCapHapParentage <- rbind(AB,Ab,aB,ab)
  return(combPerCapHapParentage)
}


# To plot the reproductive values calculated from the two-mutant scenarios
plotAlleleDescPerHap <- function(combDescPerAllele) {
  # Reshape the haplotype data
  plotData <- select(combDescPerAllele,c(Haplotype,
                                         DescendantsPerAllele,
                                         DescendantsPerFemaleAllele,
                                         DescendantsPerMaleAllele))
  names(plotData) <- c("Haplotype","Both","Female","Male")
  plotData <- gather(plotData,
                     key = "Sex",
                     value = "DescendantsPerAllele",
                     -Haplotype)
  plotData$Haplotype <- factor(plotData$Haplotype,
                               levels = c("ab", "aB", "Ab", "AB" ))
  plotData$Sex <- factor(plotData$Sex,
                               levels = c("Male","Female","Both"))
  # Plot
  library(ggplot2)
  p <-ggplot(plotData, aes(Haplotype, DescendantsPerAllele, width=.5))
  p <- p +geom_col(aes(fill = Sex),position = "dodge")
  p <- p +geom_hline(aes(yintercept=1.0), linetype=2)
  p <- p +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p +scale_fill_manual(values = c(rep(c("gray55", "gray75", "lightblue"))))
  # p <- p +scale_fill_manual(values = c(rep(c('#228833', '#4477AA', '#CCBB44'))))
  p <- p +scale_y_continuous(breaks=c(0.75,1.0,1.25))
  p <- p +ylab("Average Descendants per Allele")
  p <- p +coord_cartesian(ylim=c(.65,1.3))
  p <- p +ggtitle("Haplotypes Exhibit Synergistic Epistasis")
  ggsave("Haplotypes Exhibit Synergistic Epistasis.png", 
         units="in", width=5, height=3, dpi=600, device = 'png')
}


reshapeBySex <- function(hapTable) {
  # print(hapTable)
  # newNames <- c("Sex",rep(c("Allele","Homozygote","Heterozygote"),3))
  newNames <- c("Sex","Allele","Homozygote","Heterozygote")
  both <- mutate(select(hapTable,!c(4:9)),Sex="Both",.before=1)
  names(both) <- newNames
  female <- mutate(select(hapTable,!c(1:3,7:9)),Sex="Female",.before=1)
  names(female) <- newNames
  male <- mutate(select(hapTable,!c(1:6)),Sex="Male",.before=1)
  names(male) <- newNames
  reshaped <- rbind(both,female,male)
  # print(reshaped)
  reshaped <- gather(reshaped,
                     key = "Zygosity",
                     value = "Count",
                     -Sex)
  # print(reshaped)
  return(reshaped)
}

combineReshapeHaps <- function(hapTable,hapNames) {
  selHap <- tibble()
  for (name in hapNames) {
    selHap <- rbind(selHap,filter(hapTable,Haplotype==name)) 
  }
  zygote <- select(filter(selHap,LifeHistoryStage=="Zygote"),
                   !c(Generation,Haplotype,
                      OffspringCount:OffspringCountInMaleHeterozygote))
  # print(zygote)
  zygoteSums <- summarise_at(zygote, -1, sum)
  # print(zygoteSums)
  zygoteSexed <- mutate(reshapeBySex(zygoteSums),LifeHistoryStage="Zygote",.before=1)
  
  adult <- select(filter(selHap,LifeHistoryStage=="Adult"),
                  !c(Generation,Haplotype,
                     OffspringCount:OffspringCountInMaleHeterozygote))
  adultSums <- summarise_at(adult, -1, sum)
  adultSexed <- mutate(reshapeBySex(adultSums),LifeHistoryStage="Adult",.before=1)
  
  offspring <- select(filter(selHap,LifeHistoryStage=="Adult"),
                  !c(Generation,Haplotype,
                     Count:MaleHeterozygote))
  offspringSums <- summarise_at(offspring, -1, sum)
  offspringSexed <- reshapeBySex(offspringSums)
  comb <- rbind(mutate(zygoteSexed, Offspring=offspringSexed$Count),
                mutate(adultSexed, Offspring=offspringSexed$Count))
  
  # Generate the per capita parentage/descendant allele count columns
  combPerCap <- mutate(comb,OffspringPerCapita=Offspring/Count)
  # print(combPerCap)
  return(combPerCap)
}


calcSummed <- function(hapTable) {
  AB <- cbind(Haplotype=c("AB"),combineReshapeHaps(hapTable,c("0:1;","0:1;0")))
  Ab <- cbind(Haplotype=c("Ab"),combineReshapeHaps(hapTable,c("1;","1;0")))
  aB <- cbind(Haplotype=c("aB"),combineReshapeHaps(hapTable,c("0;","0;0")))
  ab <- cbind(Haplotype=c("ab"),combineReshapeHaps(hapTable,c(";",";0")))
  combPerCapHapParentage <- rbind(AB,Ab,aB,ab)
  return(combPerCapHapParentage)
}

calcZygoteFreqs <- function(hapTable) {
  freqHapData <- select(filter(hapTable,
                            LifeHistoryStage=="Zygote",
                            Sex=="Both",
                            Zygosity=="Allele"),
                     -Sex)
  freqHapData <- mutate(freqHapData, Frequency = Count/sum(Count),.after=Count)
  # print(freqHapData)
  # print(sum(freqHapData$Count))
  return(freqHapData)
}

# To plot the frequencies of each haplotype calculated from the two-mutant scenarios
plotFrequencies <- function(combHapData) {
  # Reshape the haplotype data, calculate frequencies
  calcZygoteFreqs(combHapData)
  
  # Plot
  library(ggplot2)
  p <-ggplot(plotData, aes(Haplotype, Frequency, width=.5))
  p <- p +geom_col()
  p <- p +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p +scale_y_continuous(breaks=c(0.25,0.5,0.75,1.0))
  p <- p +ylab("Average Frequency")
  p <- p +coord_cartesian(ylim=c(0,1))
  p <- p +ggtitle("Intermediate Haplotypes are Removed")
  ggsave("Intermediate Haplotypes are Removed.png", 
         units="in", width=4, height=3, dpi=600, device = 'png')
}


calcR2 <- function(hapData) {
  pAB <- filter(hapData,Haplotype=="AB")$Frequency
  print(pAB)
  pAb <- filter(hapData,Haplotype=="Ab")$Frequency
  print(pAb)
  paB <- filter(hapData,Haplotype=="aB")$Frequency
  print(paB)
  pab <- filter(hapData,Haplotype=="ab")$Frequency
  print(pab)
  pA <- pAB+pAb
  print(pA)
  pB <- pAB+paB
  print(pB)
  D <- pAB - pA*pB
  print(D)
  r2 <- D^2/(pA*(1-pA)*pB*(1-pB))
  return(r2)
}


# Run the plot script for the two-mutation analysis

# Default variables
wd = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/Paternity/"
resultsDir = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Paternity/"
twoMutRes = "TwoMutPaternity/"
spacing = c(0.5375)
singMutRes = "SingMutPaternity/"
numReps = 1000

# Set working directory
setwd(wd)
twoMutResDir = paste0(resultsDir,twoMutRes)


if (!exists("twoMutHapPatTable")){
  if (file.exists(paste0(resultsDir,"twoMutHapPatTable.Rdata"))){
    load(paste0(resultsDir,"twoMutHapPatTable.Rdata"))
  } else {
    twoMutHapPatTable <- genHapTable(twoMutResDir,spacing,numReps,0.5,TRUE)
    
    save(twoMutHapPatTable,file=paste0(resultsDir,"twoMutHapPatTable.Rdata"))
  }
}

# Plot the fitnesses of haplotypes by sex
plotAlleleDescPerHap(calcDescendantsPerAllelle(twoMutHapPatTable))
# Plot the frequencies of each haplotype
plotFrequencies(calcSummed(twoMutHapPatTable))



# Run the plot script for the single-mutation analysis
library(tidyverse)

# Default variables
wd = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/Paternity/"
resultsDir = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Paternity/"
singMutRes = "SingMutPaternity/"

# Varied inputs to the single mutation sims
numReps = 1000
# repNums = seq(0,999,1)
encNums = c(2,3,4,10,100)
repQuals = seq(0,2,.1)
freqs = seq(.05,1,.05)

#test params
encNums = c(2,10)
repQuals = c(1,1.5)
freqs = c(.4,.6)

# Non-varied input variables
surQuals = 1
ph = 1
popSizes = 1000
nmc = FALSE
nfc = FALSE

# Parameter input list for parallel reading of sim tables
inParLs = list(singMutResDir,surQuals,repQuals,freqs,ph,nmc,nfc,encNums,popSizes,numReps)

# Set working directory
setwd(wd)
singMutResDir = paste0(resultsDir,singMutRes)


#Collect the single mutation sim data in parallel
library(parallel)

coreProportion = 0.6
minCores = 2

checkSingMutPatFullTable <- function(inParLs){
  if (!exists("singMutPatFullTable")){
    if (file.exists(paste0(resultsDir,"singMutPatFullTable.Rdata"))){
      load(paste0(resultsDir,"singMutPatFullTable.Rdata"))
    } else {
      singMutPatFullTable <- parallelExtract(getMutDataAcrossReps,inParLs,finishingFunc=NULL)
      
      save(singMutPatFullTable,file=paste0(resultsDir,"singMutPatFullTable.Rdata"))
    }
  }
  return(singMutPatFullTable)
}

singMutPatFullTable <- checkSingMutPatFullTable(inParLs)


# Plot a heatmap of paternity proportion values across tested parameters



### Functions associated with the single mutation simulation analyses

# For reading a SAIsim mutation file and returning only the final adult and zygote generation data
readFinalGenMutDataFile <- function(fileName) {
  mutCounts <- read.table(fileName,sep='\t',header=TRUE,na.strings = c("None","NA"))
  finalGen <- tail(mutCounts$Generation,1)
  finalGenData <- filter(mutCounts,Generation==finalGen)
  finalAdultGen <- finalGen - 1
  finalAdultGenData <- filter(mutCounts,Generation==finalAdultGen,LifeHistoryStage=="Adult")
  mutData <- rbind(finalAdultGenData,finalGenData)
  return(mutData)
}

# For reading a SAIsim mutation file and returning all data
readMutDataFile <- function(fileName) {
  mutCounts <- read.table(fileName,sep='\t',header=TRUE,na.strings = c("None","NA"))
  return(mutCounts)
}

# For reading all replicate SAIsim mutation files for a given parameter combination
getMutDataAcrossReps <- function(filePrefix,survEffect,reprEffect,frequency,proportionHeterozygous,
                       noMaleCosts,noFemaleCosts,encounterNum,popSize,numReps) {
  s <- format(survEffect,nsmall = 3)
  r <- format(reprEffect,nsmall = 3)
  f <- format(frequency,nsmall = 3)
  d <- format(proportionHeterozygous,nsmall = 3)
  nmc <- as.integer(noMaleCosts)
  nfc <- as.integer(noFemaleCosts)
  m <- formatC(encounterNum, width = 3, format = "d", flag = "0")
  N <- formatC(popSize, width = 3, format = "d", flag = "0")
  #print(s)
  #print(r)
  fileBase <- paste(filePrefix,"s",s,"r",r,"f",f,"d",d,"nmc",nmc,"nfc",nfc,"m",m,"N",N,"n",sep="")
  print(fileBase)
  mutTable <- tibble()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,".txt",sep="")
      fileData <- readMutDataFile(fileName)
      mutTable <- rbind(mutTable,
                        mutate(fileData,
                               SurvivalQuality=survEffect,ReproductiveQuality=reprEffect,Frequency=frequency,
                               ProportionHeterozygous=proportionHeterozygous,MaleCosts=!noMaleCosts,
                               FemaleCosts=!noFemaleCosts,EncounterNumber=encounterNum,PopulationSize=popSize,
                               Replicate=num,.before=1))
    })
  }
  #print(mutTable)
  return(mutTable)
}

# For reading all SAIsim mutation files across given parameter combinations
genMutDataTable <- function(filePrefix,surEffects,repEffects,frequencies,proportionsHeterozygous,
                            noMaleCosts,noFemaleCosts,encounterNumbers,popSizes,numReps) {
  mutTable <- tibble()
  for (s in surEffects) {
    for (r in repEffects) {
      for (f in frequencies) {
        for (d in proportionsHeterozygous) {
          for (nmc in noMaleCosts) {
            for (nfc in noFemaleCosts) {
              for (m in encounterNumbers) {
                for (N in popSizes) {
                  mutTable <- rbind(mutTable,
                                    getMutDataAcrossReps(filePrefix,s,r,f,d,nmc,nfc,m,N,numReps))
                }
              }
            }
          }
        }
      }
    }
  }
  return(mutTable)
}

extractPaternities <- function(mutTable) {
  
  return(propPoly)
}


# Running data collection and plotting in parallel


parallelPlot <- function(plotFunc,inParLs,finishingFunc=NULL) {
  maxReqCores <- length(inParLs[[1]])
  print(maxReqCores)
  numGuessCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  numCores <- min(numGuessCores,maxReqCores)
  parLs <- c(plotFunc,inParLs,
             SIMPLIFY = FALSE,
             mc.cores = numCores)
  print(parLs)
  plotReturns <- do.call(mcmapply,parLs)
  print(plotReturns)
  if (!is.null(finishingFunc)){
    finishingFunc(plotReturns)
  }
  return()
}





distParLsRecursion <- function(i,inParLs,outParLs,currPars,singlePars) {
  if (i > length(inParLs)) {
    for (j in seq(length(currPars))) {
      if (j > length(outParLs)) {
        outParLs <- append(outParLs,c(currPars[[j]]))
      } else {
        if (!singlePars[[j]]){
          outParLs[[j]] <- c(outParLs[[j]],currPars[[j]])
        }
      }
    }
  } else {
    for (p in inParLs[[i]]) {
      outParLs <- distParLsRecursion(i+1,inParLs,outParLs,c(currPars,list(p)),singlePars)
    }
  }
  return(outParLs)
}

distributedParamList <- function(inParLs) {
  if (length(inParLs) == 0){
    return(inParLs)
  }
  singlePars <- list()
  for (par in inParLs) {
    singlePars <- c(singlePars,length(par)==1)
  }
  outParLs <- distParLsRecursion(1,inParLs,list(),list(),singlePars)
  return(outParLs)
}

parallelExtract <- function(extractFunc,inParLs,finishingFunc=NULL) {
  maxReqCores <- length(inParLs[[1]])
  numGuessCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  numCores <- min(numGuessCores,maxReqCores)
  parLs <- c(extractFunc,distributedParamList(inParLs),
             SIMPLIFY = FALSE,
             mc.cores = numCores)
  print(parLs)
  extractResults <- do.call(mcmapply,parLs)
  # print(extractResults)
  dataTable <- Reduce(rbind,extractResults)
  # print(dataTable)
  if (!is.null(finishingFunc)){
    dataTable <- finishingFunc(dataTable)
  }
  return(dataTable)
}


getTableParamLs <- function(statTable,parVect){
  parLs <- list()
  for (name in parVect){
    parLs <- c(parLs,list(unique(statTable[[name]])))
  }
  return(parLs)
}

getSimParamLs <- function(statTable,includeReps=FALSE,includeGens=FALSE){
  expectedNames <- c("SurvivalQuality","ReproductiveQuality","Frequency",
                     "ProportionHeterozygous","MaleCosts","FemaleCosts",
                     "EncounterNumber","PopulationSize")
  if (includeReps) {
    expectedNames <- c(expectedNames,"Replicate")
  }
  if (includeGens) {
    expectedNames <- c(expectedNames,"Generation")
  }
  inParLs <- getTableParamLs(statTable,expectedNames)
  return(inParLs)
}

getFilteredTableLs <- function(tableLs,groupByVect){
  groupParLs <- getTableParamLs(tableLs[[1]],groupByVect)
  # print(groupParLs)
  groupParGrid <- do.call(expand.grid,groupParLs)
  # print(groupParGrid)
  filteredTableLs <- list()
  for (i in seq(length(tableLs))){
    filteredTableLs <- c(filteredTableLs,list(list()))
  }
  # print(filteredTableLs)
  for (i in seq(length(groupParGrid[[1]]))){
    subLs <- list()
    for (j in seq(length(groupByVect))){
      subLs <- c(subLs,list(substitute(name==val,
                                       list(name=sym(groupByVect[j]),
                                            val=groupParGrid[i,j]))))
      # print(groupParGrid[i,j])
    }
    # print(subLs)
    for (k in seq(length(tableLs))){
      specParFilteredTable <- do.call(filter,c(list(tableLs[[k]]),subLs))
      filteredTableLs[[k]] <- c(filteredTableLs[[k]],list(specParFilteredTable))
    }
  }
  return(filteredTableLs)
}

parallelAnalysis <- function(calcFunc, statTableLs, inParLs=list(), finishingFunc=NULL) {
  numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  parLs <- c(calcFunc,statTableLs,distributedParamList(inParLs),
             SIMPLIFY = FALSE,
             mc.cores = numCores)
  # print(parLs)
  extractResults <- do.call(mcmapply,parLs)
  
  calcTable <- Reduce(rbind,extractResults)
  if (!is.null(finishingFunc)){
    calcTable <- finishingFunc(calcTable)
  }
  return(calcTable)
}



