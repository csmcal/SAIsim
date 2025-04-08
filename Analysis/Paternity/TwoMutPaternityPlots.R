# PLotting Per-capita or relative Paternity in SAIsim results
#   for simulations with two variants linked by an inversion
# For generating plots of two-mutation simulations
library(tidyverse)



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
plotAlleleDescPerHap <- function(combDescPerAllele,filename) {
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
  p <- p +scale_fill_manual(values = c(rep(c("gray75", "gray55", "lightblue"))))
  # p <- p +scale_fill_manual(values = c(rep(c('#228833', '#4477AA', '#CCBB44'))))
  p <- p +scale_y_continuous(breaks=c(0.75,1.0,1.25))
  p <- p +ylab("Average Descendants per Allele")
  p <- p +coord_cartesian(ylim=c(.65,1.3))
  p <- p +ggtitle("Haplotypes Exhibit Synergistic Epistasis")
  ggsave(filename, 
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


calcSummed <- function(hapTable,hapList) {
  AB <- cbind(Haplotype=hapList[[1]][[1]],combineReshapeHaps(hapTable,hapList[[1]][[2]]))
  Ab <- cbind(Haplotype=hapList[[2]][[1]],combineReshapeHaps(hapTable,hapList[[2]][[2]]))
  aB <- cbind(Haplotype=hapList[[3]][[1]],combineReshapeHaps(hapTable,hapList[[3]][[2]]))
  ab <- cbind(Haplotype=hapList[[4]][[1]],combineReshapeHaps(hapTable,hapList[[4]][[2]]))
  combPerCapHapParentage <- rbind(AB,Ab,aB,ab)
  return(combPerCapHapParentage)
}

# calcSummed <- function(hapTable) {
#   AB <- cbind(Haplotype=c("AB"),combineReshapeHaps(hapTable,c("0:1;","0:1;0")))
#   Ab <- cbind(Haplotype=c("Ab"),combineReshapeHaps(hapTable,c("1;","1;0")))
#   aB <- cbind(Haplotype=c("aB"),combineReshapeHaps(hapTable,c("0;","0;0")))
#   ab <- cbind(Haplotype=c("ab"),combineReshapeHaps(hapTable,c(";",";0")))
#   combPerCapHapParentage <- rbind(AB,Ab,aB,ab)
#   return(combPerCapHapParentage)
# }

calcFreqs <- function(hapTable,LHS="Zygote") {
  freqHapData <- select(filter(hapTable,
                               LifeHistoryStage==LHS,
                               Sex=="Both",
                               Zygosity=="Allele"),
                        -Sex)
  freqHapData <- mutate(freqHapData, Frequency = Count/sum(Count),.after=Count)
  # print(freqHapData)
  # print(sum(freqHapData$Count))
  return(freqHapData)
}

# To plot the frequencies of each haplotype calculated from the two-mutant scenarios
plotFrequencies <- function(combHapData,filename) {
  # Reshape the haplotype data, calculate frequencies
  plotData <- calcFreqs(combHapData,LHS="Zygote")
  
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
  ggsave(filename,
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
twoMutNoInvRes = "TwoMutPatNoInv/"
spacingNoInv = c(0.25)
singMutRes = "SingMutPaternity/"
numReps = 1000


# Set working directory
setwd(wd)
twoMutResDir = paste0(resultsDir,twoMutRes)
twoMutNoInvResDir = paste0(resultsDir,twoMutNoInvRes)

#Inv/non parameters
hapInv=list(list(c("AB"),c("0:1;","0:1;0")),
            list(c("Ab"),c("1;","1;0")),
            list(c("aB"),c("0;","0;0")),
            list(c("ab"),c(";",";0")))

hapNoInv=list(list(c("AB"),c("0:1;")),
            list(c("Ab"),c("1;")),
            list(c("aB"),c("0;")),
            list(c("ab"),c(";")))

if (!exists("twoMutHapPatTable")){
  if (file.exists(paste0(resultsDir,"twoMutHapPatTable.Rdata"))){
    load(paste0(resultsDir,"twoMutHapPatTable.Rdata"))
  } else {
    twoMutHapPatTable <- genHapTable(twoMutResDir,spacing,numReps,0.5,TRUE)
    
    save(twoMutHapPatTable,file=paste0(resultsDir,"twoMutHapPatTable.Rdata"))
  }
}

# Plot the fitnesses of haplotypes by sex
plotAlleleDescPerHap(calcDescendantsPerAllelle(twoMutHapPatTable),"Haplotypes Exhibit Synergistic Epistasis.png")
# Plot the frequencies of each haplotype
plotFrequencies(calcSummed(twoMutHapPatTable,hapInv),"Intermediate Haplotypes are Removed.png")


# Without the inversion
if (!exists("twoMutHapPatNoInvTable")){
  if (file.exists(paste0(resultsDir,"twoMutHapPatNoInvTable.Rdata"))){
    load(paste0(resultsDir,"twoMutHapPatNoInvTable.Rdata"))
  } else {
    twoMutHapPatNoInvTable <- genHapTable(twoMutNoInvResDir,spacingNoInv,numReps,0.5,FALSE)
    
    save(twoMutHapPatNoInvTable,file=paste0(resultsDir,"twoMutHapPatNoInvTable.Rdata"))
  }
}

# Plot the fitnesses of haplotypes by sex
descPerAlleleNoInv <- calcDescendantsPerAllelle(twoMutHapPatNoInvTable)
plotAlleleDescPerHap(descPerAlleleNoInv,"Haplotypes Exhibit Synergistic Epistasis Without Inversions.png")
# Plot the frequencies of each haplotype
plotFrequencies(calcSummed(twoMutHapPatNoInvTable,hapInv),"Intermediate Haplotypes are Removed Without Inversions.png")

