# SAIsim Mutation Frequency by Spacing Script (with and w/o Inversions)
# For generating plots of two-mutation simulations

library(tidyverse)


getFinalFreqMutInv <- function(spacing,mutFileName,invFileName,hasInv) {
  mutFreq <- read.table(mutFileName,sep='\t',header=TRUE)
  lastMutFreqs <- unlist(tail(mutFreq[,4:5],n=1))
  finalGenData <- c(spacing,lastMutFreqs)
  if (hasInv) {
    invFreq <- read.table(invFileName,sep='\t',header=TRUE)
    lastInvFreq <- unlist(tail(invFreq[,3],n=1))
    finalGenData <- c(spacing,lastMutFreqs,lastInvFreq)
  }
  return(finalGenData)
}

getFreqsAtSpa <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  s <- format(spacing,nsmall = 4)
  f <- format(initFreq,nsmall = 3)
  i <- as.integer(hasInv)
  # print(s)
  # fileBase <- paste(filePrefix,"spacing",s,"n", sep="")
  fileBase <- paste(filePrefix,'pos',s,'inv',i,'freq',f,'n', sep="")
  # print(fileBase)
  finalFreqs <- cbind()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      # print(paste0("Trying ",fileBase,n))
      invFileName <- paste(fileBase,n,"Inv.txt",sep="")
      mutFileName <- paste(fileBase,n,"Mut.txt",sep="")
      # hapFileName <- paste(fileBase,n,"Hap.txt",sep="")
      finalFreqs <- rbind(finalFreqs,getFinalFreqMutInv(spacing,mutFileName,invFileName,hasInv))
    })
  }
  # print(finalFreqs)
  return(finalFreqs)
}

genFreqSpaTable <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  # spacing <- seq(0.0,0.35,by=0.025)
  # repEffects <- seq(0.00,1.0,by=0.05)
  freqSpaTable <- tibble()
  # summaryTable <- cbind()
  for (s in spacing) {
    freqSpaTable <- rbind(freqSpaTable,getFreqsAtSpa(filePrefix,s,numReps,initFreq,hasInv))
    # for m in 2:ncol(freqSpaTable) {
    #   
    # }
  }
  # freqSpaTable <- as.data.frame(freqSpaTable)
  # freqSpaTable <- as_tibble(freqSpaTable)
  colnames(freqSpaTable)[1]<-"Spacing"
  colnames(freqSpaTable)[2]<-"M1"
  colnames(freqSpaTable)[3]<-"M2"
  if (hasInv) {
    colnames(freqSpaTable)[4]<-"I1"
  }
  # print(freqSpaTable)
  return(freqSpaTable)
}



getFinalHapFreqs <- function(spacing,hapFileName) {
  hapFreqTable <- read.table(hapFileName,sep='\t',header=TRUE)
  finalGen <- hapFreqTable$Generation[-1]
  finalGenHapData <- filter(hapFreqTable,Generation==finalGen)
  return(finalGenHapData)
}

getHapFreqsAtSpa <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  s <- format(spacing,nsmall = 4)
  f <- format(initFreq,nsmall = 3)
  i <- as.integer(hasInv)
  fileBase <- paste(filePrefix,'pos',s,'inv',i,'freq',f,'n', sep="")
  finalHapFreqs <- tibble()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      # print(paste0("Trying ",fileBase,n))
      hapFileName <- paste(fileBase,n,"Hap.txt",sep="")
      finalHapFreqs <- rbind(finalHapFreqs,getFinalHapFreqs(spacing,hapFileName))
    })
  }
  return(finalHapFreqs)
}

genHapFreqSpaTable <- function(filePrefix,spacing,numReps,initFreq,hasInv) {
  freqSpaTable <- tibble()
  for (s in spacing) {
    freqSpaTable <- rbind(freqSpaTable,getHapFreqsAtSpa(filePrefix,s,numReps,initFreq,hasInv))
  }
  # freqSpaTable <- as.tibble(freqSpaTable)
  return(freqSpaTable)
}




# To retrieve polymorphism data from
summaryFixLoss <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                           .groups="drop") {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  dataSumm <- data %>% group_by(across(all_of(groupvars))) %>%
    summarize(N     = length2(.data[[measurevar]], na.rm=na.rm),
              Nfix  = length2(subset(.data[[measurevar]],.data[[measurevar]]==2000), na.rm=na.rm),
              Nloss = length2(subset(.data[[measurevar]],.data[[measurevar]]==0), na.rm=na.rm),
              .groups=.groups)
  
  
  dataSumm$Pfix <- dataSumm$Nfix/dataSumm$N
  dataSumm$Ploss <- dataSumm$Nloss/dataSumm$N
  dataSumm$Peither <- dataSumm$Pfix + dataSumm$Ploss
  dataSumm$Ppoly <- rep(1.0,length(dataSumm$Peither)) - dataSumm$Peither
  
  return(dataSumm)
}


# Plot generator for a single-panel plot of the proportion remaining polymorphic
plotPolyCombined <- function(freqSpaTable,freqSpaInvTable){
  
  #Combining the data to one table
  polyMutNoInv <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                 measurevar = "M1")
  polyMutInv <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                               measurevar = "M1")
  polyInv <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                            measurevar = "I1")
  
  polyMutNoInv <- cbind(polyMutNoInv,Locus=c("M"))
  polyMutInv <- cbind(polyMutInv,Locus=c("MI"))
  polyInv <- cbind(polyInv,Locus=c("I"))
  
  polyData <- as.data.frame(rbind(polyMutNoInv,polyMutInv,polyInv))[,c("Spacing","Ppoly","Locus")]
  # print(polyData)
  
  polyData$Locus <- factor(polyData$Locus,levels=c("M","MI","I"))
  
  # Plotting using ggplot2
  library(ggplot2)
  p <- ggplot() + 
    geom_line(data=polyData,aes(x=Spacing, y=Ppoly, group=Locus, linetype=Locus)) +
    xlab("Chromosomal Position of the Lesser SA Locus (M)") +
    ylab("Proportion of Sims Retaining Polymorphism") + 
    scale_x_continuous(breaks=seq(0, 1, .1), minor_breaks=seq(0, 1, .05)) + 
    scale_y_continuous(breaks=seq(0, 1, .25), minor_breaks=seq(0, 1, .05)) + 
    # coord_cartesian(ylim=c(0,1.0), xlim=c(0,.35)) +
    scale_linetype_manual(values = c(3,4,1),
                     labels=c("Smaller SA Locus (w/o Inv)",
                              "Smaller SA Locus (w/ Inv)",
                              "Inversion (Breakpoints 0.21, 0.55)")) +
    ggtitle("Persistence of Polymoprhism in SAI Simulations") +
    theme_bw() +
    theme(legend.justification=c(1,0),
      legend.position=c(.97,0.6),  # Position legend in bottom right
      legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
      legend.title = element_text(size=10),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=1.9, size=17, face="bold"))
  ggsave("Small Effect Locus & Inversion Poly Freq Fine.png", units="in", width=6, height=4, dpi=600, device = 'png')
}

spacingPropFixPlot <- function(freqSpaTable,
                               fileName="Polymorphism in Simulations of Two SA Mutations.png",
                               hasInv=FALSE,
                               plotLegend=FALSE) {
  polyTableM1 <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                  measurevar = "M1")
  polyTableM2 <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                  measurevar = "M2")
  polyTableM1 <- cbind(polyTableM1,Locus=c("M1"))
  polyTableM2 <- cbind(polyTableM2,Locus=c("M2"))
  
  polyData <- as_tibble(rbind(polyTableM1,polyTableM2))[,c("Spacing","Ppoly","Locus")]
  polyData$Locus <- factor(polyData$Locus,levels=c("M1","M2"))
  
  if (hasInv) {
    polyTableI <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                     measurevar = "I1")
    polyTableI <- cbind(polyTableI,Locus=c("I"))
    polyData <- as_tibble(rbind(polyTableM1,polyTableM2,polyTableI))[,c("Spacing","Ppoly","Locus")]
    polyData$Locus <- factor(polyData$Locus,levels=c("M1","M2","I"))
  }
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  p <- ggplot() + 
    geom_line(data=polyData,aes(x=Spacing, y=Ppoly, group=Locus, colour=Locus)) + #, colour=c("lightcoral","red4","steelblue2")) +
    xlab("Chromosomal Position of the Lesser SA Locus (M)") +
    ylab("Sims Retaining Polymorphism") +
    # coord_cartesian(ylim = c(0, 1.0),xlim = c(0, 1.0)) +
    scale_x_continuous(breaks=seq(0, 1, .1), minor_breaks=seq(0, 1, .05)) + 
    scale_y_continuous(breaks=seq(0, 1, .25), minor_breaks=seq(0, 1, .05)) + 
    # scale_colour_hue(name="Inversion",    # Legend label, use darker colors
    #                  breaks=c("1", "2"),
    #                  labels=c("Inv 0.21-0.55", "Mut 2"),
    #                  l=40) +                    # Use darker colors, lightness=40
    # ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    # expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw()
  if (hasInv) {
    p <- p +
      scale_colour_manual(values = c("lightcoral","red4","steelblue2"),
                          labels=c("Less SA Locus (Variable)",
                                   "More SA Locus (at 0.225M)",
                                   "Inversion (0.21 to 0.55M)"))
    if (plotLegend) {
      p <- p + theme(legend.justification=c(1,0),
                     legend.position=c(.97,0.5),  # Position legend in bottom right
                     legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
                     legend.title = element_text(size=10),
                     legend.text = element_text(size=7),
                     plot.title = element_text(hjust=1.9, size=17, face="bold"))
    } else {
      p <- p + theme(legend.position="none")
    }
    
  } else {
    p <- p +
      scale_colour_manual(values = c("lightcoral","red4"),
                          labels=c("Less SA Locus (Variable)",
                                   "More SA Locus (at 0.225M)"))
    if (plotLegend) {
      p <- p + theme(legend.justification=c(1,0),
                     legend.position=c(.97,0.5),  # Position legend in bottom right
                     legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
                     legend.title = element_text(size=10),
                     legend.text = element_text(size=7),
                     plot.title = element_text(hjust=1.9, size=17, face="bold"))
    } else {
      p <- p + theme(legend.position="none")
    }
  }
  ggsave(fileName, units="in", width=5, height=3, dpi=300, device = 'png')
}


spacingInvPropFixPlot <- function(freqSpaInvTable) {
  fixLossTable <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                                 measurevar = "I1")
  # print(fixLossTable)
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  p <- ggplot(fixLossTable) + 
    geom_line(aes(x=Spacing, y=Ppoly, group=1), colour="steelblue2") +
    # geom_line(aes(x=Spacing, y=Pfix, group=1), colour="blue") +
    # geom_line(aes(x=Spacing, y=Ploss, group=1), colour="red") +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Percent of Simulations") +
    coord_cartesian(ylim = c(0, 1.0),xlim = c(0, 1.0)) +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Inv 0.21-0.55", "Mut 2"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
  ggsave("Inversion Polymorphism.png", units="in", width=4, height=3, dpi=300, device = 'png')
}




# Run the plot script

# Default variables
wd = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/TwoMutSpacing/"
resultsDir = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/TwoMutSpacing/"
stdRes = "SpacingNoInv/"
invRes = "SpacingInv/"
rareStdRes = "SpacingRareNoInv/"
rareInvRes = "SpacingRareInv/"
partRareStdRes = "SpacingEqBRareNoInv/"
partRareInvRes = "SpacingEqBRareInv/"
numReps = 1000
spacing = seq(0, 1, .005)
  

# Collect input arguments
args = commandArgs(trailingOnly=TRUE)
# Replace position-defined input arguments
if (length(args)>0) {
  wd = args[1]
  print(paste0("Defaulting to using resultsDir = ",resultsDir))
  print(paste0("Defaulting to using numReps = ",numReps))
} else if (length(args)>1) {
} else {
  resultsDir = args[1]
  numReps = args[2]
  print(paste0("Using input resultsDir = ",resultsDir))
  print(paste0("Using input numReps = ",numReps))
}

if (substr(resultsDir,nchar(resultsDir),nchar(resultsDir))!="/") {
  resultsDir <- paste0(resultsDir,'/')
}


# Set working directory
setwd(wd)

stdResDir = paste0(resultsDir,stdRes)
invResDir = paste0(resultsDir,invRes)


if (!exists("twoMutNoInvEqFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutNoInvEqFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutNoInvEqFreqSpaTable.Rdata"))
  } else {
    twoMutNoInvEqFreqSpaTable <- genFreqSpaTable(stdResDir,spacing,numReps,0.5,FALSE)
    
    save(twoMutNoInvEqFreqSpaTable,file=paste0(resultsDir,"twoMutNoInvEqFreqSpaTable.Rdata"))
  }
}

if (!exists("twoMutInvEqFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutInvEqFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutInvEqFreqSpaTable.Rdata"))
  } else {
    twoMutInvEqFreqSpaTable <- genFreqSpaTable(invResDir,spacing,numReps,0.5,TRUE)
    
    save(twoMutInvEqFreqSpaTable,file=paste0(resultsDir,"twoMutInvEqFreqSpaTable.Rdata"))
  }
}

# Plot a combined figure
plotPolyCombined(twoMutNoInvEqFreqSpaTable,twoMutInvEqFreqSpaTable)

spacingPropFixPlot(twoMutNoInvEqFreqSpaTable)
spacingPropFixPlot(twoMutInvEqFreqSpaTable,
                   fileName = "Polymorphism in Simulations of Two SA Mutations and an Inversion.png",
                   hasInv = TRUE)

# spacingInvPropFixPlot(twoMutInvEqFreqSpaTable)



rareStdResDir = paste0(resultsDir,rareStdRes)
rareInvResDir = paste0(resultsDir,rareInvRes)


if (!exists("twoMutNoInvRareFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutNoInvRareFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutNoInvRareFreqSpaTable.Rdata"))
  } else {
    twoMutNoInvRareFreqSpaTable <- genFreqSpaTable(rareStdResDir,spacing,numReps,0.05,FALSE)
    
    save(twoMutNoInvRareFreqSpaTable,file=paste0(resultsDir,"twoMutNoInvRareFreqSpaTable.Rdata"))
  }
}

if (!exists("twoMutInvRareFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutInvRareFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutInvRareFreqSpaTable.Rdata"))
  } else {
    twoMutInvRareFreqSpaTable <- genFreqSpaTable(rareInvResDir,spacing,numReps,0.05,TRUE)
    
    save(twoMutInvRareFreqSpaTable,file=paste0(resultsDir,"twoMutInvRareFreqSpaTable.Rdata"))
  }
}


spacingPropFixPlot(twoMutNoInvRareFreqSpaTable,
                   fileName = "Polymorphism in Simulations of Two SA Mutations at Low Initial Freq.png")
spacingPropFixPlot(twoMutInvRareFreqSpaTable,
                   fileName = "Polymorphism in Simulations of Two SA Mutations and an Inversion at Low Initial Freq.png",
                   hasInv = TRUE)
                   # plotLegend = TRUE)


partRareStdResDir = paste0(resultsDir,partRareStdRes)
partRareInvResDir = paste0(resultsDir,partRareInvRes)


if (!exists("twoMutNoInvEqBRareFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutNoInvEqBRareFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutNoInvEqBRareFreqSpaTable.Rdata"))
  } else {
    twoMutNoInvEqBRareFreqSpaTable <- genFreqSpaTable(partRareStdResDir,spacing,numReps,0.05,FALSE)
    
    save(twoMutNoInvEqBRareFreqSpaTable,file=paste0(resultsDir,"twoMutNoInvEqBRareFreqSpaTable.Rdata"))
  }
}

if (!exists("twoMutInvEqBRareFreqSpaTable")){
  if (file.exists(paste0(resultsDir,"twoMutInvEqBRareFreqSpaTable.Rdata"))){
    load(paste0(resultsDir,"twoMutInvEqBRareFreqSpaTable.Rdata"))
  } else {
    twoMutInvEqBRareFreqSpaTable <- genFreqSpaTable(partRareInvResDir,spacing,numReps,0.05,TRUE)
    
    save(twoMutInvEqBRareFreqSpaTable,file=paste0(resultsDir,"twoMutInvEqBRareFreqSpaTable.Rdata"))
  }
}


spacingPropFixPlot(twoMutNoInvEqBRareFreqSpaTable,
                   fileName = "Polymorphism in Simulations of Two SA Mutations with the Weaker at Low Initial Freq.png")
spacingPropFixPlot(twoMutInvEqBRareFreqSpaTable,
                   fileName = "Polymorphism in Simulations of Two SA Mutations with the Weaker and an Inversion at Low Initial Freq.png",
                   hasInv = TRUE)
                   # plotLegend = TRUE)

