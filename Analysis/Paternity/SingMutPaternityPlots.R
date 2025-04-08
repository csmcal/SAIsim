# PLotting Per-capita or relative Paternity in SAIsim results for single variant populations
# For generating plots of two-mutation simulations
library(tidyverse)
library(parallel)


# Example run: Rscript SingMutPaternityPlots.R ./SingMutPaternity/


# Default variables
wd = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/Paternity/"
resultsBase = "/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Paternity/"
singMutRes = "SingMutPaternity/"
resultsDir = paste0(resultsBase,singMutRes)

# Varied inputs to the single mutation sims
numReps = 1000
# repNums = seq(0,999,1)
encNums = c(2,3,4,10,100)
repQuals = seq(0,2,.1)
freqs = seq(.05,1,.05)

# #test params
# encNums = c(2,10)
# repQuals = c(1,1.5)
# freqs = c(.4,.6)

# Non-varied input variables
surQuals = 1
ph = 1
popSizes = 1000
nmc = FALSE
nfc = FALSE


# Collect input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  print("Defaulting to using:")
  print(paste0("resultsDir = ",resultsDir))
} else if (length(args)==1) {
  resultsDir <- args[1]
  print("Using input results directory:")
  print(paste0("resultsDir = ",resultsDir))
  print("For sim parameters, defaulting to using:")
} else {
  stop(paste0("Argument interpretation not implemented yet for > 2 args; ",
              length(args)," args given"), call.=FALSE)
}
print("Survival Qualities = ")
print(surQuals)
print("Reproductive Qualities = ")
print(repQuals)
print("Initial Frequencies = ")
print(freqs)
print("Proportions Heterozygous = ")
print(ph)
print("No Male Costs = ")
print(nmc)
print("No Female Costs = ")
print(nfc)
print("Encounter Numbers = ")
print(encNums)
print(paste0("Population Sizes = ",popSizes))
print(paste0("Number of Replicates = ",numReps))

if (substr(resultsDir,nchar(resultsDir),nchar(resultsDir))!="/") {
  resultsDir <- paste0(resultsDir,'/')
}

# Parameter input list for parallel reading of sim tables
inParLs = list(resultsDir,surQuals,repQuals,freqs,ph,nmc,nfc,encNums,popSizes,numReps)


# Variables for collecting the single mutation sim data in parallel
library(parallel)

coreProportion = 0.6
minCores = 2


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
    # try({
    fileName <- paste(fileBase,n,".txt",sep="")
    fileData <- readMutDataFile(fileName)
    mutTable <- rbind(mutTable,
                      mutate(fileData,
                             SurvivalQuality=survEffect,ReproductiveQuality=reprEffect,Frequency=frequency,
                             ProportionHeterozygous=proportionHeterozygous,MaleCosts=!noMaleCosts,
                             FemaleCosts=!noFemaleCosts,EncounterNumber=encounterNum,PopulationSize=popSize,
                             Replicate=num,.before=1))
    # })
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
  maxReqCores <- max(sapply(inParLs,FUN=length))
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


extractPatPerCapOnly <- function(mutTable) {
  adultMaleAlleleCountTable <- select(filter(mutTable,
                                             Sex=="Male",
                                             LifeHistoryStage=="Adult",
                                             CountType=="Allele"),
                                      !c(CountType,Homozygote, Heterozygote))
  offspringCountTable <- select(filter(mutTable,
                                       Sex=="Male",
                                       LifeHistoryStage=="Adult",
                                       CountType=="Offspring"),
                                !c(CountType,Homozygote, Heterozygote))
  patPerCapTable <- mutate(adultMaleAlleleCountTable,
                           Offspring=offspringCountTable$Count,
                           OffspringPerCap=offspringCountTable$Count/Count)
  return(patPerCapTable)
}



patGridHeatmap <- function(avgPatPerCapTable) {
  p <- ggplot(avgPatPerCapTable,aes(Frequency, ReproductiveQuality, fill=AvgOffspringPerCap)) + 
    geom_tile() +
    labs(title=NULL,x=NULL,y=NULL) +
    theme(legend.position="none") + 
    facet_wrap(vars(EncounterNumber)) +
    scale_fill_gradientn(colors=paul_tol_modified_rev) +
    guides(colour = guide_colourbar(title.position = "right"))+
    theme(legend.title = element_text(size = 10, angle = 90),
          legend.title.align = 0.5,
          # legend.key.height = unit(height/10, "in"),
          # legend.key.align = 0.5,
          legend.position = "right",
          legend.direction = "vertical"
    )
  return(p)
}


patGridLine <- function(avgPatPerCapTable,repQual=1) {
  p <- ggplot(filter(avgPatPerCapTable, 
                     EncounterNumber %in% c(2,3,4,10,100),
                     ReproductiveQuality %in% c(0,.1,.2,.3,.4,.5)),
              aes(Frequency, AvgOffspringPerCap, group = ReproductiveQuality, color = ReproductiveQuality)) +
    geom_line()+
    geom_point() +
    geom_tile() +
    labs(title=NULL,x=NULL,y=NULL) +
    # theme(legend.position="none") + 
    facet_wrap(vars(EncounterNumber)) + 
    scale_color_gradientn(colors=paul_tol_modified_rev) +
    guides(colour = guide_colourbar(title.position = "right"))+
    theme(legend.title = element_text(size = 10, angle = 90),
          legend.title.align = 0.5,
          # legend.key.height = unit(height/10, "in"),
          # legend.key.align = 0.5,
          legend.position = "right",
          legend.direction = "vertical"
    )
  return(p)
}



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


### Collect the data from results file, or skip if already tabulated
singMutPatFullTable <- checkSingMutPatFullTable(inParLs)

### Process the data to per-capita paternity values
patPerCapTable <- extractPatPerCapOnly(singMutPatFullTable)
avgPatPerCapTable <- summarize(group_by(patPerCapTable,across(c(SurvivalQuality:PopulationSize,Generation:Sex))),
                               AvgCount=mean(Count),
                               AvgOffspring=mean(Offspring),
                               AvgOffspringPerCap=mean(OffspringPerCap))

### Plot heatmaps of paternity proportion values across tested parameters


freq_palette <- colorRampPalette(c("cornflowerblue","green", "yellow"))(n = 3999)
basic_freq_palette <- c("cornflowerblue","green", "yellow")
poly_palette <- colorRampPalette(c("white","cornflowerblue"))(n = 3999)
paul_tol_default <- c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499')
paul_tol_ordered <- c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499')
paul_tol_modified <- c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 'white')
paul_tol_modified_rev <- rev(paul_tol_modified)
paul_tol_NA <- c('#DDDDDD')


