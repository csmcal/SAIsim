# SAIsim EffectGrid Heatmap Generator Script 
#


# Example run: Rscript Heatmap_SingleLocus_marula.R ./dmelN3/


# Collect input arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  # stop("At least one argument must be supplied (input file).n", call.=FALSE)
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


# Source the plotting and data collection functions

source("Heatmap_SingleLocus.R")


# Running data collection and plotting in parallel

inParLs <- list(dataDirList,prefixList,list(surEffects),list(repEffects),popSize,numReps)
print("Plotting proportion polymorphic")
parallelPlot(plotPolyHeatmap,inParLs)
print("Plotting average final frequency")
parallelPlot(plotFreqHeatmap,inParLs)


