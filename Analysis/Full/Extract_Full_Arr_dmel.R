# For extracting and analyzing full SAIsim model data, sourced from other scripts
#
# Example run: Rscript Extract_Full_Arr_dmel.R ./dmelN3arrangements/


# Collect input arguments
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  # stop("At least one argument must be supplied (input file).n", call.=FALSE)
  resultsDir = './'
  numReps = 1000
  print(paste0("Defaulting to using resultsDir = ",resultsDir))
  print(paste0("Defaulting to using numReps = ",numReps))
} else if (length(args)==1) {
  resultsDir = args[1]
  numReps = 1000
  print(paste0("Using input resultsDir = ",resultsDir))
  print(paste0("Defaulting to using numReps = ",numReps))
} else {
  resultsDir = args[1]
  numReps = args[2]
  print(paste0("Using input resultsDir = ",resultsDir))
  print(paste0("Using input numReps = ",numReps))
}

if (substr(resultsDir,nchar(resultsDir),nchar(resultsDir))!="/") {
  resultsDir <- paste0(resultsDir,'/')
}
# 
# resultsDir = '~/cMcAl/UWM/Pool/SAIsim/Results/Full/dmelN3subset/'
# numReps = 10
# setwd('~/cMcAl/UWM/Pool/SAIsim/Analysis/Full/')

# Source the extraction functions
source("Extract_Full.R")

#
#
# Extract, process, and save specific data when run
#


checkAllGenAllArrTable <- function(){
  if (!exists("allGenAllArrTable")){
    if (file.exists(paste0(resultsDir,"allGenAllArrTable.Rdata"))){
      load(paste0(resultsDir,"allGenAllArrTable.Rdata"))
    } else {
      # partialAllGenN3L1MbpRareMutAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                        c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e1),c(1.295e-1),
      #                                                                                        c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpRareMutAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                           c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e1),c(1.295e-1),
      #                                                                                           c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpRareMutAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                           c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e1),c(1.295e-1),
      #                                                                                           c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                         c(1.042e-6),c(1.042e-8,1.042e-6),c(10,100),c(4.36e-1),c(1.295e-1),
      #                                                                                         c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                            c(1.042e-6),c(1.042e-8,1.042e-6),c(10,100),c(4.36e-1),c(1.295e-1),
      #                                                                                            c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
      #                                                                                            c(1.042e-6),c(1.042e-8,1.042e-6),c(10,100),c(4.36e-1),c(1.295e-1),
      #                                                                                            c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L1MbpAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                      c(1.042e-2),c(1.042e-4,1.042e-2),c(10,100),c(4.36e1),c(1.295e-1),
                                                                                      c(1e3),c(4e4),c(1),c(0),c(0),seq(numReps)-1))
      print(paste0("partialAllGenN3L1MbpAllArrTable: ",ncol(partialAllGenN3L1MbpAllArrTable)," columns"))
      partialAllGenN3L1MbpAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                         c(1.042e-2),c(1.042e-4,1.042e-2),c(10,100),c(4.36e1),c(1.295e-1),
                                                                                         c(1e3),c(4e4),c(1),c(1),c(0),seq(numReps)-1))
      print(paste0("partialAllGenN3L1MbpAllArrTableNMC: ",ncol(partialAllGenN3L1MbpAllArrTableNMC)," columns"))
      partialAllGenN3L1MbpAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                         c(1.042e-2),c(1.042e-4,1.042e-2),c(10,100),c(4.36e1),c(1.295e-1),
                                                                                         c(1e3),c(4e4),c(1),c(0),c(1),seq(numReps)-1))
      print(paste0("partialAllGenN3L1MbpAllArrTableNFC: ",ncol(partialAllGenN3L1MbpAllArrTableNFC)," columns"))
      partialAllGenN3L10KbpAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                       c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e-1),c(1.295e-1),
                                                                                       c(1e3),c(4e4),c(1),c(0),c(0),seq(numReps)-1))
      print(paste0("partialAllGenN3L10KbpAllArrTable: ",ncol(partialAllGenN3L10KbpAllArrTable)," columns"))
      partialAllGenN3L10KbpAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                          c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e-1),c(1.295e-1),
                                                                                          c(1e3),c(4e4),c(1),c(1),c(0),seq(numReps)-1))
      print(paste0("partialAllGenN3L10KbpAllArrTableNMC: ",ncol(partialAllGenN3L10KbpAllArrTableNMC)," columns"))
      partialAllGenN3L10KbpAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                          c(1.042e-4),c(1.042e-6,1.042e-4),c(10,100),c(4.36e-1),c(1.295e-1),
                                                                                          c(1e3),c(4e4),c(1),c(0),c(1),seq(numReps)-1))
      print(paste0("partialAllGenN3L10KbpAllArrTableNFC: ",ncol(partialAllGenN3L10KbpAllArrTableNFC)," columns"))
      # allGenAllArrTable <- rbind(partialAllGenN3L1MbpRareMutAllArrTable,
      #                            partialAllGenN3L1MbpRareMutAllArrTableNMC,
      #                            partialAllGenN3L1MbpRareMutAllArrTableNFC,
      #                            partialAllGenN3L10KbpRareMutAllArrTable,
      #                            partialAllGenN3L10KbpRareMutAllArrTableNMC,
      #                            partialAllGenN3L10KbpRareMutAllArrTableNFC,
      allGenAllArrTable <- rbind(partialAllGenN3L1MbpAllArrTable,
                                 partialAllGenN3L1MbpAllArrTableNMC,
                                 partialAllGenN3L1MbpAllArrTableNFC,
                                 partialAllGenN3L10KbpAllArrTable,
                                 partialAllGenN3L10KbpAllArrTableNMC,
                                 partialAllGenN3L10KbpAllArrTableNFC)
      save(allGenAllArrTable,file=paste0(resultsDir,"allGenAllArrTable.Rdata"))
    }
  }
  return(allGenAllArrTable)
}

allGenAllArrTable <- checkAllGenAllArrTable()



checkAllGenAllArrPopStatsTable <- function(allGenAllArrTable){
  if (!exists("allGenAllArrPopStatsTable")){
    if (file.exists(paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))){
      load(paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))
    } else {
      filteredTableLs <- varCoreFilteredTableLs(list(allGenAllArrTable),maxSplitInd=12)
      allGenAllArrPopStatsTable <- parallelAnalysis(calcAllArrPopStats,filteredTableLs,list())
      # allGenAllArrPopStatsTable <- parallelAnalysis(calcAllArrPopStats,
      #                                               list(allGenAllArrTable),
      #                                               getSimParamLs(allGenAllArrTable,includeReps=TRUE))
      save(allGenAllArrPopStatsTable,file=paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))
    }
  }
  return(allGenAllArrPopStatsTable)
}

allGenAllArrPopStatsTable <- checkAllGenAllArrPopStatsTable(allGenAllArrTable)



if (!exists("finalGenAllArrPopStatsTable")){
  if (file.exists(paste0(resultsDir,"finalGenAllArrPopStatsTable.Rdata"))){
    load(paste0(resultsDir,"finalGenAllArrPopStatsTable.Rdata"))
  } else {
    finalGenAllArrPopStatsTable <- rbind(filter(allGenAllArrPopStatsTable,
                                                   Generation==max(allGenAllArrPopStatsTable$Generation)),
                                            filter(allGenAllArrPopStatsTable,
                                                   Generation==max(allGenAllArrPopStatsTable$Generation)-1,
                                                   LifeHistoryStage!="Embryo"))
    save(finalGenAllArrPopStatsTable,file=paste0(resultsDir,"finalGenAllArrPopStatsTable.Rdata"))
  }
}


rm(allGenAllArrPopStatsTable)
rm(finalGenAllArrPopStatsTable)

