# For extracting and analyzing full SAIsim model data, sourced from other scripts
#
# Example run: Rscript Extract_Full_dmel.R ./dmelN3/


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


# if (!exists("fullFinalGenStatTable")){
#   if (file.exists(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))){
#     load(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))
#   } else {
#     if (!exists("fullFinalGenN3StatTable")){
#       fullFinalGenN3StatTable <- extractFinalGenInvStatTable(N3dir,c(1e-3,1e-5),1e-4,c(10,100),1,
#                                                      c(1e-2,1e-5),1e3,2e4,1,1000)
#     }
#     if (!exists("partE10FinalGenN4StatTable")){
#       partE10FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-3,1e-5),1e-4,c(10),1,
#                                                         c(1e-2,1e-5),1e4,2e4,1,100)
#     }
#     if (!exists("partE100m5FinalGenN4StatTable")){
#       partE100m5FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-5),1e-4,c(100),1,
#                                                            c(1e-2,1e-5),1e4,2e4,1,100)
#     }
#     if (!exists("partE100m3c2FinalGenN4StatTable")){
#       partE100m3c2FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-3),1e-4,c(100),1,
#                                                              c(1e-2),1e4,2e4,1,100)
#     }
#     if (!exists("partE10FinalGenN2StatTable")){
#       partE10FinalGenN2StatTable <- extractFinalGenInvStatTable(N2dir,c(1e-3,1e-5),1e-4,c(10),1,
#                                                         c(1e-2,1e-5),1e2,2e4,1,1000)
#     }
#     fullFinalGenStatTable <- rbind(fullFinalGenN3StatTable,
#                                    partE10FinalGenN2StatTable,
#                                    partE10FinalGenN4StatTable,
#                                    partE100m5FinalGenN4StatTable,
#                                    partE100m3c2FinalGenN4StatTable)
#     save(fullFinalGenStatTable,file=paste0(resultsDir,"fullFinalGenStatTable.Rdata"))
#   }
# }


# if (!exists("finalGenDiffsTable")){
#   if (file.exists(paste0(resultsDir,"finalGenDiffsTable.Rdata"))){
#     load(paste0(resultsDir,"finalGenDiffsTable.Rdata"))
#   } else {
#     finalGenDiffsTable <- calcDiffsTable(fullFinalGenStatTable)
#     save(finalGenDiffsTable,file=paste0(resultsDir,"finalGenDiffsTable.Rdata"))
#   }
# }
# if (!exists("finalGenAvgDiffsTable")){
#   if (file.exists(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))){
#     load(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))
#   } else {
#     finalGenAvgDiffsTable <- calcAvgDiffsTable(fullFinalGenStatTable)
#     save(finalGenAvgDiffsTable,file=paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))
#   }
# }


checkAllGenInvStatTable <- function(){
  if (!exists("allGenInvStatTable")){
    if (file.exists(paste0(resultsDir,"allGenInvStatTable.Rdata"))){
      load(paste0(resultsDir,"allGenInvStatTable.Rdata"))
    } else {
      partialAllGenN3L1MbpRareMutStatTable <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                              c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                              c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1),
                                                              finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L1MbpRareMutStatTableNMC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                                 c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                 c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1),
                                                                 finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L1MbpRareMutStatTableNFC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                                 c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                 c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1),
                                                                 finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpRareMutStatTable <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                               c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                               c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1),
                                                               finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpRareMutStatTableNMC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                                  c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                  c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1),
                                                                  finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpRareMutStatTableNFC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                                  c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                  c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1),
                                                                  finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L1MbpStatTable <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                       c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                       c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1),
                                                       finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L1MbpStatTableNMC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                          c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                          c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1),
                                                          finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L1MbpStatTableNFC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                          c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                          c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1),
                                                          finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpStatTable <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                        c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                        c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1),
                                                        finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpStatTableNMC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                           c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                           c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1),
                                                           finishingFunc=mutInvStatTableAddFreq)
      partialAllGenN3L10KbpStatTableNFC <- parallelExtract(extractInvStatRow,list(resultsDir,
                                                           c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                           c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1),
                                                           finishingFunc=mutInvStatTableAddFreq)
      allGenInvStatTable <- rbind(partialAllGenN3L1MbpRareMutStatTable,
                                   partialAllGenN3L1MbpRareMutStatTableNMC,
                                   partialAllGenN3L1MbpRareMutStatTableNFC,
                                   partialAllGenN3L10KbpRareMutStatTable,
                                   partialAllGenN3L10KbpRareMutStatTableNMC,
                                   partialAllGenN3L10KbpRareMutStatTableNFC,
                                   partialAllGenN3L1MbpStatTable,
                                   partialAllGenN3L1MbpStatTableNMC,
                                   partialAllGenN3L1MbpStatTableNFC,
                                   partialAllGenN3L10KbpStatTable,
                                   partialAllGenN3L10KbpStatTableNMC,
                                   partialAllGenN3L10KbpStatTableNFC)
      save(allGenInvStatTable,file=paste0(resultsDir,"allGenInvStatTable.Rdata"))
    }
  }
  return(allGenInvStatTable)
}

allGenInvStatTable <- checkAllGenInvStatTable()

# extractInvStatTable(resultsDir,
#                     c(1.042e-2),c(1.042e-4,0),c(10,100),c(4.36e1),c(1.295e-1),
#                     c(1e3),c(2e4),c(1),c(0),c(0),1)


if (!exists("allGenDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenDiffsTable.Rdata"))
  } else {
    allGenDiffsTable <- parallelAnalysis(calcSpecParamDiffsTable,
                                         list(list(allGenInvStatTable)),
                                         getSimParamLs(allGenInvStatTable,includeGens = TRUE))
    save(allGenDiffsTable,file=paste0(resultsDir,"allGenDiffsTable.Rdata"))
  }
}
if (!exists("allGenAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))
  } else {
    allGenAvgDiffsTable <- parallelAnalysis(calcParamGenSpecAvgDiffsRow,
                                            list(list(allGenInvStatTable)),
                                            getSimParamLs(allGenInvStatTable,includeGens = TRUE))
    save(allGenAvgDiffsTable,file=paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))
  }
}

rm(allGenInvStatTable)
rm(allGenDiffsTable)
rm(allGenAvgDiffsTable)

# testallGenDiffsTable <- parallelStatCalc(calcSpecParamDiffsTable,filter(fullAllGenStatTable,
#                                                                         Generation %in% seq(0,20000,by=2000),
#                                                                         ReplicateNum %in% seq(0,10)))





checkAllGenInvFreqTable <- function(){
  if (!exists("allGenInvFreqTable")){
    if (file.exists(paste0(resultsDir,"allGenInvFreqTable.Rdata"))){
      load(paste0(resultsDir,"allGenInvFreqTable.Rdata"))
    } else {
      partialAllGenN3L1MbpRareMutInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                  c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                  c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      partialAllGenN3L1MbpRareMutInvFreqTableNMC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                     c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                     c(1e3),c(2e4),c(1),c(1),c(0),numReps))
      partialAllGenN3L1MbpRareMutInvFreqTableNFC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                     c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                     c(1e3),c(2e4),c(1),c(0),c(1),numReps))
      partialAllGenN3L10KbpRareMutInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                   c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                   c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      partialAllGenN3L10KbpRareMutInvFreqTableNMC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                      c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                      c(1e3),c(2e4),c(1),c(1),c(0),numReps))
      partialAllGenN3L10KbpRareMutInvFreqTableNFC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                      c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                      c(1e3),c(2e4),c(1),c(0),c(1),numReps))
      
      partialAllGenN3L1MbpInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                           c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                           c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      partialAllGenN3L1MbpInvFreqTableNMC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                              c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                              c(1e3),c(2e4),c(1),c(1),c(0),numReps))
      partialAllGenN3L1MbpInvFreqTableNFC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                              c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                              c(1e3),c(2e4),c(1),c(0),c(1),numReps))
      partialAllGenN3L10KbpInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                            c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                            c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      partialAllGenN3L10KbpInvFreqTableNMC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                               c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                               c(1e3),c(2e4),c(1),c(1),c(0),numReps))
      partialAllGenN3L10KbpInvFreqTableNFC <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                               c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                               c(1e3),c(2e4),c(1),c(0),c(1),numReps))
      
      partialAllGenN3L1MbpNoMutInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                c(0),c(1.042e-6,1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      partialAllGenN3L10KbpNoMutInvFreqTable <- parallelExtract(extractSpecParamInvSFSTable,list(resultsDir,
                                                                                                 c(0),c(1.042e-8,1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                 c(1e3),c(2e4),c(1),c(0),c(0),numReps))
      allGenInvFreqTable <- rbind(partialAllGenN3L1MbpRareMutInvFreqTable,
                                  partialAllGenN3L1MbpRareMutInvFreqTableNMC,
                                  partialAllGenN3L1MbpRareMutInvFreqTableNFC,
                                  partialAllGenN3L10KbpRareMutInvFreqTable,
                                  partialAllGenN3L10KbpRareMutInvFreqTableNMC,
                                  partialAllGenN3L10KbpRareMutInvFreqTableNFC,
                                  partialAllGenN3L1MbpInvFreqTable,
                                  partialAllGenN3L1MbpInvFreqTableNMC,
                                  partialAllGenN3L1MbpInvFreqTableNFC,
                                  partialAllGenN3L10KbpInvFreqTable,
                                  partialAllGenN3L10KbpInvFreqTableNMC,
                                  partialAllGenN3L10KbpInvFreqTableNFC,
                                  partialAllGenN3L1MbpNoMutInvFreqTable,
                                  partialAllGenN3L10KbpNoMutInvFreqTable)
      # allGenInvFreqTable <- extractInvSFSTable(N3dir,c(0,1e-3),c(1e-4),c(100),1,c(1e-2),1e3,2e4,1,1000)
      save(allGenInvFreqTable,file=paste0(resultsDir,"allGenInvFreqTable.Rdata"))
    }
  }
  return(allGenInvFreqTable)
}

allGenInvFreqTable <- checkAllGenInvFreqTable()

rm(allGenInvFreqTable)




checkAllGenPopStatTable <- function(){
  if (!exists("allGenPopStatTable")){
    if (file.exists(paste0(resultsDir,"allGenPopStatTable.Rdata"))){
      load(paste0(resultsDir,"allGenPopStatTable.Rdata"))
    } else {
      partialAllGenN3L1MbpRareMutPopStatTable <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                 c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                 c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutPopStatTableNMC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                    c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                    c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutPopStatTableNFC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                    c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                    c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutPopStatTable <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                  c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                  c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutPopStatTableNMC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                     c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                     c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutPopStatTableNFC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                                     c(1.042e-6),c(1.042e-8,1.042e-6,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                     c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L1MbpPopStatTable <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                          c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                          c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpPopStatTableNMC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                             c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                             c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpPopStatTableNFC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                             c(1.042e-2),c(1.042e-4,1.042e-2,0),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                             c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpPopStatTable <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                           c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                           c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpPopStatTableNMC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                              c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                              c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpPopStatTableNFC <- parallelExtract(extractPopStatRow,list(resultsDir,
                                                              c(1.042e-4),c(1.042e-6,1.042e-4,0),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                              c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      allGenPopStatTable <- rbind(partialAllGenN3L1MbpRareMutPopStatTable,
                                  partialAllGenN3L1MbpRareMutPopStatTableNMC,
                                  partialAllGenN3L1MbpRareMutPopStatTableNFC,
                                  partialAllGenN3L10KbpRareMutPopStatTable,
                                  partialAllGenN3L10KbpRareMutPopStatTableNMC,
                                  partialAllGenN3L10KbpRareMutPopStatTableNFC,
                                  partialAllGenN3L1MbpPopStatTable,
                                  partialAllGenN3L1MbpPopStatTableNMC,
                                  partialAllGenN3L1MbpPopStatTableNFC,
                                  partialAllGenN3L10KbpPopStatTable,
                                  partialAllGenN3L10KbpPopStatTableNMC,
                                  partialAllGenN3L10KbpPopStatTableNFC)
      allGenPopStatTable$LifeHistoryStage <- sub("'",'',sub("'",'',allGenPopStatTable[["LifeHistoryStage"]]))
      save(allGenPopStatTable,file=paste0(resultsDir,"allGenPopStatTable.Rdata"))
    }
  }
  return(allGenPopStatTable)
}

allGenPopStatTable <- checkAllGenPopStatTable()

# allGenPopStatTable$LifeHistoryStage <- sub("'",'',sub("'",'',allGenPopStatTable[["LifeHistoryStage"]]))


# allGenAvgPopStatTable <- averagePopStatsByRep(allGenPopStatTable)
if (!exists("allGenAvgPopStatTable")){
  if (file.exists(paste0(resultsDir,"allGenAvgPopStatTable.Rdata"))){
    load(paste0(resultsDir,"allGenAvgPopStatTable.Rdata"))
  } else {
    allGenAvgPopStatTable <- averagePopStatsByRep(allGenPopStatTable)
    save(allGenAvgPopStatTable,file=paste0(resultsDir,"allGenAvgPopStatTable.Rdata"))
  }
}

rm(allGenAvgPopStatTable)




checkAllGenAllInvTable <- function(){
  if (!exists("allGenAllInvTable")){
    if (file.exists(paste0(resultsDir,"allGenAllInvTable.Rdata"))){
      load(paste0(resultsDir,"allGenAllInvTable.Rdata"))
    } else {
      partialAllGenN3L1MbpRareMutAllInvTable <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                     c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                     c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutAllInvTableNMC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                        c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                        c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutAllInvTableNFC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                        c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                        c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllInvTable <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                      c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                      c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllInvTableNMC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                         c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                         c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllInvTableNFC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                         c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                         c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L1MbpAllInvTable <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                              c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                              c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpAllInvTableNMC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                 c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                 c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpAllInvTableNFC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                 c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                 c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpAllInvTable <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                               c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                               c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpAllInvTableNMC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                  c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                  c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpAllInvTableNFC <- parallelExtract(extractSpecParamAllInvs,list(resultsDir,
                                                                                  c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                  c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      allGenAllInvTable <- rbind(partialAllGenN3L1MbpRareMutAllInvTable,
                                  partialAllGenN3L1MbpRareMutAllInvTableNMC,
                                  partialAllGenN3L1MbpRareMutAllInvTableNFC,
                                  partialAllGenN3L10KbpRareMutAllInvTable,
                                  partialAllGenN3L10KbpRareMutAllInvTableNMC,
                                  partialAllGenN3L10KbpRareMutAllInvTableNFC,
                                  partialAllGenN3L1MbpAllInvTable,
                                  partialAllGenN3L1MbpAllInvTableNMC,
                                  partialAllGenN3L1MbpAllInvTableNFC,
                                  partialAllGenN3L10KbpAllInvTable,
                                  partialAllGenN3L10KbpAllInvTableNMC,
                                  partialAllGenN3L10KbpAllInvTableNFC)
      save(allGenAllInvTable,file=paste0(resultsDir,"allGenAllInvTable.Rdata"))
    }
  }
  return(allGenAllInvTable)
}

allGenAllInvTable <- checkAllGenAllInvTable()



checkAllGenAllInvPopAvgDiffsTable <- function(allGenAllInvTable,allGenPopStatTable){
  if (!exists("allGenAllInvPopAvgDiffsTable")){
    if (file.exists(paste0(resultsDir,"allGenAllInvPopAvgDiffsTable.Rdata"))){
      load(paste0(resultsDir,"allGenAllInvPopAvgDiffsTable.Rdata"))
    } else {
      # partialAllGenN3L1MbpRareMutAllInvPopDiffsTable <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                    list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                    list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                         c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpRareMutAllInvPopDiffsTableNMC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                       list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                       list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                            c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpRareMutAllInvPopDiffsTableNFC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                       list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                       list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                            c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllInvPopDiffsTable <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                     list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                     list(c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                          c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllInvPopDiffsTableNMC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                        list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                        list(c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                             c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpRareMutAllInvPopDiffsTableNFC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                        list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                        list(c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                             c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      # partialAllGenN3L1MbpAllInvPopDiffsTable <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                             list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                             list(c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                  c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpAllInvPopDiffsTableNMC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                list(c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                     c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L1MbpAllInvPopDiffsTableNFC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                list(c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
      #                                                                     c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      # partialAllGenN3L10KbpAllInvPopDiffsTable <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                              list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                              list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                   c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpAllInvPopDiffsTableNMC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                 list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                 list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                      c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      # partialAllGenN3L10KbpAllInvPopDiffsTableNFC <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,
      #                                                                 list(list(allGenAllInvTable),list(allGenPopStatTable)),
      #                                                                 list(c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
      #                                                                      c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      # allGenAllInvPopAvgDiffsTable <- rbind(partialAllGenN3L1MbpRareMutAllInvPopDiffsTable,
      #                            partialAllGenN3L1MbpRareMutAllInvPopDiffsTableNMC,
      #                            partialAllGenN3L1MbpRareMutAllInvPopDiffsTableNFC,
      #                            partialAllGenN3L10KbpRareMutAllInvPopDiffsTable,
      #                            partialAllGenN3L10KbpRareMutAllInvPopDiffsTableNMC,
      #                            partialAllGenN3L10KbpRareMutAllInvPopDiffsTableNFC,
      #                            partialAllGenN3L1MbpAllInvPopDiffsTable,
      #                            partialAllGenN3L1MbpAllInvPopDiffsTableNMC,
      #                            partialAllGenN3L1MbpAllInvPopDiffsTableNFC,
      #                            partialAllGenN3L10KbpAllInvPopDiffsTable,
      #                            partialAllGenN3L10KbpAllInvPopDiffsTableNMC,
      #                            partialAllGenN3L10KbpAllInvPopDiffsTableNFC)
      filteredTableLs <- varCoreFilteredTableLs(list(allGenAllInvTable,allGenPopStatTable),maxSplitInd=11)
      allGenAllInvPopAvgDiffsTable <- parallelAnalysis(calcAllInvPopAvgDiffs,filteredTableLs,list())
      # allGenAllInvPopAvgDiffsTable <- parallelAnalysis(calcSpecParamAllInvPopAvgDiffs,filteredTableLs,list())
      save(allGenAllInvPopAvgDiffsTable,file=paste0(resultsDir,"allGenAllInvPopAvgDiffsTable.Rdata"))
    }
  }
  return(allGenAllInvPopAvgDiffsTable)
}


allGenAllInvPopAvgDiffsTable <- checkAllGenAllInvPopAvgDiffsTable(allGenAllInvTable,allGenPopStatTable)


if (!exists("finalGenAllInvPopAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"finalGenAllInvPopAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"finalGenAllInvPopAvgDiffsTable.Rdata"))
  } else {
    finalGenAllInvPopAvgDiffsTable <- rbind(filter(allGenAllInvPopAvgDiffsTable,
                                                   Generation==max(allGenAllInvPopAvgDiffsTable$Generation)),
                                            filter(allGenAllInvPopAvgDiffsTable,
                                                   Generation==max(allGenAllInvPopAvgDiffsTable$Generation)-1,
                                                   LifeHistoryStage=="Adult"))
    save(finalGenAllInvPopAvgDiffsTable,file=paste0(resultsDir,"finalGenAllInvPopAvgDiffsTable.Rdata"))
  }
}

rm(allGenAllInvPopAvgDiffsTable)
rm(finalGenAllInvPopAvgDiffsTable)


# if (!exists("i4e100c2N3InvFreqTable")){
#   if (file.exists(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))){
#     load(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))
#   } else {
#     i4e100c2N3InvFreqTable <- extractInvSFSTable(N3dir,c(0,1e-3),c(1e-4),c(100),1,
#                                             c(1e-2),1e3,2e4,1,1000)
#     save(i4e100c2N3InvFreqTable,file=paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))
#   }
# }
# 
# 
# if (!exists("m3e100c2N3PopStatTable")){
#   if (file.exists(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))){
#     load(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
#   } else {
#     m3e100c2N3PopStatTable <- extractPopStatTable(N3dir,c(1e-3),c(0,1e-4),c(100),1,
#                                             c(1e-2),1e3,2e4,1,1000)
#     save(m3e100c2N3PopStatTable,file=paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
#   }
# }

# repAvgm3e100c2N3PopStatTable <- averagePopStatsByRep(m3e100c2N3PopStatTable)





checkAllGenAllArrTable <- function(){
  if (!exists("allGenAllArrTable")){
    if (file.exists(paste0(resultsDir,"allGenAllArrTable.Rdata"))){
      load(paste0(resultsDir,"allGenAllArrTable.Rdata"))
    } else {
      partialAllGenN3L1MbpRareMutAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                             c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                             c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                                c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpRareMutAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                                c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                                c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                              c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                              c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                                 c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                 c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpRareMutAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                                 c(1.042e-6),c(1.042e-8,1.042e-6),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                                 c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L1MbpAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                      c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                      c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                         c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                         c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L1MbpAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                         c(1.042e-2),c(1.042e-4,1.042e-2),c(2,10,100),c(4.36e1),c(1.295e-1),
                                                                                         c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      partialAllGenN3L10KbpAllArrTable <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                       c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                       c(1e3),c(2e4),c(1),c(0),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpAllArrTableNMC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                          c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                          c(1e3),c(2e4),c(1),c(1),c(0),seq(numReps)-1))
      partialAllGenN3L10KbpAllArrTableNFC <- parallelExtract(extractSpecParamAllArrs,list(resultsDir,
                                                                                          c(1.042e-4),c(1.042e-6,1.042e-4),c(2,10,100),c(4.36e-1),c(1.295e-1),
                                                                                          c(1e3),c(2e4),c(1),c(0),c(1),seq(numReps)-1))
      allGenAllArrTable <- rbind(partialAllGenN3L1MbpRareMutAllArrTable,
                                 partialAllGenN3L1MbpRareMutAllArrTableNMC,
                                 partialAllGenN3L1MbpRareMutAllArrTableNFC,
                                 partialAllGenN3L10KbpRareMutAllArrTable,
                                 partialAllGenN3L10KbpRareMutAllArrTableNMC,
                                 partialAllGenN3L10KbpRareMutAllArrTableNFC,
                                 partialAllGenN3L1MbpAllArrTable,
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
      filteredTableLs <- varCoreFilteredTableLs(list(allGenAllArrTable),maxSplitInd=11)
      allGenAllArrPopStatsTable <- parallelAnalysis(calcAllArrPopStats,filteredTableLs,list())
      save(allGenAllArrPopStatsTable,file=paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))
    }
  }
  return(allGenAllArrPopStatsTable)
}

allGenAllArrPopStatsTable <- checkAllGenAllArrPopStatsTable(allGenAllArrTable,allGenPopStatTable)


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