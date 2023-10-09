# For extracting and analyzing full SAIsim model data, sourced from other scripts
#

library(tidyverse)
library(parallel)

# Global parameters
coreProportion = 0.6
minCores = 2

readLastLineToTable <- function(file) {
  read_tsv(pipe(paste('(head -n 1; tail -n 1) <', file)), show_col_types = FALSE)
}

extractNonFixedInvs <- function(invFreqs,alleleSize){
  nonZeroFreqs <- invFreqs[-1:-3][which(invFreqs[-1:-3]!=0)]
  # print(nonZeroFreqs)
  # Allow for populations with only Std karyotype
  if (colnames(nonZeroFreqs)[1] == "Std"){
    nonFixedInvFreqs <- nonZeroFreqs[-1][which(nonZeroFreqs[-1]!=alleleSize)]
    # print(nonFixedInvFreqs)
    # print(ncol(nonFixedInvFreqs)+1)
    if (ncol(nonFixedInvFreqs) > 0){
      nonFixedFreqs <- cbind(nonZeroFreqs[1],nonFixedInvFreqs)
      # print(nonFixedFreqs)
      maxInvInd <- which.max(nonFixedFreqs)
      maxInvID <- colnames(nonFixedFreqs)[maxInvInd]
      otherIDs <- colnames(nonFixedFreqs[-maxInvInd])
      # print(otherIDs)
      # if (maxInvID != "Std"){
      #   maxInvID <- sub('.', '', maxInvID)
      # }
      # if (otherIDs[1] == "Std"){
      #   fixedIDs <- c(otherIDs[1],sub('.', '', otherIDs[-1]))
      #   otherIDs <- fixedIDs
      # } else {
      #   fixedIDs <- sub('.', '', otherIDs)
      #   otherIDs <- fixedIDs
      # }
      orderedIDs <- sub('I','',c(maxInvID,otherIDs))
      return(orderedIDs)
    } else{
      if (nonZeroFreqs$Std != alleleSize){
        print("WARNING - how did this happen? Std arrangement the only non-fixed, but not at count = alleleSize")
      }
      return("Std")
    }
  } else{
    # Should be invariant that there are at least two non-fixed inversions if count of Std = 0
    nonFixedInvFreqs <- nonZeroFreqs[which(nonZeroFreqs!=alleleSize)]
    maxInvInd <- which.max(nonFixedInvFreqs)
    maxInvID <- colnames(nonFixedInvFreqs)[maxInvInd]
    otherIDs <- colnames(nonFixedInvFreqs[-maxInvInd])
    if ("Std" %in% colnames(nonFixedInvFreqs)){
      print("WARNING - Std arrangement occurred after the second column, unaccounted")
    }
    orderedIDs <- sub('.', '', c(maxInvID,otherIDs))
    return(orderedIDs)
  }
}

extractFinalGenInv <- function(simDir,expFinalGen,alleleSize){
  invFreqFile = paste0(simDir,"InvFreqs.txt")
  finalGenFreqs = readLastLineToTable(invFreqFile)
  # print(finalGenFreqs)
  if (finalGenFreqs[1,1] != expFinalGen){
    print(past0("WARNING - final generation in ",invFreqFile," is ",
                finalGenFreqs[1,1]," rather than the expected ",expFinalGen))
  }
  return(extractNonFixedInvs(finalGenFreqs,alleleSize))
}



getAvgStatsFinalGenInvs <- function(simDir,IDs){
  avgStatsInvTable <- readLastLineToTable(paste0(simDir,"Inv",IDs[1],".txt"))
  for (ID in IDs[-1]){
    avgStatsInvTable <- rbind(avgStatsInvTable,readLastLineToTable(paste0(simDir,"Inv",ID,".txt")))
  }
  avgStatsInvTable <- mutate(avgStatsInvTable,InvID=IDs,.before=2)
  return(avgStatsInvTable)
}

avgNonMaxInvs <- function(invStats){
  avgStats <- invStats
  if (nrow(avgStats) > 2){
    otherInvs <- avgStats[-1,]
    gen = otherInvs$Generation[1]
    lhs = otherInvs$LifeHistoryStage[1]
    if (any(otherInvs$Generation != gen)){
      print(paste0("Warning: not all generations are the same: ",otherInvs$Generation))
    }
    if (any(otherInvs$LifeHistoryStage != lhs)){
      print(paste0("Warning: not all life history stages are the same: ",otherInvs$LifeHistoryStage))
    }
    countData = subset(otherInvs, select=c(Count,Heterozygote,Homozygote,
                                           FemaleCount,FemaleHeterozygote,FemaleHomozygote,
                                           MaleCount,MaleHeterozygote,MaleHomozygote))
    countTotals = colSums(countData)
    weights = otherInvs$Count/sum(otherInvs$Count)
    # print(otherInvs)
    # print(weights)
    avgOthers <- colSums(subset(otherInvs, select=-c(InvID,Generation,LifeHistoryStage,
                                                     Count,Heterozygote,Homozygote,
                                                     FemaleCount,FemaleHeterozygote,FemaleHomozygote,
                                                     MaleCount,MaleHeterozygote,MaleHomozygote))*weights)
    avgOtherRow <- mutate(cbind(as_tibble_row(countTotals),as_tibble_row(avgOthers)),
                          Generation=gen,LifeHistoryStage=lhs,InvID="Other",.before=1)
    # print(avgOtherRow)
    # print(avgStats[1,])
    avgStats <- rbind(avgStats[1,],avgOtherRow)
  }
  if (nrow(avgStats) > 1){
    avgStats$InvID <- c("Max","Other")
  } else {
    avgStats$InvID <- c("Max")
  }
  return(avgStats)
}

getSimDir <- function(resultsDir,mutRate,mutRateInv,encNum,
                      crossRate,convRate,popSize,numGens,
                      wrightFisher,noMaleCost,noFemaleCost,
                      repNum){
  m <- formatC(mutRate, digits = 3, width = 4, format = "e", flag = "0")
  tempSimDir <- paste0('mut',m)
  i <- formatC(mutRateInv, digits = 3, width = 4, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'inv',i)
  e <- formatC(encNum, width = 3, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'enc',e)
  r <- formatC(crossRate, digits = 3, width = 4, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'cro',r)
  c <- formatC(convRate, digits = 3, width = 4, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'con',c)
  N <- formatC(popSize, digits = 1, width = 2, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'pop',N)
  g <- formatC(numGens, digits = 1, width = 2, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'gen',g)
  wf <- formatC(wrightFisher, digits = 1, width = 1, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'wf',wf)
  nmc <- formatC(noMaleCost, digits = 1, width = 1, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'nmc',nmc)
  nfc <- formatC(noFemaleCost, digits = 1, width = 1, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'nfc',nfc)
  n <- formatC(repNum, width = 4, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'rep',n)
  return(paste0(resultsDir,tempSimDir,'/'))
}

extractFinalGenInvStatRow <- function(resultsDir,mutRate,invRate,encNum,
                                  recRate,conRate,popSize,numGens,
                                  wrightFisher,noMaleCost,noFemaleCost,
                                  repNum) {
  alleleSize <- 2*popSize
  simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                     recRate,conRate,popSize,numGens,
                     wrightFisher,noMaleCost,noFemaleCost,
                     repNum)
  print(simDir)
  finalGenInvData <- extractFinalGenInv(simDir,numGens,alleleSize)
  statData <- getAvgStatsFinalGenInvs(simDir,finalGenInvData)
  maxAltStatData <- avgNonMaxInvs(statData)
  statRow <- mutate(maxAltStatData,
                    MutRate=mutRate,
                    InvRate=invRate,
                    EncounterNum=encNum,
                    CrossoverMapLength=recRate,
                    ConversionRate=conRate,
                    PopSize=popSize,
                    MaxGeneration=numGens,
                    WrightFisherSurvival=as.logical(wrightFisher),
                    NoMaleCost=as.logical(noMaleCost),
                    NoFemaleCost=as.logical(noFemaleCost),
                    ReplicateNum=repNum,
                    .before=1)
  return(statRow)
}

mutInvStatTableAddFreq <- function(invStatTable){
  return(mutate(invStatTable,Frequency=invStatTable$Count/(2*invStatTable$PopSize),.before=15))
}


extractFinalGenInvStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,nmcs,nfcs,repN){
  invStatTable <- tibble()
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  for (noMaleCost in nmcs){
                    for (noFemaleCost in nfcs){
                      # popSize <- 10^size
                      # alleleSize <- 2*popSize
                      for (repNum in seq(repN)-1){
                        statRow <- extractFinalGenInvStatRow(resultsDir,mutRate,invRate,encNum,
                                                             recRate,conRate,popSize,numGens,
                                                             wrightFisher,noMaleCost,noFemaleCost,
                                                             repNum)
                        invStatTable <- rbind(invStatTable,statRow)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  invStatTable <- mutInvStatTableAddFreq(invStatTable)
  return(invStatTable)
}

extractAvgStatsAllGenInv <- function(simDir,alleleSize,lhs="Zygote",gens=seq(0,20000,by=200)){
  invFreqFile = paste0(simDir,"InvFreqs.txt")
  invFreqs = read_tsv(invFreqFile, show_col_types = FALSE)
  # print(invFreqs)
  invFreqs <- filter(invFreqs,Generation %in% gens,LifeHistoryStage==lhs)
  # print(invFreqs[1,])
  nonPolyCols <- c()
  for (colName in colnames(invFreqs)[-1:-3]) {
    if (all(invFreqs[colName] == 0)) {
      nonPolyCols <- c(nonPolyCols,colName)
    }
  }
  # print(nonPolyCols)
  if (!is.null(nonPolyCols)){
    invFreqs <- subset(invFreqs,select=-which(colnames(invFreqs) %in% nonPolyCols))
  }
  invIDs <- sub('I','',colnames(invFreqs[1,-1:-3]))
  # print(invIDs)
  # print(invFreqs[1,])
  invStatTables <- list()
  for (invID in invIDs) {
    specInvStatTable <- read_tsv(paste0(simDir,"Inv",invID,".txt"),
                                 show_col_types = FALSE)
    # names(specInvStatTable) <- c('Generation','LifeHistoryStage',
    #                              'Count','Homozygote','Heterozygote',
    #                              'FemaleCount','FemaleHomozygote','FemaleHeterozygote',
    #                              'MaleCount','MaleHomozygote','MaleHeterozygote',
    #                              'AvgNumMut','VarNumMut','ModeNumMut',
    #                              'AvgSurEff','VarSurEff','ModeSurEff',
    #                              'AvgRepEff','VarRepEff','ModeRepEff')
    invStatTables <- c(invStatTables,
                       list(specInvStatTable))
  }
  names(invStatTables) <- invIDs
  # print(invStatTables)
  avgInvStatsAllGens <- tibble()
  for (gen in invFreqs$Generation){
    # print(filter(invFreqs,
    #        Generation == gen,
    #        LifeHistoryStage==lhs))
    nonFixedInvs <- extractNonFixedInvs(filter(invFreqs,
                                               Generation == gen,
                                               LifeHistoryStage==lhs),
                                        alleleSize)
    invStats <- tibble()
    for (ID in nonFixedInvs){
      # print(ID)
      # print(filter(invStatTables[[ID]],Generation == gen))
      invStats <- rbind(invStats,
                        filter(invStatTables[[ID]],
                               Generation == gen,
                               LifeHistoryStage==lhs))
    }
    invStats <- mutate(invStats,InvID=nonFixedInvs,.before=3)
    avgStatsThisGen <- avgNonMaxInvs(invStats)
    avgInvStatsAllGens <- rbind(avgInvStatsAllGens,avgStatsThisGen)
  }
  return(avgInvStatsAllGens)
}

extractInvStatRow <- function(resultsDir,mutRate,invRate,encNum,
                              recRate,conRate,popSize,numGens,
                              wrightFisher,noMaleCost,noFemaleCost,
                              repNum) {
  alleleSize <- 2*popSize
  simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                     recRate,conRate,popSize,numGens,
                     wrightFisher,noMaleCost,noFemaleCost,
                     repNum)
  print(simDir)
  statData <- extractAvgStatsAllGenInv(simDir,alleleSize)
  statRow <- mutate(statData,
                    MutRate=mutRate,
                    InvRate=invRate,
                    EncounterNum=encNum,
                    CrossoverMapLength=recRate,
                    ConversionRate=conRate,
                    PopSize=popSize,
                    MaxGeneration=numGens,
                    WrightFisherSurvival=as.logical(wrightFisher),
                    NoMaleCost=as.logical(noMaleCost),
                    NoFemaleCost=as.logical(noFemaleCost),
                    ReplicateNum=repNum,
                    .before=1)
  return(statRow)
}

extractInvStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,nmcs,nfcs,repN){
  invStatTable <- tibble()
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  for (noMaleCost in nmcs){
                    for (noFemaleCost in nfcs){
                      # alleleSize <- 2*popSize
                      for (repNum in seq(repN)-1){
                        statRow <- extractInvStatRow(resultsDir,mutRate,invRate,encNum,
                                                     recRate,conRate,popSize,numGens,
                                                     wrightFisher,noMaleCost,noFemaleCost,
                                                     repNum)
                        invStatTable <- rbind(invStatTable,statRow)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  invStatTable <- mutInvStatTableAddFreq(invStatTable)
  return(invStatTable)
}

calcSpecParamDiffsTable <- function(statTable, m, i, e, r, c, N, g, w, nmc, nfc, gen){
  # print(c(m, i, e, r, c, N, g, w, nmc, nfc, gen))
  # print(statTable)
  # e <- as.integer(e)
  # w <- as.logical(w)
  # nmc <- as.logical(nmc)
  # nfc <- as.logical(nfc)
  specParamDiffsTable <- tibble()
  diffColNames <- paste0('AbsDiff',colnames(select(statTable,
                                                   !c(MutRate:MaleHeterozygote))))
  specificParamData <- filter(statTable,
                              MutRate == m,
                              InvRate == i,
                              EncounterNum == e,
                              CrossoverMapLength == r,
                              ConversionRate == c,
                              PopSize == N,
                              MaxGeneration == g,
                              WrightFisherSurvival == w,
                              NoMaleCost == nmc,
                              NoFemaleCost == nfc,
                              Generation == gen)
  for (num in unique(specificParamData$ReplicateNum)){
    specificReplicateGen <- filter(specificParamData,
                                   ReplicateNum == num)
    # print(specificReplicateGen)
    specificParams <- select(specificReplicateGen[1,],
                             c(MutRate:LifeHistoryStage))
    # print(specificParams)
    countFreqData <- select(specificReplicateGen[1,],
                            c(Count,Frequency))
    colnames(countFreqData) <- c("MaxArrangementCount","MaxArrangementFrequency")
    specificReplicateData <- select(specificReplicateGen,
                                    !c(MutRate:MaleHeterozygote))
    diffs <- tibble()
    if (nrow(specificReplicateData) == 1){
      # alleleNum <- specificReplicateData$Count
      diffs <- abs(specificReplicateData - specificReplicateData)
      # diffs$Count <- alleleNum
      colnames(diffs) <- diffColNames
    } else if (nrow(specificReplicateData) == 2){
      diffs <- abs(specificReplicateData[1,] - specificReplicateData[2,])
      colnames(diffs) <- diffColNames
    } else if (nrow(specificReplicateData) == 0){
      print(paste0("WARNING - NO entries for a specific sim replicate in given Max/Alt inversion data table"))
    } else {
      print(paste0("WARNING - >2 entries for a specific sim replicate in given Max/Alt inversion data table"))
      print(specificReplicateData)
    }
    # print(countFreqData)
    # print(diffs)
    diffsRow <- cbind(specificParams,countFreqData,diffs)
    specParamDiffsTable <- rbind(specParamDiffsTable,diffsRow)
  }
  # print(specParamDiffsTable)
  return(specParamDiffsTable)
}

calcDiffsTable <- function(statTable){
  # diffColNames <- paste0('AbsDiff',colnames(select(statTable,
  #                                                     !c(MutRate:Frequency))))
  diffsTable <- tibble()
  for (m in unique(statTable$MutRate)){
    for (i in unique(statTable$InvRate)){
      for (e in unique(statTable$EncounterNum)){
        for (r in unique(statTable$CrossoverMapLength)){
          for (c in unique(statTable$ConversionRate)){
            for (N in unique(statTable$PopSize)){
              for (g in unique(statTable$MaxGeneration)){
                for (w in unique(statTable$WrightFisherSurvival)){
                  for (nmc in unique(statTable$NoMaleCost)){
                    for (nfc in unique(statTable$NoFemaleCost)){
                      for (gen in unique(statTable$Generation)){
                        # specificParamData <- filter(statTable,
                        #                             MutRate == m,
                        #                             InvRate == i,
                        #                             EncounterNum == e,
                        #                             CrossoverMapLength == r,
                        #                             ConversionRate == c,
                        #                             PopSize == N,
                        #                             MaxGeneration == g,
                        #                             WrightFisherSurvival == w,
                        #                             NoMaleCost == nmc,
                        #                             NoFemaleCost == nfc,
                        #                             Generation == gen)
                        # for (num in unique(specificParamData$ReplicateNum)){
                        #   specificReplicateGen <- filter(specificParamData,
                        #                                  ReplicateNum == num)
                        #   specificParams <- select(specificReplicateGen[1,],
                        #                            c(MutRate:Generation))
                        #   countFreqData <- select(specificReplicateGen[1,],
                        #                           c(Count,Frequency))
                        #   colnames(countFreqData) <- c("MaxArrangementCount","MaxArrangementFrequency")
                        #   specificReplicateData <- select(specificReplicateGen,
                        #                                   !c(MutRate:Frequency))
                        #   diffs <- tibble()
                        #   if (nrow(specificReplicateData) == 1){
                        #     # alleleNum <- specificReplicateData$Count
                        #     diffs <- abs(specificReplicateData - specificReplicateData)
                        #     # diffs$Count <- alleleNum
                        #     colnames(diffs) <- diffColNames
                        #   } else if (nrow(specificReplicateData) == 2){
                        #     diffs <- abs(specificReplicateData[1,] - specificReplicateData[2,])
                        #     colnames(diffs) <- diffColNames
                        #   } else if (nrow(specificReplicateData) == 0){
                        #     print(paste0("WARNING - NO entries for a specific sim replicate in given Max/Alt inversion data table"))
                        #   } else {
                        #     print(paste0("WARNING - >2 entries for a specific sim replicate in given Max/Alt inversion data table"))
                        #     print(specificReplicateData)
                        #   }
                        #   # print(countFreqData)
                        #   # print(diffs)
                        #   diffsRow <- cbind(specificParams,countFreqData,diffs)
                        #   diffsTable <- rbind(diffsTable,diffsRow)
                        #   print(specificParams)
                        # }
                        specParamDiffsTable <- calcSpecParamDiffsTable(statTable, m, i, e, r, c, N, g, w, nmc, nfc, gen)
                        diffsTable <- rbind(diffsTable,specParamDiffsTable)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(diffsTable)
}

calcParamGenSpecAvgDiffsRow <- function(statTable, m, i, e, r, c, N, g, w, nmc, nfc, gen){
  diffColNames <- paste0('AvgAbsDiff',colnames(select(statTable[1,],
                                                      !c(MutRate:MaleHeterozygote))))
  specificParamGenData <- filter(statTable,
                                 MutRate == m,
                                 InvRate == i,
                                 EncounterNum == e,
                                 CrossoverMapLength == r,
                                 ConversionRate == c,
                                 PopSize == N,
                                 MaxGeneration == g,
                                 WrightFisherSurvival == w,
                                 NoMaleCost == nmc,
                                 NoFemaleCost == nfc,
                                 Generation == gen)
  if (nrow(specificParamGenData) > 0){
    specificParams <- select(specificParamGenData[1,],
                             c(MutRate:NoFemaleCost,Generation:LifeHistoryStage))
    specificParamDiffs <- tibble()
    numMonomorphic <- 0
    numPolymorphic <- 0
    for (num in unique(specificParamGenData$ReplicateNum)){
      specificReplicate <- filter(specificParamGenData,
                                  ReplicateNum == num)
      specificReplicateData <- select(specificReplicate,
                                      !c(MutRate:MaleHeterozygote))
      countFreqData <- select(specificReplicate[1,],
                              c(Count,Frequency))
      colnames(countFreqData) <- c("AvgMaxArrangementCount","AvgMaxArrangementFrequency")
      if (nrow(specificReplicateData) == 1){
        numMonomorphic <- numMonomorphic + 1
        diffs <- abs(specificReplicateData - specificReplicateData)
        colnames(diffs) <- diffColNames
        diffsWCountFreq <- cbind(countFreqData,diffs)
        specificParamDiffs <- rbind(specificParamDiffs,diffsWCountFreq)
      } else if (nrow(specificReplicateData) == 2){
        numPolymorphic <- numPolymorphic + 1
        diffs <- abs(specificReplicateData[1,] - specificReplicateData[2,])
        colnames(diffs) <- diffColNames
        diffsWCountFreq <- cbind(countFreqData,diffs)
        specificParamDiffs <- rbind(specificParamDiffs,diffsWCountFreq)
      } else if (nrow(specificReplicateData) == 0){
        print(paste0("WARNING - NO entries for a specific sim replicate in given Max/Alt inversion data table"))
      } else {
        print(paste0("WARNING - >2 entries for a specific sim replicate in given Max/Alt inversion data table"))
        print(specificReplicateData)
      }
    }
    totalReps <- sum(numMonomorphic,numPolymorphic)
    specificParamAvgDiffs <- colSums(specificParamDiffs)/totalReps
    avgDiffsRow <- cbind(specificParams,
                         NumMonomorphic=numMonomorphic,
                         NumPolymorphic=numPolymorphic,
                         NumReplicates=totalReps,
                         as_tibble_row(specificParamAvgDiffs))
    print(specificParams)
    return(avgDiffsRow)
  } else{
    return(NULL)
  }
}

calcAvgDiffsTable <- function(statTable){
  # diffColNames <- paste0('AvgAbsDiff',colnames(select(statTable[1,],
  #                                                     !c(MutRate:Frequency))))
  avgDiffsTable <- tibble()
  for (m in unique(statTable$MutRate)){
    for (i in unique(statTable$InvRate)){
      for (e in unique(statTable$EncounterNum)){
        for (r in unique(statTable$CrossoverMapLength)){
          for (c in unique(statTable$ConversionRate)){
            for (N in unique(statTable$PopSize)){
              for (g in unique(statTable$MaxGeneration)){
                for (w in unique(statTable$WrightFisherSurvival)){
                  for (nmc in unique(statTable$NoMaleCost)){
                    for (nfc in unique(statTable$NoFemaleCost)){
                      for (gen in unique(statTable$Generation)){
                        # specificParamGenData <- filter(statTable,
                        #                                MutRate == m,
                        #                                InvRate == i,
                        #                                EncounterNum == e,
                        #                                CrossoverMapLength == r,
                        #                                ConversionRate == c,
                        #                                PopSize == N,
                        #                                MaxGeneration == g,
                        #                                WrightFisherSurvival == w,
                        #                                NoMaleCost == nmc,
                        #                                NoFemaleCost == nfc,
                        #                                Generation == gen)
                        # if (nrow(specificParamGenData) > 0){
                        #   specificParams <- select(specificParamGenData[1,],
                        #                            c(MutRate:WrightFisherSurvival,Generation))
                        #   specificParamDiffs <- tibble()
                        #   numMonomorphic <- 0
                        #   numPolymorphic <- 0
                        #   for (num in unique(specificParamGenData$ReplicateNum)){
                        #     specificReplicate <- filter(specificParamGenData,
                        # #                                 ReplicateNum == num)
                        #     specificReplicateData <- select(specificReplicate,
                        #                                     !(MutRate:Frequency))
                        #     countFreqData <- select(specificReplicate[1,],
                        #                             c(Count,Frequency))
                        #     colnames(countFreqData) <- c("AvgMaxArrangementCount","AvgMaxArrangementFrequency")
                        #     if (nrow(specificReplicateData) == 1){
                        #       numMonomorphic <- numMonomorphic + 1
                        #       diffs <- abs(specificReplicateData - specificReplicateData)
                        #       colnames(diffs) <- diffColNames
                        #       diffsWCountFreq <- cbind(countFreqData,diffs)
                        #       specificParamDiffs <- rbind(specificParamDiffs,diffsWCountFreq)
                        #     } else if (nrow(specificReplicateData) == 2){
                        #       numPolymorphic <- numPolymorphic + 1
                        #       diffs <- abs(specificReplicateData[1,] - specificReplicateData[2,])
                        #       colnames(diffs) <- diffColNames
                        #       diffsWCountFreq <- cbind(countFreqData,diffs)
                        #       specificParamDiffs <- rbind(specificParamDiffs,diffsWCountFreq)
                        #     } else if (nrow(specificReplicateData) == 0){
                        #       print(paste0("WARNING - NO entries for a specific sim replicate in given Max/Alt inversion data table"))
                        #     } else {
                        #       print(paste0("WARNING - >2 entries for a specific sim replicate in given Max/Alt inversion data table"))
                        #       print(specificReplicateData)
                        #     }
                        #   }
                        #   totalReps <- sum(numMonomorphic,numPolymorphic)
                        #   specificParamAvgDiffs <- colSums(specificParamDiffs)/totalReps
                        #   avgDiffsRow <- cbind(specificParams,
                        #                        NumMonomorphic=numMonomorphic,
                        #                        NumPolymorphic=numPolymorphic,
                        #                        NumReplicates=totalReps,
                        #                        as_tibble_row(specificParamAvgDiffs))
                        #   avgDiffsTable <- rbind(avgDiffsTable,avgDiffsRow)
                        #   print(specificParams)
                        # }
                        avgDiffsRow <- calcParamGenSpecAvgDiffsRow(statTable,m,i,e,r,c,N,g,w,nmc,nfc,gen)
                        # if (!is.null(avgDiffsRow)) {
                        #   avgDiffsTable <- rbind(avgDiffsTable,avgDiffsRow)
                        # }
                        avgDiffsTable <- rbind(avgDiffsTable,avgDiffsRow)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(avgDiffsTable)
}


#
# Generating inversion SFS data for neutral inversion dynamic comparisons
#

extractSpecParamInvSFSTable <- function(resultsDir,mutRate,invRate,encNum,
                                        recRate,conRate,popSize,numGens,
                                        wrightFisher,noMaleCost,noFemaleCost,
                                        repN) {
  `%!in%` <- Negate(`%in%`)
  specificParamFreqCounts <- tibble()
  invFreqTables <- list()
  for (repNum in seq(repN)-1) {
    simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                       recRate,conRate,popSize,numGens,
                       wrightFisher,noMaleCost,noFemaleCost,
                       repNum)
    # print(paste0(simDir,"InvFreqs.txt"))
    # print("\n")
    invFreqTables <- c(invFreqTables,
                       list(read_tsv(paste0(simDir,"InvFreqs.txt"),
                                     show_col_types = FALSE)))
  }
  for (gen in unique(invFreqTables[[1]]$Generation)) {
    # print(gen)
    for (lhs in filter(invFreqTables[[1]],Generation==gen)$LifeHistoryStage){
      # print(filter(invFreqTables[[1]],Generation==gen))
      # print(lhs)
      genlhsFreqCounts <- tibble()
      # colnames(genlhsFreqCounts) <- c("Freq","Count","FreqRank","Ancestral","FreqRankIncAnc","OccurrenceCount")
      genInd <- which(invFreqTables[[1]]$Generation==gen & invFreqTables[[1]]$LifeHistoryStage==lhs)
      # print(invFreqTables[[1]]$Generation==gen)
      # print(invFreqTables[[1]]$LifeHistoryStage==lhs)
      # print(genInd)
      for (invFreqTable in invFreqTables) {
        alleleSize <- invFreqTable[[genInd,3]]*2
        # print(alleleSize)
        polyInds <- which(invFreqTable[genInd,-1:-4]%!in%c(0,alleleSize))+4
        if (length(polyInds)>0){
          # print(invFreqTable,n=200)
          # print(polyInds)
          # print(invFreqTable[genInd,polyInds])
          # print(as.numeric(invFreqTable[genInd,polyInds]))
          polyInvCounts <- as.numeric(invFreqTable[genInd,polyInds])
          # print(polyInvCounts)
          polyInvFreqs <- polyInvCounts/alleleSize
          polyInvFreqRank <- rank(-polyInvCounts)
          # print(polyInvFreqRank)
          numPolyInv <- length(polyInvCounts)
          countStd <- NULL
          if (invFreqTable[genInd,4]%!in%c(0,alleleSize)) {
            countStd <- invFreqTable[[genInd,4]]
          }
          # print(countStd)
          polyArrCounts <- c(polyInvCounts,countStd)
          # print(polyArrCounts)
          polyArrFreqRank <- rank(-polyArrCounts)
          # print(polyArrFreqRank)
          for (invInd in seq(numPolyInv)) {
            if (dim(genlhsFreqCounts)[1]==0) {
              genlhsFreqCounts <- tibble(Freq = polyInvFreqs[invInd],
                                         Count = polyInvCounts[invInd],
                                         FreqRank = polyInvFreqRank[invInd],
                                         Ancestral = FALSE,
                                         FreqRankIncAnc = polyArrFreqRank[invInd],
                                         OccurrenceCount = 1)
              
            } else{
              matchInd <- which(genlhsFreqCounts[,1]==polyInvFreqs[invInd] & 
                                  genlhsFreqCounts[,2]==polyInvCounts[invInd] & 
                                  genlhsFreqCounts[,3]==polyInvFreqRank[invInd] & 
                                  genlhsFreqCounts[,4]==FALSE & 
                                  genlhsFreqCounts[,5]==polyArrFreqRank[invInd])
              if (length(matchInd)==0)  {
                genlhsFreqCounts <- add_row(genlhsFreqCounts,
                                            Freq = polyInvFreqs[invInd],
                                            Count = polyInvCounts[invInd],
                                            FreqRank = polyInvFreqRank[invInd],
                                            Ancestral = FALSE,
                                            FreqRankIncAnc = polyArrFreqRank[invInd],
                                            OccurrenceCount = 1)
              } else if (length(matchInd)==1) {
                genlhsFreqCounts[matchInd,6] <- genlhsFreqCounts[matchInd,6]+1
              } else {
                stop("More than one SFS entry for the same count, etc:",str(matchInd))
              }
            }
          }
          if (!is.null(countStd)) {
            freqStd <- countStd/alleleSize
            matchInd <- which(genlhsFreqCounts[,1]==countStd & 
                                genlhsFreqCounts[,2]==freqStd & 
                                genlhsFreqCounts[,3]==-1 & 
                                genlhsFreqCounts[,4]==TRUE & 
                                genlhsFreqCounts[,5]==polyArrFreqRank[numPolyInv+1])
            if (length(matchInd)==0)  {
              genlhsFreqCounts <- add_row(genlhsFreqCounts,
                                          Freq = freqStd,
                                          Count = countStd,
                                          FreqRank = -1,
                                          Ancestral = TRUE,
                                          FreqRankIncAnc = polyArrFreqRank[numPolyInv+1],
                                          OccurrenceCount = 1)
            } else if (length(matchInd)==1) {
              genlhsFreqCounts[matchInd,6] <- genlhsFreqCounts[matchInd,6]+1
            } else {
              stop("More than one SFS entry for the same count, etc:",str(matchInd))
            }
          }
        }
      }
      # print(genlhsFreqCounts)
      # print(dim(genlhsFreqCounts))
      # genlhsFreqCounts <- cbind(Generation=gen,as_tibble_row(genlhsFreqCounts))
      if (dim(genlhsFreqCounts)[1]!=0){
        genlhsFreqCounts <- cbind(Generation=gen,LifeHistoryStage=lhs,genlhsFreqCounts)
        # print(genlhsFreqCounts)
        # print(specificParamFreqCounts)
        specificParamFreqCounts <- rbind(specificParamFreqCounts,genlhsFreqCounts)
      }
      # print(specificParamFreqCounts)
      # specificParamFreqCounts <- rbind(specificParamFreqCounts,genlhsFreqCounts)
    }
  }
  # print(specificParamFreqCounts)
  specificParamFreqCounts <- mutate(specificParamFreqCounts,
                                    MutRate=mutRate,
                                    InvRate=invRate,
                                    EncounterNum=encNum,
                                    CrossoverMapLength=recRate,
                                    ConversionRate=conRate,
                                    PopSize=popSize,
                                    MaxGeneration=numGens,
                                    WrightFisherSurvival=as.logical(wrightFisher),
                                    NoMaleCost=as.logical(noMaleCost),
                                    NoFemaleCost=as.logical(noFemaleCost),
                                    .before=1)
  return(specificParamFreqCounts)
}

extractInvSFSTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,nmcs,nfcs,repN){
  invSFSTable <- tibble()
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  for (noMaleCost in nmcs){
                    for (noFemaleCost in nfcs){
                      specificParamFreqCounts <- extractSpecParamInvSFSTable(resultsDir,mutRate,invRate,encNum,
                                                                             recRate,conRate,popSize,numGens,
                                                                             wrightFisher,noMaleCost,noFemaleCost,
                                                                             repN)
                      invSFSTable <- rbind(invSFSTable,specificParamFreqCounts)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(invSFSTable)
}




#
# Generating whole population statistic data for neutral inversion dynamic comparisons
#
extractPopStatRow <- function(resultsDir,mutRate,invRate,encNum,
                              recRate,conRate,popSize,numGens,
                              wrightFisher,noMaleCost,noFemaleCost,
                              repNum){
  simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                     recRate,conRate,popSize,numGens,
                     wrightFisher,noMaleCost,noFemaleCost,
                     repNum)
  print(simDir)
  specPopStatTable <- read_tsv(paste0(simDir,"PopStats.txt"),
                               show_col_types = FALSE)
  colnames(specPopStatTable)[colnames(specPopStatTable)=="PopSize"]<-"NumSelected"
  specPopStatTable <- mutate(specPopStatTable,
                             MutRate=mutRate,
                             InvRate=invRate,
                             EncounterNum=encNum,
                             CrossoverMapLength=recRate,
                             ConversionRate=conRate,
                             PopSize=popSize,
                             MaxGeneration=numGens,
                             WrightFisherSurvival=as.logical(wrightFisher),
                             NoMaleCost=as.logical(noMaleCost),
                             NoFemaleCost=as.logical(noFemaleCost),
                             ReplicateNum=repNum,
                             .before=1)
  return(specPopStatTable)
}


extractPopStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,nmcs,nfcs,repN){
  popStatTable <- tibble()
  # `%!in%` <- Negate(`%in%`)
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  for (noMaleCost in nmcs){
                    for (noFemaleCost in nfcs){
                      for (repNum in seq(repN)-1){
                        specPopStatTable <- extractPopStatRow(resultsDir,mutRate,invRate,encNum,
                                                              recRate,conRate,popSize,numGens,
                                                              wrightFisher,noMaleCost,noFemaleCost,
                                                              repNum)
                        popStatTable <- rbind(popStatTable,specPopStatTable)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(popStatTable)
}

averagePopStatsByRep <- function(popStatTable) {
  avgPopStatTable <- popStatTable %>% 
    group_by(MutRate, InvRate, EncounterNum, CrossoverMapLength, ConversionRate,
             PopSize, MaxGeneration, WrightFisherSurvival, NoMaleCost, NoFemaleCost, 
             Generation, LifeHistoryStage) %>% 
    summarise(across(SurMean:RepVar, mean), .groups="drop")
  colnames(avgPopStatTable) <- c(head(colnames(avgPopStatTable),-4),paste0("Avg",tail(colnames(avgPopStatTable),4)))
  return(avgPopStatTable)
}




#
# Generating a table of all inversions in every replicate, for further analysis
#

extractSpecParamAllInvs <- function(resultsDir,mutRate,invRate,encNum,
                                    recRate,conRate,popSize,numGens,
                                    wrightFisher,noMaleCost,noFemaleCost,
                                    repNum) {
  specificRepInvData <- tibble()
  simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                     recRate,conRate,popSize,numGens,
                     wrightFisher,noMaleCost,noFemaleCost,
                     repNum)
  invFreqFile = paste0(simDir,"InvFreqs.txt")
  invFreqs = read_tsv(invFreqFile, show_col_types = FALSE)
  invIDs <- sub('I','',colnames(invFreqs[1,-1:-3]))
  # invStatTables <- list()
  for (invID in invIDs) {
    invStatTable <- read_tsv(paste0(simDir,"Inv",invID,".txt"),
                             show_col_types = FALSE)
    specificRepInvData <- rbind(specificRepInvData,mutate(invStatTable,InvID=invID,.before=3))
  }
  specificRepInvData <- filter(specificRepInvData,Count!=0) # Std arrangements are still being added when frequency = 0, as they can 'revive' when other arrangements fix
  specificRepInvData <- mutate(specificRepInvData,
                               MutRate=mutRate,
                               InvRate=invRate,
                               EncounterNum=encNum,
                               CrossoverMapLength=recRate,
                               ConversionRate=conRate,
                               PopSize=popSize,
                               MaxGeneration=numGens,
                               WrightFisherSurvival=as.logical(wrightFisher),
                               NoMaleCost=as.logical(noMaleCost),
                               NoFemaleCost=as.logical(noFemaleCost),
                               ReplicateNum=repNum,
                               .before=1)
  return(specificRepInvData)
}

#
# Calculating the difference of each inversion from population averages
#

calcSpecParamAllInvPopAvgDiffs <- function(specParAllInv,specParPopStat) {
  # print(specParAllInv[,-seq(1,12)],n=40)
  # print(specParPopStat[,-seq(1,12)],n=40)
  specificRepInvData <- tibble()
  if (nrow(specParAllInv)==0){
    return(specificRepInvData)
  } else if (nrow(specParPopStat)==0){
    stop("No PopStat data for given Inv data: ",str(specParAllInv))
  }
  for (gen in unique(specParAllInv$Generation)) {
    genPopStat <- filter(specParPopStat,Generation == gen)
    genAllInv <- filter(specParAllInv,Generation == gen)
    # print(genPopStat)
    # print(genAllInv)
    for (lhs in unique(genAllInv$LifeHistoryStage)) {
      lhsPopStat <- filter(genPopStat,LifeHistoryStage == lhs)
      if (nrow(lhsPopStat) >= 2){
        stop("More than one population statistic entry for the same generation and life history stage:",str(lhsPopStat))
      }
      lhsAllInv <- filter(genAllInv,LifeHistoryStage == lhs)
      if (nrow(lhsAllInv) > 0){
        # print(lhsAllInv)
        # print(lhsPopStat)
        # specGenPopAvgNum <- lhsPopStat$AvgNumMut
        specGenPopAvgSur <- lhsPopStat$SurMean
        specGenPopAvgRep <- lhsPopStat$RepMean
        # specGenInvNumDiff <- lhsAllInv$AvgNumMut - specGenPopAvgNum
        specGenInvSurDiff <- lhsAllInv$AvgSurEff - specGenPopAvgSur
        specGenInvRepDiff <- lhsAllInv$AvgRepEff - specGenPopAvgRep
        # print(specGenPopAvgSur)
        # print(specGenPopAvgRep)
        # print(specGenInvSurDiff)
        # print(specGenInvRepDiff)
        
        # print(lhsAllInv$Count)
        # print((lhsPopStat$NumSelected*2))
        # print(lhsAllInv$Count/(lhsPopStat$NumSelected*2))
        specGenFreqs <- lhsAllInv$Count/(lhsPopStat$NumSelected*2)
        specGenFemaleFreqs <- lhsAllInv$FemaleCount/(lhsPopStat$NumFemales*2)
        specGenMaleFreqs <- lhsAllInv$MaleCount/(lhsPopStat$NumMales*2)
        
        specGenHomFreqs <- lhsAllInv$Homozygote/lhsPopStat$NumSelected
        specGenFemaleHomFreqs <- lhsAllInv$FemaleHomozygote/lhsPopStat$NumFemales
        specGenMaleHomFreqs <- lhsAllInv$MaleHomozygote/lhsPopStat$NumMales
        
        specGenHWExpHomFreqs <- specGenFreqs^2
        specGenHWExpHomFemaleFreqs <- specGenFemaleFreqs^2
        specGenHWExpHomMaleFreqs <- specGenMaleFreqs^2
        
        specGenHWExpHom <- specGenHWExpHomFreqs*lhsPopStat$NumSelected
        specGenHWExpHomFemale <- specGenHWExpHomFemaleFreqs*lhsPopStat$NumFemales
        specGenHWExpHomMale <- specGenHWExpHomMaleFreqs*lhsPopStat$NumMales
        
        
        specGenHeterFreqs <- lhsAllInv$Heterozygote/lhsPopStat$NumSelected
        specGenFemaleHeterFreqs <- lhsAllInv$FemaleHeterozygote/lhsPopStat$NumFemales
        specGenMaleHeterFreqs <- lhsAllInv$MaleHeterozygote/lhsPopStat$NumMales
        
        specGenHWExpHeterFreqs <- 2*specGenFreqs*(1-specGenFreqs)
        specGenHWExpHeterFemaleFreqs <- 2*specGenFemaleFreqs*(1-specGenFemaleFreqs)
        specGenHWExpHeterMaleFreqs <- 2*specGenMaleFreqs*(1-specGenMaleFreqs)
        
        specGenHWExpHeter <- specGenHWExpHeterFreqs*lhsPopStat$NumSelected
        specGenHWExpHeterFemale <- specGenHWExpHeterFemaleFreqs*lhsPopStat$NumFemales
        specGenHWExpHeterMale <- specGenHWExpHeterMaleFreqs*lhsPopStat$NumMales
        
        
        specGenFreqRank <- rank(-lhsAllInv$Count)
        specGenSurRank <- rank(-lhsAllInv$AvgSurEff)
        specGenRepRank <- rank(-lhsAllInv$AvgRepEff)
        specGenNumArrangements <- length(specGenRepRank)
        specGenNormSurRank <- -specGenSurRank+((specGenNumArrangements+1)/2)
        specGenNormRepRank <- -specGenRepRank+((specGenNumArrangements+1)/2)
        
        # specGenArrScen <- "Other"
        # intermediateArr <- specGenFreqs > .1 & specGenFreqs < .9
        # numIntermediateArr <- sum(intermediateArr)
        # hasMaleSpec <- sum(intermediateArr & specGenFemaleFreqs < .1) == 1
        # hasFemaleSpec <- sum(intermediateArr & specGenMaleFreqs < .1) == 1
        # if (numIntermediateArr == 2) {
        #   if (hasMaleSpec) {
        #     specGenArrScen <- "Male Specific"
        #   } else if (hasFemaleSpec) {
        #     specGenArrScen <- "Female Specific"
        #   } else {
        #     specGenArrScen <- "Other Two Arrangement"
        #   }
        # } else if (numIntermediateArr == 3) {
        #   if (hasMaleSpec && hasFemaleSpec) {
        #     specGenArrScen <- "Both Male and Female Specific"
        #   } else{
        #     specGenArrScen <- "Other Three Arrangement"
        #   }
        # }
        
        
        # specGenInvUpdate <- mutate(lhsAllInv,
        #                            DiffAvgNumMutFromPop=specGenInvNumDiff,
        #                            DiffAvgSurEffFromPop=specGenInvSurDiff,
        #                            DiffAvgRepEffFromPop=specGenInvRepDiff)
        specGenInvUpdate <- mutate(lhsAllInv,
                                   Frequency=specGenFreqs,
                                   FemaleFreq=specGenFemaleFreqs,
                                   MaleFreq=specGenMaleFreqs,
                                   HomozygoteFreq=specGenHomFreqs,
                                   FemaleHomozygoteFreq=specGenFemaleHomFreqs,
                                   MaleHomozygoteFreq=specGenMaleHomFreqs,
                                   HWExpectedHomozygoteFreq=specGenHWExpHomFreqs,
                                   HWExpectedFemaleHomozygoteFreq=specGenHWExpHomFemaleFreqs,
                                   HWExpectedMaleHomozygoteFreq=specGenHWExpHomMaleFreqs,
                                   HWExpectedHomozygote=specGenHWExpHom,
                                   HWExpectedFemaleHomozygote=specGenHWExpHomFemale,
                                   HWExpectedMaleHomozygote=specGenHWExpHomMale,
                                   HeterozygoteFreq=specGenHeterFreqs,
                                   FemaleHeterozygoteFreq=specGenFemaleHeterFreqs,
                                   MaleHeterozygoteFreq=specGenMaleHeterFreqs,
                                   HWExpectedHeterozygoteFreq=specGenHWExpHeterFreqs,
                                   HWExpectedFemaleHeterozygoteFreq=specGenHWExpHeterFemaleFreqs,
                                   HWExpectedMaleHeterozygoteFreq=specGenHWExpHeterMaleFreqs,
                                   HWExpectedHeterozygote=specGenHWExpHeter,
                                   HWExpectedFemaleHeterozygote=specGenHWExpHeterFemale,
                                   HWExpectedMaleHeterozygote=specGenHWExpHeterMale,
                                   FrequencyRank=specGenFreqRank,
                                   SurvivalRank=specGenSurRank,
                                   ReproductiveRank=specGenRepRank,
                                   TotalNumPolyArrangements=specGenNumArrangements,
                                   NormalizedSurvivalRank=specGenNormSurRank,
                                   NormalizedReproductiveRank=specGenNormRepRank,
                                   # ArrangementScenario=specGenArrScen,
                                   .before=15)
        specificRepInvData <- rbind(specificRepInvData,specGenInvUpdate)
      }
    }
  }
  # specificRepInvData <- mutate(specificRepInvData,
  #                              DiffAvgNumMutFromPop=invNumDiffs,
  #                              DiffAvgRepEffFromPop=invRepDiffs,
  #                              DiffAvgSurEffFromPop=invSurDiffs)
  return(specificRepInvData)
}

calcAllInvPopAvgDiffs <- function(allInvTable,popStatTable){
  invData <- tibble()
  for (m in unique(allInvTable$MutRate)){
    for (i in unique(allInvTable$InvRate)){
      for (e in unique(allInvTable$EncounterNum)){
        for (r in unique(allInvTable$CrossoverMapLength)){
          for (c in unique(allInvTable$ConversionRate)){
            for (N in unique(allInvTable$PopSize)){
              for (g in unique(allInvTable$MaxGeneration)){
                for (wfs in unique(allInvTable$WrightFisherSurvival)){
                  for (nmc in unique(allInvTable$NoMaleCost)){
                    for (nfc in unique(allInvTable$NoFemaleCost)){
                      for (repNum in unique(allInvTable$ReplicateNum)){
                        specParAllInv <- filter(allInvTable,
                                                MutRate == m,
                                                InvRate == i,
                                                EncounterNum == e,
                                                CrossoverMapLength == r,
                                                ConversionRate == c,
                                                PopSize == N,
                                                MaxGeneration == g,
                                                WrightFisherSurvival == wfs,
                                                NoMaleCost == nmc,
                                                NoFemaleCost == nfc,
                                                ReplicateNum == repNum)
                        specParPopStat <- filter(popStatTable,
                                                 MutRate == m,
                                                 InvRate == i,
                                                 EncounterNum == e,
                                                 CrossoverMapLength == r,
                                                 ConversionRate == c,
                                                 PopSize == N,
                                                 MaxGeneration == g,
                                                 WrightFisherSurvival == wfs,
                                                 NoMaleCost == nmc,
                                                 NoFemaleCost == nfc,
                                                 ReplicateNum == repNum)
                        # print(list(specParAllInv,specParPopStat))
                        specificRepInvData <- calcSpecParamAllInvPopAvgDiffs(specParAllInv,
                                                                             specParPopStat)
                        invData <- rbind(invData,specificRepInvData)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(invData)
}




#
# Generating a table of all inversions in every replicate, for further analysis
#

extractSpecParamAllArrs <- function(resultsDir,mutRate,invRate,encNum,
                                            recRate,conRate,popSize,numGens,
                                            wrightFisher,noMaleCost,noFemaleCost,
                                            repNum) {
  specificRepArrData <- tibble()
  simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                     recRate,conRate,popSize,numGens,
                     wrightFisher,noMaleCost,noFemaleCost,
                     repNum)
  arrFile = paste0(simDir,"ArrangementData.txt")
  # print(arrFile)
  # arrFile = paste0(resultsDir,"ArrangementData.txt")
  # arrFileData = read_tsv(arrFile, show_col_types = FALSE)
  arrFileData = tryCatch(read_tsv(arrFile, show_col_types = FALSE),
                         error = function(e) {
                           print(paste0("Tried opening ",arrFile," experienced error ",e))
                           return(specificRepArrData)
                           })
  
  # specificRepArrData <- rbind(specificRepArrData,arrFileData)
  
  # specificRepArrData <- filter(specificRepArrData,Count!=0) # Arrangements can 'revive' when variants fix (and so are no longer tracked) or recombination occurs
  specificRepArrData <- mutate(arrFileData,
                               MutRate=mutRate,
                               InvRate=invRate,
                               EncounterNum=encNum,
                               CrossoverMapLength=recRate,
                               ConversionRate=conRate,
                               PopSize=popSize,
                               MaxGeneration=numGens,
                               WrightFisherSurvival=as.logical(wrightFisher),
                               NoMaleCost=as.logical(noMaleCost),
                               NoFemaleCost=as.logical(noFemaleCost),
                               ReplicateNum=repNum,
                               .before=1)
  return(specificRepArrData)
}


calcSpecParamAllArrPopStats <- function(specParAllArr) {
  # print(specParAllArr[,-seq(1,12)],n=40)
  specificRepArrData <- tibble()
  if (nrow(specParAllArr)==0){
    return(specificRepArrData)
  }
  for (gen in unique(specParAllArr$Generation)) {
    genAllArr <- filter(specParAllArr,Generation == gen)
    # print(genAllArr)
    for (lhs in unique(genAllArr$LifeHistoryStage)) {
      lhsAllArr <- filter(genAllArr,LifeHistoryStage == lhs)
      if (nrow(lhsAllArr) > 0){
        # print(lhsAllArr)
        numAlleles <- sum(lhsAllArr$Count)
        numFemAlleles <- sum(lhsAllArr$FemaleCount)
        numMalAlleles <- sum(lhsAllArr$MaleCount)
        
        specGenPopAvgSur <- sum(lhsAllArr$AvgSurEff*lhsAllArr$Count/numAlleles)
        specGenPopAvgRep <- sum(lhsAllArr$AvgRepEff*lhsAllArr$Count/numAlleles)
        
        specGenArrSurDiff <- lhsAllArr$AvgSurEff - specGenPopAvgSur
        specGenArrRepDiff <- lhsAllArr$AvgRepEff - specGenPopAvgRep
        
        print(specGenPopAvgSur)
        print(specGenPopAvgRep)
        print(specGenArrSurDiff)
        print(specGenArrRepDiff)
        
        # print((lhsPopStat$NumSelected*2))
        # print(lhsAllArr$Count/(lhsPopStat$NumSelected*2))
        specGenFreqs <- lhsAllArr$Count/numAlleles
        specGenFemaleFreqs <- lhsAllArr$FemaleCount/numFemAlleles
        specGenMaleFreqs <- lhsAllArr$MaleCount/numMalAlleles
        
        specGenHomFreqs <- lhsAllArr$Homozygote/(numAlleles/2)
        specGenFemaleHomFreqs <- lhsAllArr$FemaleHomozygote/(numFemAlleles/2)
        specGenMaleHomFreqs <- lhsAllArr$MaleHomozygote/(numMalAlleles/2)
        
        specGenHWExpHomFreqs <- specGenFreqs^2
        specGenHWExpHomFemaleFreqs <- specGenFemaleFreqs^2
        specGenHWExpHomMaleFreqs <- specGenMaleFreqs^2
        
        specGenHWExpHom <- specGenHWExpHomFreqs*numAlleles/2
        specGenHWExpHomFemale <- specGenHWExpHomFemaleFreqs*numFemAlleles/2
        specGenHWExpHomMale <- specGenHWExpHomMaleFreqs*numMalAlleles/2
        
        
        specGenHeterFreqs <- lhsAllArr$Heterozygote/(numAlleles/2)
        specGenFemaleHeterFreqs <- lhsAllArr$FemaleHeterozygote/(numFemAlleles/2)
        specGenMaleHeterFreqs <- lhsAllArr$MaleHeterozygote/(numMalAlleles/2)
        
        specGenHWExpHeterFreqs <- 2*specGenFreqs*(1-specGenFreqs)
        specGenHWExpHeterFemaleFreqs <- 2*specGenFemaleFreqs*(1-specGenFemaleFreqs)
        specGenHWExpHeterMaleFreqs <- 2*specGenMaleFreqs*(1-specGenMaleFreqs)
        
        specGenHWExpHeter <- specGenHWExpHeterFreqs*(numAlleles/2)
        specGenHWExpHeterFemale <- specGenHWExpHeterFemaleFreqs*(numFemAlleles/2)
        specGenHWExpHeterMale <- specGenHWExpHeterMaleFreqs*(numMalAlleles/2)
        
        
        specGenFreqRank <- rank(-lhsAllArr$Count)
        specGenSurRank <- rank(-lhsAllArr$AvgSurEff)
        specGenRepRank <- rank(-lhsAllArr$AvgRepEff)
        specGenNumArrangements <- length(specGenRepRank)
        specGenNormSurRank <- -specGenSurRank+((specGenNumArrangements+1)/2)
        specGenNormRepRank <- -specGenRepRank+((specGenNumArrangements+1)/2)
        
        
        specGenArrScen <- singleArrangementState(specGenFreqs,
                                                 specGenFemaleFreqs,
                                                 specGenMaleFreqs)
        specGenArrScen <- factor(specGenArrScen,
                                 levels=c("Other",
                                          "Favors Female",
                                          "Female-Specific",
                                          "Male-Specific",
                                          "Both Male- and Female-Specific",
                                          "Three Arrangement"))
        
        specGenArrUpdate <- mutate(lhsAllArr,
                                   Frequency=specGenFreqs,
                                   FemaleFreq=specGenFemaleFreqs,
                                   MaleFreq=specGenMaleFreqs,
                                   HomozygoteFreq=specGenHomFreqs,
                                   FemaleHomozygoteFreq=specGenFemaleHomFreqs,
                                   MaleHomozygoteFreq=specGenMaleHomFreqs,
                                   HWExpectedHomozygoteFreq=specGenHWExpHomFreqs,
                                   HWExpectedFemaleHomozygoteFreq=specGenHWExpHomFemaleFreqs,
                                   HWExpectedMaleHomozygoteFreq=specGenHWExpHomMaleFreqs,
                                   HWExpectedHomozygote=specGenHWExpHom,
                                   HWExpectedFemaleHomozygote=specGenHWExpHomFemale,
                                   HWExpectedMaleHomozygote=specGenHWExpHomMale,
                                   HeterozygoteFreq=specGenHeterFreqs,
                                   FemaleHeterozygoteFreq=specGenFemaleHeterFreqs,
                                   MaleHeterozygoteFreq=specGenMaleHeterFreqs,
                                   HWExpectedHeterozygoteFreq=specGenHWExpHeterFreqs,
                                   HWExpectedFemaleHeterozygoteFreq=specGenHWExpHeterFemaleFreqs,
                                   HWExpectedMaleHeterozygoteFreq=specGenHWExpHeterMaleFreqs,
                                   HWExpectedHeterozygote=specGenHWExpHeter,
                                   HWExpectedFemaleHeterozygote=specGenHWExpHeterFemale,
                                   HWExpectedMaleHeterozygote=specGenHWExpHeterMale,
                                   FrequencyRank=specGenFreqRank,
                                   SurvivalRank=specGenSurRank,
                                   ReproductiveRank=specGenRepRank,
                                   TotalNumPolyArrangements=specGenNumArrangements,
                                   NormalizedSurvivalRank=specGenNormSurRank,
                                   NormalizedReproductiveRank=specGenNormRepRank,
                                   ArrangementScenario=specGenArrScen,
                                   .before=15)
        specificRepArrData <- rbind(specificRepArrData,specGenArrUpdate)
      }
    }
  }
  return(specificRepArrData)
}

calcAllArrPopStats <- function(allArrTable){
  arrData <- tibble()
  for (m in unique(allArrTable$MutRate)){
    for (i in unique(allArrTable$InvRate)){
      for (e in unique(allArrTable$EncounterNum)){
        for (r in unique(allArrTable$CrossoverMapLength)){
          for (c in unique(allArrTable$ConversionRate)){
            for (N in unique(allArrTable$PopSize)){
              for (g in unique(allArrTable$MaxGeneration)){
                for (wfs in unique(allArrTable$WrightFisherSurvival)){
                  for (nmc in unique(allArrTable$NoMaleCost)){
                    for (nfc in unique(allArrTable$NoFemaleCost)){
                      for (repNum in unique(allArrTable$ReplicateNum)){
                        specParAllArr <- filter(allArrTable,
                                                MutRate == m,
                                                InvRate == i,
                                                EncounterNum == e,
                                                CrossoverMapLength == r,
                                                ConversionRate == c,
                                                PopSize == N,
                                                MaxGeneration == g,
                                                WrightFisherSurvival == wfs,
                                                NoMaleCost == nmc,
                                                NoFemaleCost == nfc,
                                                ReplicateNum == repNum)
                        specificRepArrData <- calcSpecParamAllArrPopStats(specParAllArr)
                        arrData <- rbind(arrData,specificRepArrData)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(arrData)
}


singleArrangementStateFromTable <- function(invPopAvgDiffTable) {
  # print(invPopAvgDiffTable)
  intermediateArr <- invPopAvgDiffTable$Frequency > .1 & invPopAvgDiffTable$Frequency < .9
  numIntermediateArr <- sum(intermediateArr)
  hasMaleSpec <- sum(intermediateArr & invPopAvgDiffTable$FemaleFreq < .1) == 1
  hasFemaleSpec <- sum(intermediateArr & invPopAvgDiffTable$MaleFreq < .1) == 1
  if (numIntermediateArr == 2) {
    if (hasMaleSpec) {
      return("Male-Specific")
    } else if (hasFemaleSpec) {
      return("Female-Specific")
    } else {
      return("Other Two Arrangement")
    }
  } else if (numIntermediateArr == 3) {
    if (hasMaleSpec && hasFemaleSpec) {
      return("Both Male and Female Specific")
    } else{
      return("Other Three Arrangement")
    }
  } else {
    return("Other")
  }
}

singleArrangementState <- function(Frequency,FemaleFreq,MaleFreq) {
  intermediateArr <- Frequency > .1 & Frequency < .9
  numIntermediateArr <- sum(intermediateArr)
  hasFemaleSpec <- sum(intermediateArr & MaleFreq < .1) == 1
  hasMaleSpec <- sum(intermediateArr & FemaleFreq < .1) == 1
  # print(paste0(numIntermediateArr,'/t',hasFemaleSpec,'/t',hasMaleSpec,'/n'))
  if (numIntermediateArr == 2) {
    if (hasMaleSpec) {
      return(as.factor("Male-Specific"))
    } else if (hasFemaleSpec) {
      return(as.factor("Female-Specific"))
    } else {
      # return(as.factor("Other Two Arrangement"))
      return(as.factor("Favors Female"))
    }
  } else if (numIntermediateArr == 3) {
    if (hasMaleSpec && hasFemaleSpec) {
      return(as.factor("Both Male- and Female-Specific"))
    } else{
      return(as.factor("Three Arrangement"))
    }
  } else {
    return(as.factor("Other"))
  }
}

classifyArrangementState <- function(invData) {
  stateTable <- invData %>% 
    group_by(MutRate, InvRate, EncounterNum, CrossoverMapLength, ConversionRate,
             PopSize, MaxGeneration, WrightFisherSurvival, NoMaleCost, NoFemaleCost, 
             ReplicateNum,Generation, LifeHistoryStage) %>% 
    summarise(ArrangementScenario=singleArrangementState(Frequency,
                                                         FemaleFreq,
                                                         MaleFreq),
                                                .groups="drop")
  # print(stateTable$ArrangementScenario)
  classifiedArrTable <- full_join(invData,stateTable,
                                  by=c("MutRate", "InvRate", "EncounterNum", 
                                       "CrossoverMapLength", "ConversionRate", 
                                       "PopSize", "MaxGeneration", 
                                       "WrightFisherSurvival", "NoMaleCost", "NoFemaleCost", 
                                       "ReplicateNum", "Generation", "LifeHistoryStage"),
                                  suffix=c("",""),
                                  relationship="many-to-one") %>%
    relocate(ArrangementScenario, .before = Count)
  classifiedArrTable$ArrangementScenario <- factor(classifiedArrTable$ArrangementScenario,
                                                   levels=c("Other",
                                                            "Favors Female",
                                                            "Female-Specific",
                                                            "Male-Specific",
                                                            "Both Male- and Female-Specific",
                                                            "Three Arrangement"))
  # print(classifiedArrTable$ArrangementScenario)
  return(classifiedArrTable)
}




distParLsRecursion <- function(i,inParLs,outParLs,currPars,singlePars) {
  if (i > length(inParLs)) {
    for (j in seq(length(currPars))) {
      if (j > length(outParLs)) {
        # newParVect <- c(outParLs[[j]],currPars[[j]])
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
  # print(inParLs)
  singlePars <- list()
  for (par in inParLs) {
    singlePars <- c(singlePars,length(par)==1)
  }
  # print(singlePars)
  outParLs <- distParLsRecursion(1,inParLs,list(),list(),singlePars)
  return(outParLs)
}

parallelTreeRbind <- function() {
  
}

parallelExtract <- function(extractFunc,inParLs,finishingFunc=NULL) {
  numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  parLs <- c(extractFunc,distributedParamList(inParLs),
             SIMPLIFY = FALSE,
             mc.cores = numCores)
  # print(parLs)
  # totLen <- 1
  # for (v in list(ms,is,es,rs,cs,sizes,gens,wfs,nmcs,nfcs,repNs)) {
  #   totLen <- totLen*length(v)
  # }
  extractResults <- do.call(mcmapply,parLs)
  # print(extractResults)
  statTable <- Reduce(rbind,extractResults)
  # print(statTable)
  if (!is.null(finishingFunc)){
    statTable <- finishingFunc(statTable)
  }
  return(statTable)
}

getTableParamLs <- function(statTable,parVect){
  # symVect <- c()
  # for (strName in parVect) {
  #   symVect <- c(symVect,sym(strName))
  # }
  # print(symVect)
  parLs <- list()
  for (name in parVect){
    parLs <- c(parLs,list(unique(statTable[[name]])))
  }
  return(parLs)
}

getSimParamLs <- function(statTable,includeReps=FALSE,includeGens=FALSE){
  expectedNames <- c("MutRate","InvRate","EncounterNum",
                     "CrossoverMapLength","ConversionRate","PopSize",
                     "MaxGeneration","WrightFisherSurvival",
                     "NoMaleCost","NoFemaleCost")
  if (includeReps) {
    expectedNames <- c(expectedNames,"ReplicateNum")
  }
  if (includeGens) {
    expectedNames <- c(expectedNames,"Generation")
  }
  inParLs <- getTableParamLs(statTable,expectedNames)
  # es <- unique(statTable$EncounterNum)
  # rs <- unique(statTable$CrossoverMapLength)
  # cs <- unique(statTable$ConversionRate)
  # Ns <- unique(statTable$PopSize)
  # gens <- unique(statTable$MaxGeneration)
  # wfs <- unique(statTable$WrightFisherSurvival)
  # nmcs <- unique(statTable$NoMaleCost)
  # nfcs <- unique(statTable$NoFemaleCost)
  # inParLs <- list(ms,is,es,rs,cs,Ns,gens,wfs,nmcs,nfcs)
  # if (includeReps) {
  #   repNums <- unique(statTable$ReplicateNum)
  #   inParLs <- c(inParLs,list(repNums))
  # }
  # if (includeGens) {
  #   t <- unique(statTable$Generation)
  #   inParLs <- c(inParLs,list(t))
  # }
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
  # # print(tableLs)
  # groupBySym <- list()
  # for (strName in groupByVect) {
  #   groupBySym <- c(groupBySym,sym(strName))
  # }
  # print(groupBySym)
  # # print(tableLs[[1]])
  # # print(c(list(tableLs[[1]]),groupBySym))
  # filteredTableLs <- list()
  # for (table in tableLs) {
  #   groupedTable <- table
  #   for (groupSym in groupBySym){
  #     groupedTable <- group_by(groupedTable,!!groupSym,.add=TRUE)
  #   }
  #   # print(groupedTable)
  #   filteredTable <- group_split(groupedTable)
  #   filteredTableLs <- c(filteredTableLs,list(filteredTable))
  # }
  return(filteredTableLs)
}

parallelAnalysis <- function(calcFunc, statTableLs, inParLs=list(), finishingFunc=NULL) {
  numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  # ms <- unique(statTable$MutRate)
  # is <- unique(statTable$InvRate)
  # es <- unique(statTable$EncounterNum)
  # rs <- unique(statTable$CrossoverMapLength)
  # cs <- unique(statTable$ConversionRate)
  # Ns <- unique(statTable$PopSize)
  # gens <- unique(statTable$MaxGeneration)
  # wfs <- unique(statTable$WrightFisherSurvival)
  # nmcs <- unique(statTable$NoMaleCost)
  # nfcs <- unique(statTable$NoFemaleCost)
  # t <- unique(statTable$Generation)
  # totLen <- 1
  # for (v in list(m,i,e,r,c,N,g,wf,nmc,nfc,t)) {
  #   totLen <- totLen*length(v)
  # }
  # parLs <- distributedParamList(list(ms,is,es,rs,cs,Ns,gens,wfs,nmcs,nfcs,t))
  # 
  # extractResults <- mcmapply(calcFunc,list(statTable),
  #                            parLs[[1]],
  #                            parLs[[2]],
  #                            parLs[[3]],
  #                            parLs[[4]],
  #                            parLs[[5]],
  #                            parLs[[6]],
  #                            parLs[[7]],
  #                            parLs[[8]],
  #                            parLs[[9]],
  #                            parLs[[10]],
  #                            parLs[[11]],
  #                            SIMPLIFY = FALSE,
  #                            mc.cores = numCores)
  
  
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


# For splitting a table into a (not easily divisible) number of cores
# WILL THROW AN ERROR if there are fewer distinct subsets than cores and maxSplitInd is not set
varCoreFilteredTableLs <- function(tableLs,numCores=NULL,maxSplitInd=NULL){
  if (is.null(numCores)){
    numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  }
  if (!is.null(maxSplitInd)){
    maxCores <- nrow(distinct(tableLs[[1]][,seq(maxSplitInd-1)]))
    print(maxCores)
    numCores <- min(numCores,maxCores)
    print(numCores)
  } else {
    cat(paste0("WARNING: \nmaxSplitInd not set, this will fail if there are fewer \n",
               "distinct subsets than cores among the parameter \n",
               "columns you are trying to split over. \n",
               "It is safer to specify maxSplitInd as the last column to split over."))
  }
  print(paste0("Cores allocated: ",numCores))
  splitIndex <- 1
  # print(tableLs[[1]][,seq(splitIndex)])
  totDistinct <- nrow(distinct(tableLs[[1]][,seq(splitIndex)]))
  while (totDistinct < numCores && splitIndex < ncol(tableLs[[1]])-1 && (is.null(maxSplitInd) || splitIndex < maxSplitInd-1)) {
    splitIndex <- splitIndex+1
    totDistinct <- nrow(distinct(tableLs[[1]][,seq(splitIndex)]))
  }
  splitChunkSize <- as.integer(totDistinct/numCores)
  splitChunkRemainder <- totDistinct%%numCores
  allDistinct <- distinct(tableLs[[1]][,seq(splitIndex)])
  groupVect <- colnames(allDistinct)
  # print(splitIndex)
  # print(totDistinct)
  # print(nrow(allDistinct))
  
  # indexedFullTable <- tableLs[[1]] %>%
  #   group_by(A, B, D) %>%
  #   mutate(index = cur_group_id()) %>%
  #   ungroup()
  
  splitTableLs <- list()
  for (t in seq(length(tableLs))){
    splitTableLs <- c(splitTableLs,list(list()))
  }
  # print(splitTableLs)
  for (i in seq(numCores)){
    # print(paste0("Preparing filtered tables for core ",i))
    if (i > splitChunkRemainder){
      distStart <- 1+splitChunkRemainder*(splitChunkSize+1)+(i-splitChunkRemainder-1)*splitChunkSize
      distEnd <- splitChunkRemainder*(splitChunkSize+1)+(i-splitChunkRemainder)*splitChunkSize
      if (i == numCores && distEnd != totDistinct){
        stop(paste0("varCoreFilteredTableLs attempted to split into ",totDistinct,
                    " parameter combinations, but the last split index was only ",distEnd))
      }
    } else{
      distStart <- 1+(i-1)*(splitChunkSize+1)
      distEnd <- i*(splitChunkSize+1)
    }
    
    specParFilteredTableLs <- list()
    for (t in seq(length(tableLs))){
      specParFilteredTableLs <- c(specParFilteredTableLs,list(list()))
    }
    # print('Init specParFilteredTableLs:')
    # print(specParFilteredTableLs)
    for (j in seq(distStart,distEnd)){
      subLs <- list()
      for (k in seq(splitIndex)){
        subLs <- c(subLs,list(substitute(name==val,
                                         list(name=sym(groupVect[k]),
                                              val=allDistinct[[j,k]]))))
      }
      # print('Substituted parameters:')
      # print(subLs)
      for (l in seq(length(tableLs))){
        # print(paste0("Filtering table ",l))
        specParFilteredTable <- do.call(filter,c(list(tableLs[[l]]),subLs))
        # print(specParFilteredTable)
        specParFilteredTableLs[[l]] <- rbind(specParFilteredTableLs[[l]],specParFilteredTable)
      }
    }
    # print('Combined Filtered Tables for this Core:')
    # print(specParFilteredTableLs)
    for (m in seq(length(tableLs))){
      # splitTableLs[[m]] <- c(splitTableLs[[m]],list(specParFilteredTableLs[[m]]))
      splitTableLs[[m]] <- c(splitTableLs[[m]],specParFilteredTableLs[m])
    }
    # print(splitTableLs)
  }
  # print(splitTableLs)
  return(splitTableLs)
}


# parallelAnalysisVarCore <- function(calcFunc, statTableLs, 
#                                     inParLs=list(), finalParIndex=NULL, 
#                                     finishingFunc=NULL) {
#   numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
#   groupParLs <- getTableParamLs(statTableLs[[1]],groupByVect)
#   filteredTableLs <- varCoreFilteredTableLs(statTableLs,numCores,maxSplitInd=finalParIndex)
#   
#   parLs <- c(calcFunc,statTableLs,distributedParamList(inParLs),
#              SIMPLIFY = FALSE,
#              mc.cores = numCores)
#   # print(parLs)
#   extractResults <- do.call(mcmapply,parLs)
#   
#   calcTable <- Reduce(rbind,extractResults)
#   if (!is.null(finishingFunc)){
#     calcTable <- finishingFunc(calcTable)
#   }
#   return(calcTable)
# }


# varCoreParallelTableAnalysis <- function(tableLs,numCores=NULL,
#                                          maxSplitInd=NULL,
#                                          finishingFunc=NULL){
#   if (is.null(numCores)){
#     numCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
#   }
#   splitIndex <- 1
#   totDistinct <- nrow(distinct(tableLs[[1]][,seq(splitIndex)]))
#   while (totDistinct < numCores && splitIndex < nrow(tableLs[[1]])-1 && (is.null(maxSplitInd)) || splitIndex < maxSplitInd-1) {
#     splitIndex <- splitIndex+1
#     totDistinct <- nrow(distinct(tableLs[[1]][,seq(splitIndex)]))
#   }
#   splitChunkSize <- as.integer(totDistinct/numCores)
#   splitChunkRemainder <- totDistinct%%numCores
#   allDistinct <- tableLs[[1]][,seq(splitIndex)]
#   groupVect <- colnames(allDistinct)
#   
#   splitTableLs <- list()
#   for (i in seq(length(tableLs))){
#     splitTableLs <- c(splitTableLs,list(list()))
#   }
#   print(splitTableLs)
#   for (i in seq(numCores)){
#     subLs <- list()
#     if (i > splitChunkRemainder){
#       distStart <- 1+splitChunkRemainder*(splitChunkSize+1)+(i-splitChunkRemainder-1)*splitChunkSize
#       distEnd <- splitChunkRemainder*(splitChunkSize+1)+(i-splitChunkRemainder)*splitChunkSize
#       if (i == numCores && distEnd != totDistinct){
#         stop(paste0("varCoreFilteredTableLs attempted to split into ",totDistinct,
#                     " parameter combinations, but the last split index was only ",distEnd))
#       }
#     } else{
#       distStart <- 1+(i-1)*(splitChunkSize+1)
#       distEnd <- i*(splitChunkSize+1)
#     }
#     
#     specParFilteredTableLs <- list()
#     for (i in seq(length(tableLs))){
#       specParFilteredTableLs <- c(specParFilteredTableLs,list(list()))
#     }
#     for (j in seq(distStart,distEnd)){
#       params <- allDistinct[j,]
#       for (k in seq(splitIndex)){
#         subLs <- c(subLs,list(substitute(name==val,
#                                          list(name=sym(groupVect[k]),
#                                               val=allDistinct[j,k]))))
#         # print(groupParGrid[i,j])
#       }
#       print(subLs)
#       for (l in seq(length(tableLs))){
#         specParFilteredTable <- do.call(filter,c(list(tableLs[[l]]),subLs))
#         specParFilteredTableLs[[l]] <- rbind(specParFilteredTableLs[[l]],specParFilteredTable)
#       }
#     }
#     for (m in seq(length(tableLs))){
#       # splitTableLs[[m]] <- c(splitTableLs[[m]],list(specParFilteredTableLs[[m]]))
#       splitTableLs[[m]] <- c(splitTableLs[[m]],specParFilteredTableLs[m])
#     }
#   }
#   return()
# }
