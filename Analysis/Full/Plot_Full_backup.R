# For analyzing and plotting full SAIsim model data 
#

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/Full")

# library(readr)
# library(ggplot2)
library(tidyverse)
library(svglite)

readLastLineToTable <- function(file) {
  read_tsv(pipe(paste('(head -n 1; tail -n 1) <', file)), show_col_types = FALSE)
}

extractNonFixedInvs <- function(invFreqs,alleleSize){
  nonZeroFreqs <- invFreqs[-1][which(invFreqs[-1]!=0)]
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
    for (testGen in otherInvs$Generation){
      if (testGen != gen){
        print(paste0("Warning: not all generations are the same: ",otherInvs$Generation))
      }
    }
    counts = otherInvs$Count
    total = sum(counts)
    weights = counts/total
    avgOthers <- colSums(subset(otherInvs, select=-c(InvID,Generation,Count))*weights)
    avgOtherRow <- mutate(as_tibble_row(avgOthers),Generation=gen,InvID="Other",Count=total,.before=1)
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
                      wrightFisher,repNum){
  m <- formatC(mutRate, digits = 2, width = 3, format = "e", flag = "0")
  tempSimDir <- paste0('mut',m)
  i <- formatC(mutRateInv, digits = 2, width = 3, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'inv',i)
  e <- formatC(encNum, width = 3, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'enc',e)
  r <- formatC(crossRate, digits = 2, width = 3, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'cro',r)
  c <- formatC(convRate, digits = 2, width = 3, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'con',c)
  N <- formatC(popSize, digits = 1, width = 2, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'pop',N)
  g <- formatC(numGens, digits = 1, width = 2, format = "e", flag = "0")
  tempSimDir <- paste0(tempSimDir,'gen',g)
  wf <- formatC(wrightFisher, digits = 1, width = 1, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'wf',wf)
  n <- formatC(repNum, width = 4, format = "d", flag = "0")
  tempSimDir <- paste0(tempSimDir,'rep',n)
  return(paste0(resultsDir,tempSimDir,'/'))
}

extractFinalGenInvStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,repN){
  invStatTable <- tibble()
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  # popSize <- 10^size
                  alleleSize <- 2*popSize
                  for (repNum in seq(repN)-1){
                    simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                                       recRate,conRate,popSize,numGens,
                                       wrightFisher,repNum)
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
                                      ReplicateNum=repNum,
                                      .before=1)
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
  invStatTable <- mutate(invStatTable,Frequency=invStatTable$Count/(2*invStatTable$PopSize),.before=13)
  return(invStatTable)
}

extractAvgStatsAllGenInv <- function(simDir,alleleSize){
  invFreqFile = paste0(simDir,"InvFreqs.txt")
  invFreqs = read_tsv(invFreqFile, show_col_types = FALSE)
  invIDs <- sub('I','',colnames(invFreqs[1,-1]))
  invStatTables <- list()
  for (invID in invIDs) {
    invStatTables <- c(invStatTables,
                       list(read_tsv(paste0(simDir,"Inv",invID,".txt"),
                                     show_col_types = FALSE)))
  }
  names(invStatTables) <- invIDs
  avgInvStatsAllGens <- tibble()
  for (gen in invFreqs$Generation){
    nonFixedInvs <- extractNonFixedInvs(filter(invFreqs,Generation == gen),alleleSize)
    invStats <- tibble()
    for (ID in nonFixedInvs){
      invStats <- rbind(invStats,filter(invStatTables[[ID]],Generation == gen))
    }
    invStats <- mutate(invStats,InvID=nonFixedInvs,.before=2)
    avgStatsThisGen <- avgNonMaxInvs(invStats)
    avgInvStatsAllGens <- rbind(avgInvStatsAllGens,avgStatsThisGen)
  }
  return(avgInvStatsAllGens)
}

extractInvStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,repN){
  invStatTable <- tibble()
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  alleleSize <- 2*popSize
                  for (repNum in seq(repN)-1){
                    simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                                       recRate,conRate,popSize,numGens,
                                       wrightFisher,repNum)
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
                                      ReplicateNum=repNum,
                                      .before=1)
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
  invStatTable <- mutate(invStatTable,Frequency=invStatTable$Count/(2*invStatTable$PopSize),.before=13)
  return(invStatTable)
}

calcDiffsTable <- function(fullFinalGenStatTable){
  diffColNames <- paste0('AbsDiff',colnames(select(fullFinalGenStatTable,
                                                      !c(MutRate:Frequency))))
  diffsTable <- tibble()
  for (m in unique(fullFinalGenStatTable$MutRate)){
    for (i in unique(fullFinalGenStatTable$InvRate)){
      for (e in unique(fullFinalGenStatTable$EncounterNum)){
        for (r in unique(fullFinalGenStatTable$CrossoverMapLength)){
          for (c in unique(fullFinalGenStatTable$ConversionRate)){
            for (N in unique(fullFinalGenStatTable$PopSize)){
              for (g in unique(fullFinalGenStatTable$MaxGeneration)){
                for (w in unique(fullFinalGenStatTable$WrightFisherSurvival)){
                  for (gen in unique(fullFinalGenStatTable$Generation)){
                    specificParamData <- filter(fullFinalGenStatTable,
                                                MutRate == m,
                                                InvRate == i,
                                                EncounterNum == e,
                                                CrossoverMapLength == r,
                                                ConversionRate == c,
                                                PopSize == N,
                                                MaxGeneration == g,
                                                WrightFisherSurvival == w,
                                                Generation == gen)
                    for (num in unique(specificParamData$ReplicateNum)){
                      specificReplicateGen <- filter(specificParamData,
                                                     ReplicateNum == num)
                      specificParams <- select(specificReplicateGen[1,],
                                               c(MutRate:Generation))
                      countFreqData <- select(specificReplicateGen[1,],
                                          c(Count,Frequency))
                      colnames(countFreqData) <- c("MaxArrangementCount","MaxArrangementFrequency")
                      specificReplicateData <- select(specificReplicateGen,
                                                      !c(MutRate:Frequency))
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
                      diffsTable <- rbind(diffsTable,diffsRow)
                      print(specificParams)
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

calcAvgDiffsTable <- function(fullFinalGenStatTable){
  diffColNames <- paste0('AvgAbsDiff',colnames(select(fullFinalGenStatTable[1,],
                                                      !c(MutRate:Frequency))))
  avgDiffsTable <- tibble()
  for (m in unique(fullFinalGenStatTable$MutRate)){
    for (i in unique(fullFinalGenStatTable$InvRate)){
      for (e in unique(fullFinalGenStatTable$EncounterNum)){
        for (r in unique(fullFinalGenStatTable$CrossoverMapLength)){
          for (c in unique(fullFinalGenStatTable$ConversionRate)){
            for (N in unique(fullFinalGenStatTable$PopSize)){
              for (g in unique(fullFinalGenStatTable$MaxGeneration)){
                for (w in unique(fullFinalGenStatTable$WrightFisherSurvival)){
                  for (gen in unique(fullFinalGenStatTable$Generation)){
                    specificParamGenData <- filter(fullFinalGenStatTable,
                                                MutRate == m,
                                                InvRate == i,
                                                EncounterNum == e,
                                                CrossoverMapLength == r,
                                                ConversionRate == c,
                                                PopSize == N,
                                                MaxGeneration == g,
                                                WrightFisherSurvival == w,
                                                Generation == gen)
                    if (nrow(specificParamGenData) > 0){
                      specificParams <- select(specificParamGenData[1,],
                                               c(MutRate:WrightFisherSurvival,Generation))
                      specificParamDiffs <- tibble()
                      numMonomorphic <- 0
                      numPolymorphic <- 0
                      for (num in unique(specificParamGenData$ReplicateNum)){
                        specificReplicate <- filter(specificParamGenData,
                                                    ReplicateNum == num)
                        specificReplicateData <- select(specificReplicate,
                                                        !(MutRate:Frequency))
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
                      avgDiffsTable <- rbind(avgDiffsTable,avgDiffsRow)
                      print(specificParams)
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


axisLabels <- c(
  "AvgAbsDiffAvgNumMut"="Magnitude of the Difference in Number of Mutants",
  "AvgAbsDiffAvgSurEff"="Magnitude of the Difference in Total Survival Effect",
  "AvgAbsDiffAvgRepEff"="Magnitude of the Difference in Total Reproductive Effect",
  "AbsDiffAvgNumMut"="Magnitude of the Difference in Number of Mutants",
  "AbsDiffAvgSurEff"="Magnitude of the Difference in Total Survival Effect",
  "AbsDiffAvgRepEff"="Magnitude of the Difference in Total Reproductive Effect",
  "AvgNumMut"="Average Number of Mutations in the Arrangement",
  "AvgSurEff"="Average Survival Effect of the Arrangement",
  "AvgRepEff"="Average Reproductive Effect of the Arrangement",
  "AbsDiffCount"="Difference in Count between the Most Dominant Arrangement and all Others",
  "AvgAbsDiffCount"="Average Difference in Count between the Most Dominant Arrangement and all Others",
  "MaxArrangementCount"="Population Count of the Most Common Arrangement",
  "AvgMaxArrangementCount"="Average Population Count of the Most Common Arrangement",
  "AbsDiffFrequency"="Difference in Frequency between the Most Dominant Arrangement and all Others",
  "AvgAbsDiffFrequency"="Average Difference in Frequency between the Most Dominant Arrangement and all Others",
  "MaxArrangementFrequency"="Population Frequency of the Most Common Arrangement",
  "AvgMaxArrangementFrequency"="Average Population Frequency of the Most Common Arrangement",
  "Generation"="Simulated Generation",
  "InvCount"="Inversion Arrangement Count within the Population",
  "InvFreq"="Frequency of the Inversion Arrangement within the Population",
  # "InvCount-Generation"="Number of Occurrences across all Simulations",
  # "Generation-InvCount"="Number of Occurrences across all Simulations",
  # "InvFreq-Generation"="Number of Occurrences across all Simulations",
  # "Generation-InvFreq"="Number of Occurrences across all Simulations",
  "InvCount-Generation"="Count",
  "Generation-InvCount"="Count",
  "InvFreq-Generation"="Count",
  "Generation-InvFreq"="Count",
  "SurMean"="Mean Survival Value",
  "SurVar"="Variance in Survival Value",
  "RepMean"="Mean Display Value",
  "RepVar"="Variance in Display Value"
  )


scatterPanelPlot <- function(simData, 
                             axis1var, 
                             axis2var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             colorvar=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  saveFile <- paste0(outPrefix,'Scatter.',axis1var,'.',axis2var)
  if (is.null(colorvar)){
    p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else {
    c <- sym(colorvar)
    p <- ggplot(simData, aes(x=!!a1, y=!!a2, color=!!c))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  
  p <- p +
    # geom_point(aes(x=log10(!!a1), y=log10(!!a2)),size=.25, shape=21) +
    geom_point(size=.25, shape=21) +
    xlab(axisLabels[axis1var]) +
    ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  if (!is.null(colorvar)){
    maxColorVarVal <- max(select(simData,colorvar))
    minColorVarVal <- min(select(simData,colorvar))
    # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
    
    p <- p +
      scale_colour_gradientn(
        name=axisLabels[colorvar],
        # name=NULL,
        limits = c(minColorVarVal,maxColorVarVal),
        colors = hcl.colors(5, palette = "Plasma"),
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      # scale_colour_gradient(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "#4B0055",
      #   high = "#FDE333",
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) +
      # scale_colour_gradient2(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "red",
      #   mid = "white",
      #   high = "blue",
      #   midpoint = midColorVarVal,
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) + 
      guides(
        # size = "none",
        colour = guide_colourbar(title.position = "right")
      )+
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  
    # ggtitle("") +
    # expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    # theme(legend.justification=c(1,0),
    #       legend.position=c(0.28,0.8))               # Position legend in bottom right
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




scatterPlotFinalGenParams <- function(fullFinalGenStatTable, 
                                      axis1var, 
                                      axis2var,
                                      colorvar=NULL,
                                      outPrefix = "./",
                                      width = 6,
                                      height = 5){
  for (m in unique(fullFinalGenStatTable$MutRate)){
    for (i in unique(fullFinalGenStatTable$InvRate)){
      for (e in unique(fullFinalGenStatTable$EncounterNum)){
        for (r in unique(fullFinalGenStatTable$CrossoverMapLength)){
          for (c in unique(fullFinalGenStatTable$ConversionRate)){
            for (N in unique(fullFinalGenStatTable$PopSize)){
              specificParamData <- filter(fullFinalGenStatTable,
                                          MutRate == m,
                                          InvRate == i,
                                          EncounterNum == e,
                                          CrossoverMapLength == r,
                                          ConversionRate == c,
                                          PopSize == N)
              specificParams <- select(specificParamData[1,],
                                       c(MutRate:PopSize))
              specificOutPrefix <- paste0(outPrefix,'m',m,'i',i,'e',e,'r',r,'c',c,'N',N,'.')
              print(specificOutPrefix)
              scatterPanelPlot(specificParamData, 
                                       axis1var, 
                                       axis2var,
                                       colorvar = colorvar,
                                       outPrefix = specificOutPrefix,
                                       width = width,
                                       height = height)
            }
          }
        }
      }
    }
  }
  return()
}




violinPanelPlot <- function(simData, 
                            axis1var,
                            axis2vars,
                            facet1vars=NULL,
                            facet2vars=NULL,
                            color=FALSE,
                            outPrefix = "./",
                            width = 6,
                            height = 5) {
  a1 <- sym(axis1var)
  a21 <- sym(axis2vars[1])
  saveFile <- paste0(outPrefix,'Violin.',axis1var,'.',axis2vars[1])
  
  dualAxis <- (length(axis2vars) > 1)
  reqVars <- c(axis1var,facet1vars,facet2vars)
  tempSimData <- simData
  scale_shift_f <- function(x, scaleCoeff, shiftVal) {
    return(x * scaleCoeff - shiftVal)
  }
  inv_scale_shift_f <- function(x, scaleCoeff, shiftVal) {
    return( (x + shiftVal)/scaleCoeff )
  }
    
  if (dualAxis) {
    saveFile <- paste0(saveFile,'.',axis2vars[2])
    # a22 <- sym(axis2vars[2])
    maxA21 <- max(simData[axis2vars[1]])
    maxA22 <- max(simData[axis2vars[2]])
    minA21 <- min(simData[axis2vars[1]])
    minA22 <- min(simData[axis2vars[2]])
    scaleCoeff <- (maxA22-minA22)/(maxA21-minA21)
    shiftVal <- minA21-minA22
    # scaledVar <- paste0(axis2vars[2],"Scaled")
    scaledValVar <- "TempScaledVals"
    sourceVarCol <- "ValSource"
    sV <- sym(sourceVarCol)
    
    scaledSimData1 <- select(simData,all_of(c(reqVars,axis2vars[1])))
    scaledSimData1 <- cbind(scaledSimData1,axisLabels[axis2vars[1]])
    colnames(scaledSimData1) <- c(reqVars,scaledValVar,sourceVarCol)
    
    scaledSimData2 <- select(simData,all_of(c(reqVars,axis2vars[2])))
    scaledSimData2 <- cbind(scaledSimData2,axisLabels[axis2vars[2]])
    colnames(scaledSimData2) <- c(reqVars,scaledValVar,sourceVarCol)
    scaledSimData2[scaledValVar] <- inv_scale_shift_f(scaledSimData2[scaledValVar], scaleCoeff, shiftVal)
    
    tempSimData <- rbind(scaledSimData1,scaledSimData2)
    a22 <- sym(scaledValVar)
  }
  
  
  # saveFile <- paste0(outPrefix,'Histogram.',axis1var,'.',axis2var)
  
  # print(head(tempSimData,20))
  # print(tail(tempSimData,20))
  
  p <- ggplot(tempSimData, aes(x=!!a1))
  
  if (dualAxis) {
    p <- p +
      geom_violin(aes(y=!!a22, 
                      group=interaction(!!a1,!!sV), 
                      color=!!sV),
                  scale="count",
                  position="dodge") +
      # geom_boxplot(aes(y=!!a22, group=interaction(!!a1,!!sV), color=!!sV)) +
      geom_boxplot(aes(y=!!a22, 
                       group=interaction(!!a1,!!sV),
                       color=!!sV),
                   outlier.shape=1) +
                   # outlier.size=0.1,
                   # width=500,
                   # position=position_dodge(width=1000)) +
      # geom_boxplot(aes(y=!!a22, group=interaction(!!a1,!!sV), color=!!sV),width=~.*0.5,position=position_dodge(width=~.*1.1)) +
      scale_y_continuous(
        name = axisLabels[axis2vars[1]],
        sec.axis = sec_axis(~scale_shift_f(., scaleCoeff, shiftVal), name = axisLabels[axis2vars[2]])
      )
    # print(ggplot_build(p))
  } else {
    p <- p +
      geom_violin(aes(y=!!a21, group=!!a1)) +
      geom_boxplot(aes(y=!!a22, group=!!a1)) +
      ylab(axisLabels[axis2vars[1]])
  }
  
  p <- p +
    xlab(axisLabels[axis1var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  
  
  if (color) {
    saveFile <- paste0(saveFile,'.Colored')
    p <- p +
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




#
# Generating inversion SFS plots for neutral inversion dynamic comparisons
#

extractInvSFSTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,repN){
  invSFSTable <- tibble()
  `%!in%` <- Negate(`%in%`)
  for (mutRate in ms){
    for (invRate in is){
      for (encNum in es){
        for (recRate in rs){
          for (conRate in cs){
            for (popSize in sizes){
              for (numGens in gens){
                for (wrightFisher in wfs){
                  alleleSize <- 2*popSize
                  specificParamFreqCounts <- tibble()
                  # colnames(specificParamFreqCounts) <- c("Generation","InvCount","InvFreq","OccurrenceCount")
                  invFreqTables <- list()
                  for (repNum in seq(repN)-1) {
                    simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                                       recRate,conRate,popSize,numGens,
                                       wrightFisher,repNum)
                    print(simDir)
                    invFreqTables <- c(invFreqTables,
                                       list(read_tsv(paste0(simDir,"InvFreqs.txt"),
                                                     show_col_types = FALSE)))
                  }
                  for (gen in invFreqTables[[1]]$Generation) {
                    genFreqCounts <- tibble()
                    genInd <- which(invFreqTables[[1]]$Generation==gen)
                    # print(genInd)
                    # print(genFreqCounts)
                    for (invFreqTable in invFreqTables) {
                      for (invCount in invFreqTable[genInd,which(invFreqTable[genInd,-1:-2]%!in%c(0,alleleSize))+2]) {
                        if (dim(genFreqCounts)[1]==0 || !any(genFreqCounts[,1]==invCount)) {
                          genFreqCounts <- rbind(genFreqCounts,c(invCount,invCount/alleleSize,1))
                        } else {
                          genFreqCounts[which(genFreqCounts[,1] == invCount),3] <- genFreqCounts[which(genFreqCounts[,1] == invCount),3]+1
                        }
                      }
                    }
                    colnames(genFreqCounts) <- c("InvCount","InvFreq","OccurrenceCount")
                    # print(genFreqCounts)
                    # print(dim(genFreqCounts))
                    # genFreqCounts <- cbind(Generation=gen,as_tibble_row(genFreqCounts))
                    if (dim(genFreqCounts)[1]!=0){
                      genFreqCounts <- cbind(Generation=gen,genFreqCounts)
                    }
                    # print(genFreqCounts)
                    specificParamFreqCounts <- rbind(specificParamFreqCounts,genFreqCounts)
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
                                                    .before=1)
                  # print(specificParamFreqCounts)
                  invSFSTable <- rbind(invSFSTable,specificParamFreqCounts)
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


freqBinPanelPlot <- function(simData, 
                             axis1var, 
                             axis2var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             weightvar=NULL,
                             binwidths=NULL,
                             # colorvar=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  saveFile <- paste0(outPrefix,'FreqPlot.',axis1var,'.',axis2var)
  
  if (is.null(weightvar)){
    p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else {
    w <- sym(weightvar)
    p <- ggplot(simData, aes(x=!!a1, y=!!a2, weight=!!w))
  }
  
  if (is.null(binwidths)) {
    p <- p +
      geom_bin2d()
  } else {
    p <- p +
      geom_bin2d(binwidth=binwidths)
  }
  p <- p +
    xlab(axisLabels[axis1var]) +
    ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  p <- p +
    scale_fill_gradientn(
      name=axisLabels[paste0(axis1var,'-',axis2var)],
      # name=NULL,
      # limits = c(minColorVarVal,maxColorVarVal),
      colors = c("#FFFFFF",rev(hcl.colors(5, palette = "Blues")),rep(hcl.colors(1, palette = "Blues"),8)),
      # colors = c("#FFFFFF",rev(hcl.colors(4, palette = "Plasma")),rep(hcl.colors(1, palette = "Plasma"),3)),
      # colors = c("#FFFFFF",rev(hcl.colors(3, palette = "Plasma"))),
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    guides(
      # size = "none",
      colour = guide_colourbar(title.position = "right")
    )+
    theme(
      legend.direction = "vertical",
      legend.position = "right",
      # legend.key.height = unit(height/20, "in"),
      # legend.key.height = unit(2, "in"),
      legend.title = element_text(size = 10, angle = 90),
      legend.title.align = 0.5
    )
    # guides(
    #   colour = guide_colourbar(title.position = "right")
    # )+
    # theme(legend.title = element_text(size = 10, angle = 90),
    #   legend.title.align = 0.5,
    #   legend.key.height = unit(height/20, "in"),
    #   # legend.key.align = 0.5,
    #   legend.position = "right",
    #   legend.direction = "vertical"
    # )
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




histogramPanelPlot <- function(simData, 
                             axis1var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             weightvar=NULL,
                             colorvar=NULL,
                             binwidths=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  # a2 <- sym(axis2var)
  # saveFile <- paste0(outPrefix,'Histogram.',axis1var,'.',axis2var)
  saveFile <- paste0(outPrefix,'Histogram.',axis1var)
  
  if (is.null(weightvar) && is.null(colorvar)){
    p <- ggplot(simData, aes(x=!!a1))
    # p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else if (!is.null(weightvar) && is.null(colorvar)) {
    w <- sym(weightvar)
    p <- ggplot(simData, aes(x=!!a1, weight=!!w))
  } else if (is.null(weightvar) && !is.null(colorvar)) {
    c <- sym(colorvar)
    p <- ggplot(simData, aes(x=!!a1, colour=!!w))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  } else {
    c <- sym(colorvar)
    w <- sym(weightvar)
    p <- ggplot(simData, aes(x=!!a1, colour=!!w, weight=!!w))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  
  if (is.null(binwidths)) {
    p <- p +
      geom_histogram()
  } else {
    p <- p +
      geom_histogram(binwidth = binwidths)
  }
  
  p <- p +
    xlab(axisLabels[axis1var]) +
    # ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  
  
  if (!is.null(colorvar)){
    maxColorVarVal <- max(select(simData,colorvar))
    minColorVarVal <- min(select(simData,colorvar))
    # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
    p <- p +
      scale_colour_gradientn(
        name=axisLabels[paste0(colorvar)],
        # name=NULL,
        limits = c(minColorVarVal,maxColorVarVal),
        colors = c("#FFFFFF",rev(hcl.colors(5, palette = "Plasma"))),
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "fill"
      ) +
      guides(
        colour = guide_colourbar(title.position = "right")
      )+
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




#
# Generating whole population statistc plots for neutral inversion dynamic comparisons
#

extractPopStatTable <- function(resultsDir,ms,is,es,rs,cs,sizes,gens,wfs,repN){
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
                  # alleleSize <- 2*popSize
                  # specificParamPopStatCounts <- tibble()
                  # popStatFileTables <- list()
                  for (repNum in seq(repN)-1) {
                    simDir = getSimDir(resultsDir,mutRate,invRate,encNum,
                                       recRate,conRate,popSize,numGens,
                                       wrightFisher,repNum)
                    print(simDir)
                    specPopStatTable <- read_tsv(paste0(simDir,"PopStats.txt"),
                                                 show_col_types = FALSE)
                    specPopStatTable <- mutate(specPopStatTable,
                                               MutRate=mutRate,
                                               InvRate=invRate,
                                               EncounterNum=encNum,
                                               CrossoverMapLength=recRate,
                                               ConversionRate=conRate,
                                               PopSize=popSize,
                                               MaxGeneration=numGens,
                                               WrightFisherSurvival=as.logical(wrightFisher),
                                               ReplicateNum=repNum,
                                               .before=1)
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
  return(popStatTable)
}

averagePopStatsByRep <- function(popStatTable) {
  avgPopStatTable <- popStatTable %>% 
    group_by(MutRate, InvRate, EncounterNum, CrossoverMapLength, ConversionRate,
             PopSize, MaxGeneration, WrightFisherSurvival, Generation) %>% 
    summarise(across(SurMean:RepVar, mean), .groups="drop")
  colnames(avgPopStatTable) <- c(head(colnames(avgPopStatTable),-4),paste0("Avg",tail(colnames(avgPopStatTable),4)))
  return(avgPopStatTable)
}


#
#
# Load and process the data
#

resultsDir = '/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Full/'
N3dir = paste0(resultsDir,'N3/')
N4dir = paste0(resultsDir,'N4/')
N2dir = paste0(resultsDir,'N2/')


if (!exists("fullFinalGenStatTable")){
  if (file.exists(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))){
    load(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))
  } else {
    if (!exists("fullFinalGenN3StatTable")){
      fullFinalGenN3StatTable <- extractFinalGenInvStatTable(N3dir,c(1e-3,1e-5),1e-4,c(10,100),1,
                                                     c(1e-2,1e-5),1e3,2e4,1,1000)
    }
    if (!exists("partE10FinalGenN4StatTable")){
      partE10FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-3,1e-5),1e-4,c(10),1,
                                                        c(1e-2,1e-5),1e4,2e4,1,100)
    }
    if (!exists("partE100m5FinalGenN4StatTable")){
      partE100m5FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-5),1e-4,c(100),1,
                                                           c(1e-2,1e-5),1e4,2e4,1,100)
    }
    if (!exists("partE100m3c2FinalGenN4StatTable")){
      partE100m3c2FinalGenN4StatTable <- extractFinalGenInvStatTable(N4dir,c(1e-3),1e-4,c(100),1,
                                                             c(1e-2),1e4,2e4,1,100)
    }
    if (!exists("partE10FinalGenN2StatTable")){
      partE10FinalGenN2StatTable <- extractFinalGenInvStatTable(N2dir,c(1e-3,1e-5),1e-4,c(10),1,
                                                        c(1e-2,1e-5),1e2,2e4,1,1000)
    }
    fullFinalGenStatTable <- rbind(fullFinalGenN3StatTable,
                                   partE10FinalGenN2StatTable,
                                   partE10FinalGenN4StatTable,
                                   partE100m5FinalGenN4StatTable,
                                   partE100m3c2FinalGenN4StatTable)
    save(fullFinalGenStatTable,file=paste0(resultsDir,"fullFinalGenStatTable.Rdata"))
  }
}


if (!exists("finalGenDiffsTable")){
  if (file.exists(paste0(resultsDir,"finalGenDiffsTable.Rdata"))){
    load(paste0(resultsDir,"finalGenDiffsTable.Rdata"))
  } else {
    finalGenDiffsTable <- calcDiffsTable(fullFinalGenStatTable)
    save(finalGenDiffsTable,file=paste0(resultsDir,"finalGenDiffsTable.Rdata"))
  }
}
if (!exists("finalGenAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))
  } else {
    finalGenAvgDiffsTable <- calcAvgDiffsTable(fullFinalGenStatTable)
    save(finalGenAvgDiffsTable,file=paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))
  }
}



scatterPlotFinalGenParams(finalGenDiffsTable,
                          "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
                          colorvar="MaxArrangementFrequency")

scatterPanelPlot(finalGenDiffsTable,
                 "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
                 c("PopSize","EncounterNum"),
                 c("MutRate","ConversionRate"),
                 colorvar="MaxArrangementFrequency",
                 width = 10,
                 height = 5)


scatterPanelPlot(filter(finalGenDiffsTable,
                        ConversionRate == 1e-2,
                        PopSize == 1e3),
                 "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
                 c("EncounterNum"),
                 c("MutRate"),
                 colorvar="MaxArrangementFrequency",
                 width = 10,
                 height = 5)





if (!exists("fullAllGenStatTable")){
  if (file.exists(paste0(resultsDir,"fullAllGenStatTable.Rdata"))){
    load(paste0(resultsDir,"fullAllGenStatTable.Rdata"))
  } else {
    if (!exists("fullAllGenN3StatTable")){
      fullAllGenN3StatTable <- extractInvStatTable(N3dir,c(1e-3,1e-5),1e-4,c(10,100),1,
                                                   c(1e-2,1e-5),1e3,2e4,1,1000)
    }
    # if (!exists("fullAllGenN4StatTable")){
    #   fullAllGenN4StatTable <- extractInvStatTable(N4dir,c(1e-3,1e-5),1e-4,c(10,100),1,
    #                                                c(1e-2,1e-5),1e4,2e4,1,100)
    # }
    # if (!exists("partE10AllGenN2StatTable")){
    #   partE10AllGenN2StatTable <- extractInvStatTable(N2dir,c(1e-3,1e-5),1e-4,c(10),1,
    #                                                     c(1e-2,1e-5),1e2,2e4,1,1000)
    # }
    # fullAllGenStatTable <- rbind(fullFinalGenN3StatTable,
    #                                partE10AllGenN2StatTable,
    #                                fullAllGenN4StatTable)
    fullAllGenStatTable <- fullAllGenN3StatTable
    save(fullAllGenStatTable,file=paste0(resultsDir,"fullAllGenStatTable.Rdata"))
  }
}

if (!exists("allGenDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenDiffsTable.Rdata"))
  } else {
    allGenDiffsTable <- calcDiffsTable(fullAllGenStatTable)
    save(allGenDiffsTable,file=paste0(resultsDir,"allGenDiffsTable.Rdata"))
  }
}
if (!exists("allGenAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))
  } else {
    allGenAvgDiffsTable <- calcAvgDiffsTable(fullAllGenStatTable)
    save(allGenAvgDiffsTable,file=paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))
  }
}

testAllGenDiffsTable <- calcDiffsTable(filter(fullAllGenStatTable,
                                              MutRate==1e-3,
                                              ConversionRate==1e-2,
                                              EncounterNum==100))
                                              # Generation %in% seq(0,20000,by=2000)))

testAllGenAvgDiffsTable <- calcAvgDiffsTable(filter(fullAllGenStatTable,
                                              MutRate==1e-3,
                                              ConversionRate==1e-2,
                                              EncounterNum==100))
                                              # Generation %in% seq(0,20000,by=2000)))

scatterPlotFinalGenParams(testAllGenDiffsTable,
                          "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
                          colorvar="Generation")

scatterPlotFinalGenParams(testAllGenAvgDiffsTable,
                          "AvgAbsDiffAvgSurEff","AvgAbsDiffAvgRepEff",
                          colorvar="Generation")

violinPanelPlot(filter(testAllGenDiffsTable,Generation %in% seq(0,20000,by=1000)),
                "Generation",
                c("AbsDiffAvgSurEff","AbsDiffAvgRepEff"),
                width = 11,
                height = 3,
                outPrefix = "./m3c2N3")



if (!exists("i4e100c2N3InvFreqTable")){
  if (file.exists(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))){
    load(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))
  } else {
    i4e100c2N3InvFreqTable <- extractInvSFSTable(N3dir,c(0,1e-3),c(1e-4),c(100),1,
                                            c(1e-2),1e3,2e4,1,1000)
    save(i4e100c2N3InvFreqTable,file=paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))
  }
}

freqBinPanelPlot(i4e100c2N3InvFreqTable,
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 8,
                 height = 2,
                 outPrefix = "./MutRateControl.")

freqBinPanelPlot(filter(i4e100c2N3InvFreqTable,
                        Generation<=4000),
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 9,
                 height = 3,
                 outPrefix = "./MutRateControl.Short.")


histogramPanelPlot(filter(i4e100c2N3InvFreqTable,
                          Generation==20000),
                   "InvFreq",
                   facet2vars = "MutRate",
                   weightvar = "OccurrenceCount",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./MutRateControl.")


if (!exists("m3e100c2N3PopStatTable")){
  if (file.exists(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))){
    load(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
  } else {
    m3e100c2N3PopStatTable <- extractPopStatTable(N3dir,c(1e-3),c(0,1e-4),c(100),1,
                                            c(1e-2),1e3,2e4,1,1000)
    save(m3e100c2N3PopStatTable,file=paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
  }
}


scatterPanelPlot(m3e100c2N3PopStatTable,
                         "SurVar",
                         "RepVar",
                         facet2vars = "InvRate",
                         colorvar = "Generation",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")


scatterPanelPlot(filter(m3e100c2N3PopStatTable, Generation==20000),
                         "SurVar",
                         "RepVar",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")



scatterPanelPlot(m3e100c2N3PopStatTable,
                         "SurMean",
                         "RepMean",
                         facet2vars = "InvRate",
                         colorvar = "Generation",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")


scatterPanelPlot(filter(m3e100c2N3PopStatTable, Generation==20000),
                         "SurMean",
                         "RepMean",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,4000,by=200)),
                "Generation",
                c("SurMean","RepMean"),
                facet2vars = "InvRate",
                width = 11,
                height = 4,
                outPrefix = "./InvRateControl.")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==0),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate0e+0")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==1e-4),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate1e-4")


repAvgm3e100c2N3PopStatTable <- averagePopStatsByRep(m3e100c2N3PopStatTable)

scatterPanelPlot(repAvgm3e100c2N3PopStatTable,
                         "AvgSurMean",
                         "AvgRepMean",
                         facet2vars = "InvRate",
                         colorvar = "Generation",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")
