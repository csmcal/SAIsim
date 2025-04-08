# SAIsim EffectGrid Heatmap Generator Script 
#


getSimFinalFreq <- function(fileName,popSize) {
  mutCounts <- read.table(fileName,sep='\t',header=TRUE)
  finalCounts <- tail(mutCounts$Count,n=3)
  finalFreq <- finalCounts[1]/(2*popSize)
  return(finalFreq)
}

getAvgFreq <- function(survEffect,reprEffect,filePrefix,popSize,numReps) {
  s <- format(survEffect,nsmall = 3)
  #print(s)
  r <- format(reprEffect,nsmall = 3)
  #print(r)
  fileBase <- paste(filePrefix,"s",s,"r",r,"n", sep="")
  #print(fileBase)
  finalFreqs <- c()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,".txt",sep="")
      finalFreqs <- c(finalFreqs,getSimFinalFreq(fileName,popSize))
    })
  }
  #print(finalFreqs)
  #print(mean(finalFreqs))
  return(mean(finalFreqs))
}

genFreqGrid <- function(filePrefix,surEffects,repEffects,popSize,numReps) {
  eqFreqGrid <- matrix(nrow=length(surEffects),
                       ncol=length(repEffects),
                       dimnames=list(surEffects,repEffects))
  for (s in surEffects) {
    for (r in repEffects) {
      eqFreqGrid[paste0(s),paste0(r)] <- getAvgFreq(s,r,filePrefix,popSize,numReps)
    }
  }
  return(eqFreqGrid)
}

# freqHeatMap <- function(eqFreqMatrix) {
#   heatmap(eqFreqMatrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="none", margins=c(5,10),
#           main="Equilibrium Allele Counts in N = 1000", xlab="Reproductive Effect",ylab="Survival Effect")
# }

getPropPoly <- function(survEffect,reprEffect,filePrefix,popSize,numReps) {
  s <- format(survEffect,nsmall = 3)
  #print(s)
  r <- format(reprEffect,nsmall = 3)
  #print(r)
  fileBase <- paste(filePrefix,"s",s,"r",r,"n", sep="")
  #print(fileBase)
  equilFreqs <- c()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,".txt",sep="")
      equilFreqs <- c(equilFreqs,getSimFinalFreq(fileName,popSize))
    })
  }
  propPoly = sum(equilFreqs>0 & equilFreqs<1)/numReps
  return(propPoly)
}

genPolyGrid <- function(filePrefix,surEffects,repEffects,popSize,numReps) {
  eqPolyGrid <- matrix(nrow=length(surEffects),
                       ncol=length(repEffects),
                       dimnames=list(surEffects,repEffects))
  for (s in surEffects) {
    for (r in repEffects) {
      eqPolyGrid[paste0(s),paste0(r)] <- getPropPoly(s,r,filePrefix,popSize,numReps)
    }
  }
  return(eqPolyGrid)
}



#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(freqGridMatrix,file="freqGridMatrix.Rdata")
#load("freqGridMatrix.Rdata")

# Heatmap script
library(tidyverse)
# library(RColorBrewer)


freq_palette <- colorRampPalette(c("cornflowerblue","green", "yellow"))(n = 3999)
basic_freq_palette <- c("cornflowerblue","green", "yellow")
poly_palette <- colorRampPalette(c("white","cornflowerblue"))(n = 3999)
paul_tol_default <- c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933', '#AA4499')
paul_tol_ordered <- c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499')
paul_tol_modified <- c('#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', 'white')
# paul_tol_modified_rev <- c('white', '#DDCC77', '#999933', '#117733', '#44AA99', '#88CCEE', '#332288')
paul_tol_modified_rev <- rev(paul_tol_modified)
paul_tol_NA <- c('#DDDDDD')

flatMatrix <- function(matrix) {
  flatMat <-
    matrix %>% 
    as_tibble(rownames = "Survival_Effect") %>%
    pivot_longer(-c(Survival_Effect), names_to = "Reproductive_Effect", values_to = "Frequency") %>%
    mutate(
      Survival_Effect = as.numeric(Survival_Effect),
      Reproductive_Effect = as.numeric(Reproductive_Effect)
    )
  return(flatMat)
}


effGridHeatmap <- function(effGridMatrix) {
  p <- flatMatrix(effGridMatrix) %>% 
    ggplot(aes(Reproductive_Effect, Survival_Effect, fill=Frequency)) + 
    geom_tile() +
    labs(title=NULL,x=NULL,y=NULL) +
    # scale_fill_gradientn(colours = c("cornflowerblue","green", "yellow")) +
    theme(legend.position="none")
  return(p)
}


# Putting it all together in single functions

plotFreqHeatmap <- function(dataDir,plotPrefix,surEffects,repEffects,popSize,numReps,recalcData=FALSE){
  if (!recalcData && file.exists(paste0(plotPrefix,"FreqMatrix.Rdata"))){
    load(paste0(plotPrefix,"FreqMatrix.Rdata"))
  } else {
    freqMatrix <- genFreqGrid(dataDir,surEffects,repEffects,popSize,numReps)
    save(freqMatrix,file=paste0(plotPrefix,"FreqMatrix.Rdata"))
  }
  freqPlot <- effGridHeatmap(freqMatrix) + scale_fill_gradientn(limits=c(0,1),colours = basic_freq_palette)
  ggsave(paste0(plotPrefix," Freq Heatmap.png"),plot=freqPlot,units="in", width=2.5, height=2.5, dpi=600, device = 'png')
}

plotPolyHeatmap <- function(dataDir,plotPrefix,surEffects,repEffects,popSize,numReps,recalcData=FALSE){
  if (!recalcData && file.exists(paste0(plotPrefix,"PolyMatrix.Rdata"))){
    load(paste0(plotPrefix,"PolyMatrix.Rdata"))
  } else {
    polyMatrix <- genPolyGrid(dataDir,surEffects,repEffects,popSize,numReps)
    save(polyMatrix,file=paste0(plotPrefix,"PolyMatrix.Rdata"))
  }
  polyPlot <- effGridHeatmap(polyMatrix) + scale_fill_gradientn(limits=c(0,1),colors=paul_tol_modified_rev)
  ggsave(paste0(plotPrefix," Poly Heatmap.png"),plot=polyPlot,units="in", width=2.5, height=2.5, dpi=600, device = 'png')
}



# Running data collection and plotting in parallel

library(parallel)
# 
# distParLsRecursion <- function(i,inParLs,outParLs,currPars,singlePars) {
#   if (i > length(inParLs)) {
#     for (j in seq(length(currPars))) {
#       if (j > length(outParLs)) {
#         # newParVect <- c(outParLs[[j]],currPars[[j]])
#         outParLs <- append(outParLs,c(currPars[[j]]))
#       } else {
#         if (!singlePars[[j]]){
#           outParLs[[j]] <- c(outParLs[[j]],currPars[[j]])
#         }
#       }
#     }
#   } else {
#     for (p in inParLs[[i]]) {
#       outParLs <- distParLsRecursion(i+1,inParLs,outParLs,c(currPars,list(p)),singlePars)
#     }
#   }
#   return(outParLs)
# }
# 
# distributedParamList <- function(inParLs) {
#   if (length(inParLs) == 0){
#     return(inParLs)
#   }
#   # print(inParLs)
#   singlePars <- list()
#   for (par in inParLs) {
#     singlePars <- c(singlePars,length(par)==1)
#   }
#   # print(singlePars)
#   outParLs <- distParLsRecursion(1,inParLs,list(),list(),singlePars)
#   return(outParLs)
# }

coreProportion = 0.6
minCores = 2

parallelPlot <- function(plotFunc,inParLs,finishingFunc=NULL) {
  maxReqCores <- length(inParLs[[1]])
  print(maxReqCores)
  numGuessCores <- max(c(as.integer(coreProportion*detectCores()),minCores))
  numCores <- min(numGuessCores,maxReqCores)
  parLs <- c(plotFunc,inParLs,
             SIMPLIFY = FALSE,
             mc.cores = numCores)
  # print(parLs)
  plotReturns <- do.call(mcmapply,parLs)
  print(plotReturns)
  if (!is.null(finishingFunc)){
    finishingFunc(plotReturns)
  }
  return()
}




# # Running the script on inputs
# 
#
# # Set working directory
# setwd("/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/MutPolyHeatmap/")
# dataDirSTD = "../../Results/SingMutGrid/EffectGrid/"
# dataDirNMC = "../../Results/SingMutGrid/EffectGridNMC/"
# dataDirNFC = "../../Results/SingMutGrid/EffectGridNFC/"
# dataDirE10 = "../../Results/SingMutGrid/EffectGridE10/"
# dataDirN4 = "../../Results/SingMutGrid/EffectGridN4/"
# surEffects <- seq(0.0,1.0-0.025,by=0.025)
# repEffects <- seq(0.025,1.0,by=0.025)
# numReps <- 500
# popSize <- 1000
# 
# 
# effectGridFreqMatrix <- genFreqGrid(dataDirSTD,surEffects,repEffects,popSize,numReps)
# effectGridFreqNMCMatrix <- genFreqGrid(dataDirNMC,surEffects,repEffects,popSize,numReps)
# effectGridFreqNFCMatrix <- genFreqGrid(dataDirNFC,surEffects,repEffects,popSize,numReps)
# effectGridFreqE10Matrix <- genFreqGrid(dataDirE10,surEffects,repEffects,popSize,numReps)
# 
# ggsave("STD Freq Heatmap.png",plot=effGridHeatmap(effectGridFreqMatrix) + scale_fill_gradientn(colours = basic_freq_palette),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("NMC Freq Heatmap.png",plot=effGridHeatmap(effectGridFreqNMCMatrix) + scale_fill_gradientn(colours = basic_freq_palette),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("NFC Freq Heatmap.png",plot=effGridHeatmap(effectGridFreqNFCMatrix) + scale_fill_gradientn(colours = basic_freq_palette),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("E10 Freq Heatmap.png",plot=effGridHeatmap(effectGridFreqE10Matrix) + scale_fill_gradientn(colours = basic_freq_palette),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# 
# 
# effectGridMatrix <- genPolyGrid(dataDirSTD,surEffects,repEffects,popSize,numReps)
# effectGridNMCMatrix <- genPolyGrid(dataDirNMC,surEffects,repEffects,popSize,numReps)
# effectGridNFCMatrix <- genPolyGrid(dataDirNFC,surEffects,repEffects,popSize,numReps)
# effectGridE10Matrix <- genPolyGrid(dataDirE10,surEffects,repEffects,popSize,snumReps)
# 
# # ggsave("STD Heatmap.png",plot=effGridHeatmap(effectGridMatrix)+scale_fill_gradient2(low="white",mid="yellow",high="magenta2",midpoint=0.6),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("STD Heatmap.png",plot=effGridHeatmap(effectGridMatrix)+scale_fill_gradientn(colors=paul_tol_modified_rev),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("NMC Heatmap.png",plot=effGridHeatmap(effectGridNMCMatrix)+scale_fill_gradientn(colors=paul_tol_modified_rev),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("NFC Heatmap.png",plot=effGridHeatmap(effectGridNFCMatrix)+scale_fill_gradientn(colors=paul_tol_modified_rev),units="in", width=2.5, height=2.5, dpi=600, device = 'png')
# ggsave("E10 Heatmap.png",plot=effGridHeatmap(effectGridE10Matrix)+scale_fill_gradientn(colors=paul_tol_modified_rev),units="in", width=2.5, height=2.5, dpi=600, device = 'png')

# ggsave("STD Heatmap.png",plot=effGridHeatmap(effectGridMatrix)+
#          scale_fill_gradientn(colors=paul_tol_modified_rev,
#                               limits = c(0,1), breaks = seq(0,1,by=.2),
#                               guide=guide_colorbar(title=NULL)
#                               )+
#          # theme(legend.position="right",legend.text = element_text(angle=90))+
#          theme(legend.position="right",legend.text = element_text(angle=90, hjust=.5),
#                legend.key.width  = unit(1, "lines"),
#                legend.key.height = unit(3, "lines")
#                ),
#        units="in", width=4, height=4, dpi=600, device = 'png')


