# SAIsim EffectGrid Heatmap Generator Script
#

# Set working directory
setwd("/Users/cmcallester/Documents/Pool Lab/SAIsim/Analysis")

getAvgLast1000 <- function(fileName) {
  mutFreq <- read.table(fileName,sep='\t',header=TRUE)
  avgTailFreq <- mean(tail(mutFreq$Count,n=1000))
  return(avgTailFreq)
}

getAvgEquil <- function(survEffect,reprEffect,filePrefix) {
  s <- format(survEffect,nsmall = 2)
  #print(s)
  r <- format(reprEffect,nsmall = 2)
  #print(r)
  fileBase <- paste(filePrefix,"s",s,"r",r,"n", sep="")
  #print(fileBase)
  equilFreqs <- c()
  for (num in seq(1000)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,".txt",sep="")
      equilFreqs <- c(equilFreqs,getAvgLast1000(fileName))
    })
  }
  #print(equilFreqs)
  #print(mean(equilFreqs))
  return(mean(equilFreqs))
}

genEquilFreqGrid <- function(filePrefix) {
  surEffects <- seq(0.05,1.0,by=0.05)
  repEffects <- seq(0.00,1.0,by=0.05)
  eqFreqGrid <- matrix(nrow=20,ncol=21,dimnames=list(surEffects,repEffects))
  for (s in surEffects) {
    for (r in repEffects) {
      eqFreqGrid[paste0(s),paste0(r)] <- getAvgEquil(s,r,filePrefix)
    }
  }
  return(eqFreqGrid)
}

equilHeatMap <- function(eqFreqMatrix) {
  heatmap(eqFreqMatrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="none", margins=c(5,10),
          main="Equilibrium Allele Counts in N = 1000", xlab="Reproductive Effect",ylab="Survival Effect")
}

freqGridMatrix <- genEquilFreqGrid("../Results/EffectGrid/")
equilHeatMap(freqGridMatrix)
#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(freqGridMatrix,file="freqGridMatrix.Rdata")
#load("freqGridMatrix.Rdata")

#fit <- lm(mutFreq$Count ~ I(1/(mutFreq$Generation+1)))
#fit <- nls(mutFreq$Count~mutFreq$Generation)
#plot(mutFreq)
#points(mutFreq$Generation, predict(fit), type="l", col="red", lwd=2)
