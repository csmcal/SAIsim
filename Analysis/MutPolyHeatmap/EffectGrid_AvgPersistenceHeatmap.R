# SAIsim EffectGrid Heatmap Generator Script
#

# Set working directory
setwd("/Users/cmcallester/Documents/Pool Lab/SAIsim/Analysis")

getPersistence <- function(fileName) {
  mutFreq <- read.table(fileName,sep='\t',header=TRUE)
  fixGen <- which(mutFreq$Count==2000)[1]
  if (!is.na(fixGen)) {
    return(mutFreq$Generation[fixGen])
  }
  lossGen <- which(mutFreq$Count==0)[1]
  if (!is.na(lossGen)) {
    return(mutFreq$Generation[lossGen])
  }
  return(tail(mutFreq$Generation, n=1))
}

# getPersistence <- function(fileName) {
#   mutFreq <- read.table(fileName,sep='\t',header=TRUE)
#   i <- 0
#   len <- length(mutFreq$Count)
#   while ((0 < mutFreq$Count[i]) & (mutFreq$Count[i] < 2000) & (i < len)) {
#     i <- i+1
#   }
#   estimatedPersistence <- mutFreq$Generation[i]
#   return(avgTailFreq)
# }

getAvgPersistence <- function(survEffect,reprEffect,filePrefix) {
  s <- format(survEffect,nsmall = 2)
  #print(s)
  r <- format(reprEffect,nsmall = 2)
  #print(r)
  fileBase <- paste(filePrefix,"s",s,"r",r,"n", sep="")
  #print(fileBase)
  persistences <- c()
  for (num in seq(1000)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,".txt",sep="")
      persistences <- c(persistences,getPersistence(fileName))
    })
  }
  #print(persistences)
  #print(mean(persistences))
  return(mean(persistences))
}

genPersistenceGrid <- function(filePrefix) {
  surEffects <- seq(0.05,1.0,by=0.05)
  repEffects <- seq(0.00,1.0,by=0.05)
  persisGrid <- matrix(nrow=20,ncol=21,dimnames=list(surEffects,repEffects))
  for (s in surEffects) {
    for (r in repEffects) {
      persisGrid[paste0(s),paste0(r)] <- getAvgPersistence(s,r,filePrefix)
    }
  }
  return(persisGrid)
}

persistenceHeatMap <- function(persisMatrix) {
  heatmap(persisMatrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="none", margins=c(5,10),
          main="Generations of Persistance in N = 1000", xlab="Reproductive Effect",ylab="Survival Effect")
}

avgPersisMatrix <- genPersistenceGrid("../Results/EffectGrid/")
persistenceHeatMap(avgPersisMatrix)
library(MASS)
write.matrix(avgPersisMatrix,file="effectGridAvgPersistenceMatrix.txt",sep=" ")
#save(avgPersisMatrix,file="avgPersisMatrix.Rdata")
#load("avgPersisMatrix.Rdata")

#fit <- lm(mutFreq$Count ~ I(1/(mutFreq$Generation+1)))
#fit <- nls(mutFreq$Count~mutFreq$Generation)
#plot(mutFreq)
#points(mutFreq$Generation, predict(fit), type="l", col="red", lwd=2)
