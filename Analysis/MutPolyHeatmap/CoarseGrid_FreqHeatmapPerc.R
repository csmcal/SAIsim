# SAIsim NoMale EffectGrid Heatmap Generator Script
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
  return(mean(equilFreqs)/2000)
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

# equilHeatMap <- function(eqFreqMatrix) {
#   heatmap(eqFreqMatrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="none", #margins=c(5,10),
#           main="Equilibrium Allele Counts in N = 1000", xlab="Reproductive Effect",ylab="Survival Effect")
# }

freqCoarseGridMatrix <- genEquilFreqGrid("../Results/EffectGrids/EffectGrid/")
save(freqCoarseGridMatrix,file="EqFreqCoarseGrid.Rdata")
# load("EqFreqCoarseGrid.Rdata")

# install.packages( pkgs= c("gplots","RColorBrewer" ) )
# Heatmap script
library(gplots)
library(RColorBrewer)

# mat=matrix(rnorm(100,0,2),nrow=10)

my_palette <- colorRampPalette(c("cornflowerblue","green", "yellow"))(n = 3999)

#custom column breaks
# col_breaks = seq(0,1,length=4000)
# col_breaks = c(seq(-5,-2,length=100),
#                seq(-2,2,length=100), 
#                seq(2,5,length=100))  

# setwd("/Users/jeremylange/Documents/Projects/pigmentation_sims/figures/") #add in the desired file path between quotes(where heat maps will be saved to)

png("CoarseGridEq1.png",    # name for heat map      
    width = 5*900,
    height = 5*900,
    res = 300,
    pointsize = 16)

heatmap.2(freqCoarseGridMatrix, 
          #cellnote = mat,  # label cells with data
          notecex=2,          # size of labels, I think
          main = "Equilibrium Allele Counts in 1000 Diploids", # heat map title
          xlab="Reproductive Effect", ylab="Survival Effect",
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          #margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",     # only draw a row dendrogram
          Rowv="NA",
          #Rowv="NA",
          Colv="NA")            # turn off column clustering

dev.off()

#library(MASS)
#write.matrix(fineGridMatrix,file="fineGridEqFreqMatrix.txt",sep=" ")
#save(fineGridMatrix,file="fineGridMatrix.Rdata")
#load("fineGridMatrix.Rdata")

