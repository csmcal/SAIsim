# SAIsim EffectGrid Heatmap Generator Script
#

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/MutPolyHeatmap/")
dataDirSTD = "../../Results/SingMutGrid/EffectGrid/"
dataDirNMC = "../../Results/SingMutGrid/EffectGridNMC/"
dataDirNFC = "../../Results/SingMutGrid/EffectGridNFC/"
dataDirE10 = "../../Results/SingMutGrid/EffectGridE10/"
dataDirN4 = "../../Results/SingMutGrid/EffectGridN4/"

getSimEquilFreq <- function(fileName,useLastThou) {
  mutFreq <- read.table(fileName,sep='\t',header=TRUE)
  if (useLastThou) {
    avgTailFreq <- mean(tail(mutFreq$Count,n=1000))
    return(avgTailFreq)
  } else {
    finalFreq <- tail(mutFreq$Count,n=1)
    return(finalFreq)
  }
}

getAvgEquil <- function(survEffect,reprEffect,filePrefix,numReps,useLastThou) {
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
      equilFreqs <- c(equilFreqs,getSimEquilFreq(fileName,useLastThou))
    })
  }
  #print(equilFreqs)
  #print(mean(equilFreqs))
  return(mean(equilFreqs))
}

genFreqGrid <- function(filePrefix,surEffects,repEffects,numReps,useLastThou=FALSE) {
  eqFreqGrid <- matrix(nrow=length(surEffects),
                       ncol=length(repEffects),
                       dimnames=list(surEffects,repEffects))
  for (s in surEffects) {
    for (r in repEffects) {
      eqFreqGrid[paste0(s),paste0(r)] <- getAvgEquil(s,r,filePrefix,numReps,useLastThou)
    }
  }
  return(eqFreqGrid)
}

equilHeatMap <- function(eqFreqMatrix) {
  heatmap(eqFreqMatrix, Rowv=NA, Colv=NA, col = heat.colors(256), scale="none", margins=c(5,10),
          main="Equilibrium Allele Counts in N = 1000", xlab="Reproductive Effect",ylab="Survival Effect")
}

surEffects <- seq(0.0,1.0-0.025,by=0.025)
repEffects <- seq(0.025,1.0,by=0.025)
numReps <- 500

effectGridMatrix <- genFreqGrid(dataDirSTD,surEffects,repEffects,numReps)
effectGridNMCMatrix <- genFreqGrid(dataDirNMC,surEffects,repEffects,numReps)
effectGridNFCMatrix <- genFreqGrid(dataDirNFC,surEffects,repEffects,numReps)
effectGridE10Matrix <- genFreqGrid(dataDirE10,surEffects,repEffects,numReps)

equilHeatMap(effectGridMatrix)
#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(freqGridMatrix,file="freqGridMatrix.Rdata")
#load("freqGridMatrix.Rdata")

# Heatmap script
# library(gplots)
library(tidyverse)
library(RColorBrewer)

# mat=matrix(rnorm(100,0,2),nrow=10)

my_palette <- colorRampPalette(c("cornflowerblue","green", "yellow"))(n = 3999)


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
    theme(legend.position="none") +
    scale_fill_gradientn(colours = c("cornflowerblue","green", "yellow"))
  return(p)
}


ggsave("STD Heatmap.png",plot=effGridHeatmap(effectGridMatrix),units="in", width=2.5, height=2.5, dpi=600, device = 'png')

ggsave("NMC Heatmap.png",plot=effGridHeatmap(effectGridNMCMatrix),units="in", width=2.5, height=2.5, dpi=600, device = 'png')

ggsave("NFC Heatmap.png",plot=effGridHeatmap(effectGridNFCMatrix),units="in", width=2.5, height=2.5, dpi=600, device = 'png')

ggsave("E10 Heatmap.png",plot=effGridHeatmap(effectGridE10Matrix),units="in", width=2.5, height=2.5, dpi=600, device = 'png')


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



