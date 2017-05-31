# SAIsim Mutation Frequency by Spacing Script (with Inversions)
#

# Set working directory
setwd("/Users/cmcallester/Documents/Pool Lab/SAIsim/Analysis")


getMeanLast100Freqs <- function(spacing,fileName) {
  mutFreq <- read.table(fileName,sep='\t',header=TRUE)
  # print(colMeans(tail(mutFreq[,-1],n=100)))
  # meanFinalFreq <- cbind(c(spacing),colMeans(tail(mutFreq[,-1],n=100)))
  meanFinalFreq <- c(spacing,colMeans(tail(mutFreq[,-1],n=100)))
  # print(meanFinalFreq)
  return(meanFinalFreq)
}

getFinalFreq <- function(spacing,fileName) {
  mutFreq <- read.table(fileName,sep='\t',header=TRUE)
  finalFreq <- c(spacing,tail(mutFreq[,-1],n=1))
  return(finalFreq)
}

getFreqsAtSpa <- function(spacing,filePrefix) {
  s <- format(spacing,nsmall = 3)
  # print(s)
  fileBase <- paste(filePrefix,"spacing",s,"n", sep="")
  # print(fileBase)
  equilFreqs <- cbind()
  for (num in seq(50)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      fileName <- paste(fileBase,n,"Inv.txt",sep="")
      equilFreqs <- rbind(equilFreqs,getMeanLast100Freqs(spacing,fileName))
    })
  }
  # print(equilFreqs)
  return(equilFreqs)
}

genFreqSpaTable <- function(filePrefix) {
  spacing <- seq(0.0,0.35,by=0.025)
  # repEffects <- seq(0.00,1.0,by=0.05)
  freqSpaTable <- cbind()
  # summaryTable <- cbind()
  for (s in spacing) {
    freqSpaTable <- rbind(freqSpaTable,getFreqsAtSpa(s,filePrefix))
    # for m in 2:ncol(freqSpaTable) {
    #   
    # }
  }
  colnames(freqSpaTable)[1]<-"Spacing"
  # print(freqSpaTable)
  return(freqSpaTable)
}


# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)

# Plotting Data
#----------------------

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



# Gets mutation frequency data in the right categorical format
convertToFactors <- function(freqSpaTable) {
  factorizedTable <- cbind()
  for (m in 2:ncol(freqSpaTable)) {
    mutationData <- freqSpaTable[,c(1,m)]
    mutationData <- cbind(m-1,mutationData)
    # print(mutationData)
    factorizedTable <- rbind(factorizedTable,mutationData)
  }
  colnames(factorizedTable)[1] <- "Mutation"
  colnames(factorizedTable)[3] <- "Count"
  factorizedTable <- as.data.frame(factorizedTable)
  # print(head(factorizedTable))
  factorizedTable$Mutation <- as.factor(factorizedTable$Mutation)
  # print(factorizedTable)
  return(factorizedTable)
}

spacingFreqPlot <- function(freqSpaTable) {
  frqs <- convertToFactors(freqSpaTable)
  frqsSE <- summarySE(frqs,measurevar = "Count",groupvars = c("Mutation","Spacing"))
  # print(frqsSE)
  
  library(ggplot2)
  pd = position_dodge(0)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot(subset(frqsSE,Mutation==1), aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
    geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), colour="black", width=.02, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Count") +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
      breaks=c("1", "2"),
      labels=c("Inv 0.02-0.38", "Mut 2"),
      l=40) +                    # Use darker colors, lightness=40
    ggtitle("Final Average Inversion Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.01))               # Position legend in bottom right
}

spacingFreqScatter <- function(freqSpaTable) {
  frqs <- convertToFactors(freqSpaTable)
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot(subset(frqs,Mutation==1), aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
    # geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), colour="black", width=.02, position=pd) +
    # geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Count") +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
      breaks=c("1", "2"),
      labels=c("Inv 0.02-0.38", "Mut 2"),
      l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Final Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
}

# freqSpaTable <- genFreqSpaTable("../Results/Spa.86.06.7.15e100InvShort/sS0.860rS0.060sB0.700rB0.150")
# sumFreqSpa <- summarySE()
# spacingFreqPlot(freqSpaTable)
spacingFreqScatter(freqSpaTable)

#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(freqGridMatrix,file="freqGridMatrix.Rdata")
#load("freqGridMatrix.Rdata")

#fit <- lm(mutFreq$Count ~ I(1/(mutFreq$Generation+1)))
#fit <- nls(mutFreq$Count~mutFreq$Generation)
#plot(mutFreq)
#points(mutFreq$Generation, predict(fit), type="l", col="red", lwd=2)
