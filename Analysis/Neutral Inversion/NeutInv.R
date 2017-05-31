# SAIsim EffectGrid Heatmap Generator Script
#

# Set working directory
setwd("/Users/cmcallester/Documents/Pool Lab/SAIsim/Analysis")


getMeanLastGenInv <- function(power,countFileName,charFileName) {
  invFreq <- read.table(countFileName,sep='\t',header=TRUE)
  invChar <- read.table(charFileName,sep='\t',header=TRUE)
  lastInvFreqs <- unlist(tail(invFreq[,-1],n=1))
  sumExtant <- 0
  sumLength <- 0
  numPoly <- 0
  sumLengthPoly <- 0
  maxLength <- 0
  maxLenCount <- 0
  for (i in seq(length(lastInvFreqs))) {
    numI <- as.integer(lastInvFreqs[i])
    if (numI != 0 && numI != 2000) {
      # print(numI)
      sumExtant <- sumExtant + numI
      invLen <- (invChar[i,2]-invChar[i,1])
      sumLength <- sumLength + numI*invLen
      sumLengthPoly <- sumLengthPoly + invLen
      numPoly <- numPoly + 1
      if (invLen > maxLength) {
        maxLength <- invLen
        maxLenCount <- numI
      }
    }
  }
  # meanFinalCount <- c(power,sum(unlist(tail(invFreq[,-1],n=1)))/2000)
  # print(meanFinalCount)
  loadData <- c(power,sumExtant/2000,sumLength/2000,numPoly,sumLengthPoly/numPoly,maxLength,maxLenCount)
  return(loadData)
}

# getFinalFreq <- function(spacing,fileName) {
#   mutFreq <- read.table(fileName,sep='\t',header=TRUE)
#   finalFreq <- c(spacing,tail(mutFreq[,-1],n=1))
#   return(finalFreq)
# }

getLoadsAtPower <- function(power,filePrefix) {
  p <- format(power,nsmall = 1)
  # print(p)
  fileBase <- paste(filePrefix,p,"n", sep="")
  # print(fileBase)
  finalLoads <- cbind()
  for (num in seq(100)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      countFileName <- paste(fileBase,n,"Count.txt",sep="")
      charFileName <- paste(fileBase,n,"Char.txt",sep="")
      finalLoads <- rbind(finalLoads,getMeanLastGenInv(power,countFileName,charFileName))
    })
  }
  # print(finalLoads)
  return(finalLoads)
}

genLoadPowerTable <- function(filePrefix) {
  power <- seq(-3.0,-8.0,by=-0.5)
  loadPowerTable <- cbind()
  for (p in power) {
    loadPowerTable <- rbind(loadPowerTable,getLoadsAtPower(p,filePrefix))
  }
  # print(loadPowerTable)
  # colnames(loadPowerTable)[1]<-"Power"
  loadPowerTable<-as.data.frame(loadPowerTable)
  colnames(loadPowerTable)<-c("Power","AvgNum","AvgLen","NumPoly","AvgLenPoly","MaxLen","MaxLenCount")
  # print(loadPowerTable)
  return(loadPowerTable)
}


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



# # Gets mutation frequency data in the right categorical format
# convertToFactors <- function(freqSpaTable) {
#   factorizedTable <- cbind()
#   for (m in 2:ncol(freqSpaTable)) {
#     mutationData <- freqSpaTable[,c(1,m)]
#     mutationData <- cbind(m-1,mutationData)
#     # print(mutationData)
#     factorizedTable <- rbind(factorizedTable,mutationData)
#   }
#   colnames(factorizedTable)[1] <- "Mutation"
#   colnames(factorizedTable)[3] <- "Count"
#   factorizedTable <- as.data.frame(factorizedTable)
#   # print(head(factorizedTable))
#   factorizedTable$Mutation <- as.factor(factorizedTable$Mutation)
#   # print(factorizedTable)
#   return(factorizedTable)
# }

powerLoadNumPlot <- function(loadPowerTable) {
  colnames(loadPowerTable)[1] <- "Power"
  colnames(loadPowerTable)[2] <- "AvgNum"
  print(head(loadPowerTable))
  # loadPowerTable$Power <- as.factor(unlist(loadPowerTable$Power))
  loadPowerTable$Power <- as.factor(loadPowerTable$Power)
  # print(head(loadPowerTable))
  
  library(ggplot2)
  ggplot() +
    geom_boxplot(data=loadPowerTable, aes(x=Power, y=AvgNum, group=Power)) +
    xlab("Mutation Rate (Powers of 10)") +
    ylab("Inversions/Chromosome") +
    # scale_colour_hue(name="Mutation",    # Legend label, use darker colors
    #                  breaks=c("1", "2"),
    #                  labels=c("Sur 0.86 Rep 0.06","Sur 0.70 Rep 0.15"),
    #                  l=40) +                    # Use darker colors, lightness=40
    ggtitle("Average Number Inversions per Chromosome") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.34,0.01))               # Position legend in bottom right
}

powerLoadLenPlot <- function(loadPowerTable) {
  colnames(loadPowerTable)[1] <- "Power"
  colnames(loadPowerTable)[2] <- "AvgNum"
  colnames(loadPowerTable)[3] <- "AvgLen"
  loadPowerTable$Power <- as.factor(loadPowerTable$Power)
  
  library(ggplot2)
  pd = position_dodge(0)
  # ggplot(frqs, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot() + 
    geom_boxplot(data=loadPowerTable, aes(x=Power, y=AvgLen, group=Power)) +
    # geom_boxplot(data=subset(frqs,Mutation==2), aes(x=Spacing, y=Count, colour=Mutation, group=Spacing)) +
    xlab("Mutation Rate (Powers of 10)") +
    ylab("Inverted Length/Chromosome") +
    # scale_colour_hue(name="Mutation",    # Legend label, use darker colors
    #                  breaks=c("1", "2"),
    #                  labels=c("Sur 0.86 Rep 0.06","Sur 0.70 Rep 0.15"),
    #                  l=40) +                    # Use darker colors, lightness=40
    ggtitle("Average Inverted Length per Chromosome") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.34,0.01))               # Position legend in bottom right
}

powerScatterPlot <- function(power,loadPowerTable) {
  colnames(loadPowerTable)[1] <- "Power"
  colnames(loadPowerTable)[2] <- "AvgNum"
  colnames(loadPowerTable)[3] <- "AvgLen"
  loadPowerTable$Power <- as.factor(loadPowerTable$Power)
  
  library(ggplot2)
  pd = position_dodge(0)
  # ggplot(loadPowerTable, aes(x=AvgNum, y=AvgLen, colour=Power, group=Power)) +
  # ggplot(subset(loadPowerTable,Power==-5|-4|-3), aes(x=AvgNum, y=AvgLen, colour=Power, group=Power)) +
  ggplot(subset(loadPowerTable,Power==power), aes(x=AvgNum, y=AvgLen, group=Power)) +
    geom_point(size=3, shape=1) + # 21 is filled circle
    geom_smooth(method=lm) +  # Add linear regression line (w/95% default)
    xlab("Average Number/Chromosome") +
    ylab("Average Length/Chromosome") +
    # scale_colour_hue(name="Mutation",    # Legend label, use darker colors
    #                  breaks=c("1", "2"),
    #                  labels=c("Sur 0.86 Rep 0.06","Sur 0.70 Rep 0.15"),
    #                  l=40) +                    # Use darker colors, lightness=40
    ggtitle(paste0("Polymorphic Inversion Averages in 100 Replicates of 10^",format(power,nsmall = 1)," Mut Rate")) +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          plot.title = element_text(size=9.5,face="bold"),
          legend.position=c(0.87,0.01))               # Position legend in bottom right
}

powerMaxScatterPlot <- function(power,loadPowerTable) {
  colnames(loadPowerTable)[1] <- "Power"
  colnames(loadPowerTable)[2] <- "AvgNum"
  colnames(loadPowerTable)[3] <- "AvgLen"
  loadPowerTable$Power <- as.factor(loadPowerTable$Power)
  
  library(ggplot2)
  pd = position_dodge(0)
  # ggplot(loadPowerTable, aes(x=AvgNum, y=AvgLen, colour=Power, group=Power)) +
  # ggplot(subset(loadPowerTable,Power==-5|-4|-3), aes(x=AvgNum, y=AvgLen, colour=Power, group=Power)) +
  ggplot(subset(loadPowerTable,Power==power), aes(x=MaxLen, y=MaxLenCount, group=Power)) +
    geom_point(size=3, shape=1) + # 21 is filled circle
    geom_smooth(method=lm) +  # Add linear regression line (w/95% default)
    xlab("Maximum Lenth Polymorphic Inversion (M)") +
    ylab("Count") +
    # scale_colour_hue(name="Mutation",    # Legend label, use darker colors
    #                  breaks=c("1", "2"),
    #                  labels=c("Sur 0.86 Rep 0.06","Sur 0.70 Rep 0.15"),
    #                  l=40) +                    # Use darker colors, lightness=40
    ggtitle(paste0("Polymorphic Inversion Maximum Length in 100 Replicates of 10^",format(power,nsmall = 1)," Mut Rate")) +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          plot.title = element_text(size=8,face="bold"),
          legend.position=c(0.87,0.01))               # Position legend in bottom right
}

# loadPowerTable <- genLoadPowerTable("../Results/NeutralInversion/iMP")
# save(loadPowerTable,file="neutInvLoads.Rdata")
# powerLoadNumPlot(loadPowerTable)
# powerLoadLenPlot(loadPowerTable)
# powerScatterPlot(-3.5,loadPowerTable)
powerMaxScatterPlot(-3.5,loadPowerTable)

#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(loadPowerTable,file="neutInvLoads.Rdata")
#load("freqGridMatrix.Rdata")

#fit <- lm(mutFreq$Count ~ I(1/(mutFreq$Generation+1)))
#fit <- nls(mutFreq$Count~mutFreq$Generation)
#plot(mutFreq)
#points(mutFreq$Generation, predict(fit), type="l", col="red", lwd=2)
