# SAIsim Mutation Frequency by Spacing Script (with Inversions)
#

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis")


# getMeanLast100Freqs <- function(spacing,fileName) {
#   mutFreq <- read.table(fileName,sep='\t',header=TRUE)
#   # print(colMeans(tail(mutFreq[,-1],n=100)))
#   # meanFinalFreq <- cbind(c(spacing),colMeans(tail(mutFreq[,-1],n=100)))
#   meanFinalFreq <- c(spacing,colMeans(tail(mutFreq[,-1],n=100)))
#   # print(meanFinalFreq)
#   return(meanFinalFreq)
# }
# 
# getFinalFreq <- function(spacing,fileName) {
#   mutFreq <- read.table(fileName,sep='\t',header=TRUE)
#   finalFreq <- c(spacing,tail(mutFreq[,-1],n=1))
#   return(finalFreq)
# }

getFinalFreqMutInvHap <- function(spacing,mutFileName,invFileName,hapFileName) {
  mutFreq <- read.table(mutFileName,sep='\t',header=TRUE)
  invFreq <- read.table(invFileName,sep='\t',header=TRUE)
  hapFreq <- read.table(hapFileName,sep='\t',header=TRUE)
  lastMutFreqs <- unlist(tail(mutFreq[,-1],n=1))
  lastInvFreq <- unlist(tail(invFreq[,2],n=1))
  lastHapFreqs <- unlist(tail(hapFreq[,-1],n=1))
  finalGenData <- c(spacing,lastMutFreqs,lastInvFreq,lastHapFreqs)
  return(finalGenData)
}

getFreqsAtSpa <- function(filePrefix,spacing,numReps) {
  s <- format(spacing,nsmall = 4)
  # print(s)
  fileBase <- paste(filePrefix,"spacing",s,"n", sep="")
  # print(fileBase)
  finalFreqs <- cbind()
  for (num in seq(numReps)-1) {
    n <- formatC(num, width = 3, format = "d", flag = "0")
    try({
      invFileName <- paste(fileBase,n,"Inv.txt",sep="")
      mutFileName <- paste(fileBase,n,"Mut.txt",sep="")
      hapFileName <- paste(fileBase,n,"Hap.txt",sep="")
      finalFreqs <- rbind(finalFreqs,getFinalFreqMutInvHap(spacing,mutFileName,invFileName,hapFileName))
    })
  }
  # print(finalFreqs)
  return(finalFreqs)
}

genFreqSpaTable <- function(filePrefix,spacing,numReps) {
  # spacing <- seq(0.0,0.35,by=0.025)
  # repEffects <- seq(0.00,1.0,by=0.05)
  freqSpaTable <- cbind()
  # summaryTable <- cbind()
  for (s in spacing) {
    freqSpaTable <- rbind(freqSpaTable,getFreqsAtSpa(filePrefix,s,numReps))
    # for m in 2:ncol(freqSpaTable) {
    #   
    # }
  }
  freqSpaTable <- as.data.frame(freqSpaTable)
  colnames(freqSpaTable)[1]<-"Spacing"
  colnames(freqSpaTable)[2]<-"M1"
  colnames(freqSpaTable)[3]<-"M2"
  colnames(freqSpaTable)[4]<-"I1"
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
    mutationData <- cbind(colnames(freqSpaTable)[m],mutationData)
    # print(mutationData)
    colnames(mutationData) <- c("Mutation","Spacing","Count")
    factorizedTable <- rbind(factorizedTable,mutationData)
  }
  # colnames(factorizedTable)[1] <- "Mutation"
  # colnames(factorizedTable)[3] <- "Count"
  factorizedTable <- as.data.frame(factorizedTable)
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

summaryFixLoss <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                           .drop=TRUE) {
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
                   c(N     = length2(xx[[col]], na.rm=na.rm),
                     Nfix  = length2(subset(xx[[col]],xx[[col]]==2000), na.rm=na.rm),
                     Nloss = length2(subset(xx[[col]],xx[[col]]==0), na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac$Pfix <- datac$Nfix/datac$N
  datac$Ploss <- datac$Nloss/datac$N
  datac$Peither <- datac$Pfix + datac$Ploss
  datac$Ppoly <- rep(1.0,length(datac$Peither)) - datac$Peither
  
  return(datac)
}

spacingInvPropFixPlot <- function(freqSpaInvTable) {
  fixLossTable <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                                 measurevar = "I1")
  # print(fixLossTable)
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  p <- ggplot(fixLossTable) + 
    geom_line(aes(x=Spacing, y=Ppoly, group=1), colour="steelblue2") +
    # geom_line(aes(x=Spacing, y=Pfix, group=1), colour="blue") +
    # geom_line(aes(x=Spacing, y=Ploss, group=1), colour="red") +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Percent of Simulations") +
    coord_cartesian(ylim = c(0, 1.0),xlim = c(0, .65)) +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Inv 0.02-0.38", "Mut 2"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
  ggsave("Inversion Polymorphism.png", units="in", width=4, height=3, dpi=300, device = 'png')
}

removeFixLoss <- function(freqSpaInvTable,measurevar){
  temp <- subset(freqSpaInvTable,freqSpaInvTable[[measurevar]]!=2000)
  temp <- subset(temp,temp[[measurevar]]!=0)
  return(temp)
}

spacingInvFreqPlot <- function(freqSpaInvTable) {
  polyTable <- removeFixLoss(freqSpaInvTable,'I1')
  polyTable$I1 <- polyTable$I1/2000
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot(polyTable) + 
    geom_boxplot(aes(x=Spacing, y=I1, group=Spacing),color="yellow3") +
    # geom_line(aes(x=Spacing, y=Peither, group=1), colour="black") +
    # geom_line(aes(x=Spacing, y=Pfix, group=1), colour="blue") +
    # geom_line(aes(x=Spacing, y=Ploss, group=1), colour="red") +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Population Frequency") +
    coord_cartesian(ylim = c(0, 1.0),xlim = c(0, .35)) +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Inv 0.02-0.38", "Mut 2"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
  ggsave("Inversion Frequency.png", units="in", width=4, height=3, dpi=300, device = 'png')
}

spacingInvHapPlot <- function(freqSpaInvTable) {
  fixLossTable <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                                 measurevar = "I1")
  print(fixLossTable)
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot(fixLossTable) + 
    geom_line(aes(x=Spacing, y=Peither, group=1), colour="black") +
    geom_line(aes(x=Spacing, y=Pfix, group=1), colour="blue") +
    geom_line(aes(x=Spacing, y=Ploss, group=1), colour="red") +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Percent of Simulations") +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Inv 0.02-0.38", "Mut 2"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
}

spacingMutPropFixPlot <- function(freqSpaTable) {
  fixLossTable0 <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                  measurevar = "M1")
  fixLossTable1 <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                  measurevar = "M2")
  
  library(ggplot2)
  pd = position_dodge(0.01)
  # ggplot(frqsSE, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  p <- ggplot() + 
    geom_line(data=fixLossTable0,aes(x=Spacing, y=Ppoly, group=1), colour="lightcoral") +
    geom_line(data=fixLossTable1,aes(x=Spacing, y=Ppoly, group=1), colour="red4") +
    # geom_line(aes(x=Spacing, y=Pfix, group=1), colour="blue") +
    # geom_line(aes(x=Spacing, y=Ploss, group=1), colour="red") +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    xlab("Spacing (M)") +
    ylab("Percent of Simulations") +
    coord_cartesian(ylim = c(0, 1.0),xlim = c(0, .65)) +
    scale_colour_hue(name="Inversion",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Inv 0.02-0.38", "Mut 2"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Inversion Fixation and Loss Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.28,0.8))               # Position legend in bottom right
  ggsave("Mutation Polymorphism Inv.png", units="in", width=4, height=3, dpi=300, device = 'png')
}

spacingFreqPlot <- function(freqSpaTable) {
  frqs <- removeFixLoss(freqSpaTable,'M1')
  frqs <- convertToFactors(frqs)
  frqs$Count <- frqs$Count/2000
  # print(head(frqs))
  # frqsSE <- summarySE(frqs,measurevar = "Count",groupvars = c("Mutation","Spacing"))
  # print(frqsSE)
  
  library(ggplot2)
  pd = position_dodge(0)
  # ggplot(frqs, aes(x=Spacing, y=Count, colour=Mutation, group=Mutation)) + 
  ggplot() + 
    # geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), colour="black", width=.02, position=pd) +
    # geom_line(position=pd) +
    # geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
    geom_boxplot(data=subset(frqs,Mutation=='M1'), aes(x=Spacing, y=Count, colour=Mutation, group=Spacing), color="red3") +
    # geom_boxplot(data=subset(subset(frqs,Mutation==1),Count!=0), aes(x=Spacing, y=Count, colour=Mutation, group=Spacing)) +
    geom_boxplot(data=subset(frqs,Mutation=='M2'), aes(x=Spacing, y=Count, colour=Mutation, group=Spacing), color="cyan4") +
    xlab("Spacing (M)") +
    ylab("Count") +
    coord_cartesian(ylim = c(0, 1.0),xlim = c(0, .35)) +
    scale_colour_hue(name="Mutation",    # Legend label, use darker colors
                     breaks=c("1", "2"),
                     labels=c("Sur 0.86 Rep 0.06","Sur 0.70 Rep 0.15"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("Final Average Mutation Frequencies by Spacing") +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(breaks=0:20*4) +         # Set tick every 4
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(0.34,0.01))               # Position legend in bottom right
  ggsave("Mutation Frequency Inv.png", units="in", width=4, height=3, dpi=300, device = 'png')
}

numReps = 500
spacing <- seq(0.0,0.65,by=0.0125)
freqSpaInvTable <- genFreqSpaTable("../Results/TwoSpacing/TwoSpacingInv/",spacing,numReps)
save(freqSpaInvTable,file="SpaInvEqHapFreqTable.Rdata")
load(file="SpaInvEqHapFreqTable.Rdata")
spacingInvPropFixPlot(freqSpaInvTable)
spacingInvFreqPlot(freqSpaInvTable)
spacingMutPropFixPlot(freqSpaInvTable)
spacingFreqPlot(freqSpaInvTable)
# spacingFreqScatter(freqSpaTable)

# testSum <- summarySE(freqSpaInvTable,measurevar = "I1",groupvars = c("Spacing"))

#library(MASS)
#write.matrix(freqGridMatrix,file="effectGridEqFreqMatrix.txt",sep=" ")
#save(freqGridMatrix,file="freqGridMatrix.Rdata")
#load("freqGridMatrix.Rdata")

#fit <- lm(mutFreq$Count ~ I(1/(mutFreq$Generation+1)))
#fit <- nls(mutFreq$Count~mutFreq$Generation)
#plot(mutFreq)
#points(mutFreq$Generation, predict(fit), type="l", col="red", lwd=2)
