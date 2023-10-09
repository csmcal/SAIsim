# SAIsim Mutation Frequency by Spacing Script (with and w/o Inversions)
# For generating the plot for Prelim B Proposal

# Set working directory
setwd("/Users/cmcallester/cMcAl/Pool Lab/SAIsim/Analysis")


# To retrieve polymorphism data from
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
 
# Plot generator for the figure in Prelim B Proposal
plotPolyPrelimB <- function(freqSpaTable,freqSpaInvTable){
  
  #Combining the data to one table
  polyMutNoInv <- summaryFixLoss(freqSpaTable,groupvars = "Spacing",
                                 measurevar = "M1")
  polyMutInv <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                               measurevar = "M1")
  polyInv <- summaryFixLoss(freqSpaInvTable,groupvars = "Spacing",
                            measurevar = "I1")
  
  # polyMutNoInv <- cbind(polyMutNoInv,Locus=c("Smaller SA Locus with Inversion"))
  # polyMutInv <- cbind(polyMutInv,Locus=c("Smaller SA Locus without Inversion"))
  # polyInv <- cbind(polyInv,Locus=c("Inversion with Breakpoints 0.03 to 0.37"))
  
  polyMutNoInv <- cbind(polyMutNoInv,Locus=c("M"))
  polyMutInv <- cbind(polyMutInv,Locus=c("MI"))
  polyInv <- cbind(polyInv,Locus=c("I"))
  
  # polyMutNoInv <- cbind(Spacing=polyMutNoInv$Spacing,Ppoly=polyMutNoInv$Ppoly,
  #                       Locus=c("M"))
  # polyMutInv <- cbind(Spacing=polyMutInv$Spacing,Ppoly=polyMutInv$Ppoly,
  #                     Locus=c("MI"))
  # polyInv <- cbind(Spacing=polyInv$Spacing,Ppoly=polyInv$Ppoly,
  #                  Locus=c("I"))
  
  polyData <- as.data.frame(rbind(polyMutNoInv,polyMutInv,polyInv))[,c("Spacing","Ppoly","Locus")]
  # print(polyData)
  
  # Plotting using ggplot2
  library(ggplot2)
  p <- ggplot() + 
    geom_line(data=polyData,aes(x=Spacing, y=Ppoly, group=Locus, linetype=Locus)) +
    xlab("Distance Between SA Locuss (M)") +
    ylab("Proportion of Simulations") + 
    scale_x_continuous(breaks=seq(0, 0.35, .05), minor_breaks=seq(0, 0.35, .025)) + 
    scale_y_continuous(breaks=seq(0, 1, .25), minor_breaks=seq(0, 1, .05)) + 
    # coord_cartesian(ylim=c(0,1.0), xlim=c(0,.35)) +
    scale_linetype_manual(values = c(3,4,1),
                     labels=c("Smaller SA Locus (w/o Inv)",
                              "Smaller SA Locus (w/ Inv)",
                              "Inversion (Breakpoints -0.02, 0.32)")) +
    ggtitle("Persistence of Polymoprhism in SAI Simulations") +
    theme_bw() +
    theme(legend.justification=c(1,0),
      legend.position=c(.38,0.13),  # Position legend in bottom right
      legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),
      legend.title = element_text(size=10),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=1.9, size=17, face="bold")) +
  ggsave("Small Effect Locus & Inversion Poly Freq.png", units="in", width=6, height=4, dpi=600, device = 'png')
}

# Run the plot script
load(file="twoMutSpaTable.Rdata")
load(file="SpaInvEqHapFreqTable.Rdata")
plotPolyPrelimB(freqSpaTable,freqSpaInvTable)



