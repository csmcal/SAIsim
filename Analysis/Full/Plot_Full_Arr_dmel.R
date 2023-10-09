# For analyzing and plotting full SAIsim model data 
#

# Set working directory
setwd("/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Analysis/Full")
# Source the plot functions
source("Plot_Full.R")


#
#
# Load and plot the data
#

genDir = '/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Full/'
speDir = 'dmelN3arrangement/'
resultsDir = paste0(genDir,speDir)



if (!exists("allGenAllArrPopStatsTable")){
  if (file.exists(paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))){
    load(paste0(resultsDir,"allGenAllArrPopStatsTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenAllArrPopStatsTable.Rdata', needs to be generated with Extract_Full_Arr_dmel.R")
  }
}


subsetFinalGenAllArrPopStatsTable <- rbind(filter(allGenAllArrPopStatsTable,
                                                     Generation==max(allGenAllArrPopStatsTable$Generation),
                                                     MutRate==1.042e-4,
                                                     InvRate==1.042e-4,
                                                     EncounterNum==100),
                                              filter(allGenAllArrPopStatsTable,
                                                     Generation==max(allGenAllArrPopStatsTable$Generation)-1,
                                                     LifeHistoryStage!="Zygote",
                                                     MutRate==1.042e-4,
                                                     InvRate==1.042e-4,
                                                     EncounterNum==100))



histogramPanelPlot(filter(subsetFinalGenAllArrPopStatsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Parents"),
                   "Frequency",
                   facet1vars = "ArrangementScenario",
                   facet2vars = "NoMaleCost",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   width = 8,
                   height = 7,
                   outPrefix = "./Arr.")

histogramPanelPlot(filter(subsetFinalGenAllArrPopStatsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==TRUE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Parents"),
                   "FemaleFreq",
                   facet1vars = "ArrangementScenario",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,1200),
                   width = 5,
                   height = 8,
                   outPrefix = "./Arr.")

histogramPanelPlot(filter(subsetFinalGenAllArrPopStatsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==TRUE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Parents"),
                   "MaleFreq",
                   facet1vars = "ArrangementScenario",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,1200),
                   width = 5,
                   height = 8,
                   outPrefix = "./Arr.")

