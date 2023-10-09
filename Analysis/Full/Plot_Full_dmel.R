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
speDir = 'dmelN3/'
resultsDir = paste0(genDir,speDir)
# N3dir = paste0(resultsDir,'N3/')
# N4dir = paste0(resultsDir,'N4/')
# N2dir = paste0(resultsDir,'N2/')


# if (!exists("fullFinalGenStatTable")){
#   if (file.exists(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))){
#     load(paste0(resultsDir,"fullFinalGenStatTable.Rdata"))
#   } else {
#     print("COULD NOT LOAD 'fullFinalGenStatTable.Rdata', needs to be generated with Extract_Full.R")
#   }
# }

# if (!exists("finalGenDiffsTable")){
#   if (file.exists(paste0(resultsDir,"finalGenDiffsTable.Rdata"))){
#     load(paste0(resultsDir,"finalGenDiffsTable.Rdata"))
#   } else {
#     print("COULD NOT LOAD 'finalGenDiffsTable.Rdata', needs to be generated with Extract_Full.R")
#   }
# }
# if (!exists("finalGenAvgDiffsTable")){
#   if (file.exists(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))){
#     load(paste0(resultsDir,"finalGenAvgDiffsTable.Rdata"))
#   } else {
#     print("COULD NOT LOAD 'finalGenAvgDiffsTable.Rdata', needs to be generated with Extract_Full.R")
#   }
# }
# 
# 
# scatterPlotFinalGenParams(finalGenDiffsTable,
#                           "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
#                           colorvar="MaxArrangementFrequency")
# 
# scatterPanelPlot(finalGenDiffsTable,
#                  "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
#                  c("PopSize","EncounterNum"),
#                  c("MutRate","ConversionRate"),
#                  colorvar="MaxArrangementFrequency",
#                  width = 10,
#                  height = 5)
# 
# 
# scatterPanelPlot(filter(finalGenDiffsTable,
#                         ConversionRate == 1e-2,
#                         PopSize == 1e3),
#                  "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
#                  c("EncounterNum"),
#                  c("MutRate"),
#                  colorvar="MaxArrangementFrequency",
#                  width = 10,
#                  height = 5)





# if (!exists("allGenInvStatTable")){
#   if (file.exists(paste0(resultsDir,"allGenInvStatTable.Rdata"))){
#     load(paste0(resultsDir,"allGenInvStatTable.Rdata"))
#   } else {
#     print("COULD NOT LOAD 'allGenInvStatTable.Rdata', needs to be generated with Extract_Full.R")
#   }
# }

if (!exists("allGenDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenDiffsTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenDiffsTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

scatterPanelPlot(filter(allGenDiffsTable),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 facet1vars = c("MutRate","InvRate","CrossoverMapLength"),
                 facet2vars = c("EncounterNum","NoMaleCost","NoFemaleCost"),
                 colorvar = "Generation",
                 # x.trans = "log10",
                 y.trans = "log10",
                 x.limits = c(0,1),
                 width = 20,
                 height = 20,
                 outPrefix = "./FullGenScatter.")

scatterPanelPlot(filter(allGenDiffsTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        EncounterNum%in%c(10,100),
                        MutRate==0.0001042,
                        InvRate%in%c(1.042e-4,1.042e-6)),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 facet1vars = c("MutRate","InvRate","CrossoverMapLength"),
                 facet2vars = c("EncounterNum"),
                 colorvar = "Generation",
                 # x.trans = "log10",
                 # y.trans = "log10",
                 x.limits = c(0,1),
                 width = 10,
                 height = 10,
                 outPrefix = "./FullGenScatter.")

scatterPanelPlot(filter(allGenDiffsTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        EncounterNum==100,
                        MutRate==0.0001042,
                        InvRate==1.042e-4,
                        ReplicateNum%in%seq(0,99),
                        Generation%in%seq(0,20000,by=1000)),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 facet1vars = c("MutRate","InvRate","CrossoverMapLength"),
                 facet2vars = c("EncounterNum"),
                 colorvar = "Generation",
                 # x.trans = "log10",
                 # y.trans = "log10",
                 x.limits = c(0,1),
                 width = 8,
                 height = 5,
                 outPrefix = "./Fig1AFullGenScatter.")

scatterAnimation(filter(allGenDiffsTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        # ReplicateNum%in%seq(0,99),
                        # Generation%in%seq(0,20000,by=1000),
                        CrossoverMapLength==43.6,
                        EncounterNum==100,
                        MutRate==0.0001042,
                        InvRate==1.042e-4),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 "Generation",
                 # facet1vars = c("MutRate","InvRate","EncounterNum"),
                 # facet2vars = c("CrossoverMapLength"),
                 # facet2rev = TRUE,
                 # facet2.labs = list(setNames(c(0.436,43.6),c("Low Recombination","High Recombination"))),
                 # colorvar = "Generation",
                 # x.trans = "log10",
                 # y.trans = "log10",
                 x.limits = c(0,1),
                 y.limits = c(0,8),
                 width = 5,
                 height = 3,
                 outPrefix = "./AnimHRScatter.")

scatterAnimation(filter(allGenDiffsTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        # ReplicateNum%in%seq(0,99),
                        # Generation%in%seq(0,20000,by=1000),
                        CrossoverMapLength==.436,
                        EncounterNum==100,
                        MutRate==0.0001042,
                        InvRate==1.042e-4),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 "Generation",
                 # facet1vars = c("MutRate","InvRate","EncounterNum"),
                 # facet2vars = c("CrossoverMapLength"),
                 # facet2rev = TRUE,
                 # facet2.labs = list(setNames(c(0.436,43.6),c("Low Recombination","High Recombination"))),
                 # colorvar = "Generation",
                 # x.trans = "log10",
                 # y.trans = "log10",
                 x.limits = c(0,1),
                 y.limits = c(0,8),
                 width = 5,
                 height = 3,
                 outPrefix = "./AnimLRScatter.")

if (!exists("allGenAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenAvgDiffsTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenAvgDiffsTable.Rdata', needs to be generated with Extract_Full.R")
  }
}


scatterPanelPlot(filter(allGenAvgDiffsTable,
                        MutRate==1e-3,
                        ConversionRate==1e-2,
                        EncounterNum==100),
                 "AvgAbsDiffAvgSurEff",
                 "AvgAbsDiffAvgRepEff",
                 colorvar = "Generation",
                 x.trans = "log10",
                 width = 5,
                 height = 3,
                 outPrefix = "./FullGenScatter.")

scatterPanelPlot(filter(allGenAvgDiffsTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        EncounterNum==100,
                        MutRate==0.0001042,
                        InvRate==1.042e-4),
                 "AvgAbsDiffAvgRepEff",
                 "AvgAbsDiffAvgSurEff",
                 # facet1vars = c("MutRate","InvRate","CrossoverMapLength"),
                 # facet2vars = c("EncounterNum"),
                 facet2vars = c("CrossoverMapLength"),
                 colorvar = "Generation",
                 # x.trans = "log10",
                 # x.coord.trans = "log10",
                 # y.trans = "log10",
                 # x.limits = c(0,1),
                 y.limits = c(0,0.6),
                 width = 6,
                 height = 3*6/8,
                 # width = 6,
                 # height = 2.8,
                 outPrefix = "./Fig1BAvgDiffs.")


# scatterPlotFinalGenParams(testAllGenDiffsTable,
#                           "AbsDiffAvgSurEff","AbsDiffAvgRepEff",
#                           colorvar="Generation")
# 
# scatterPlotFinalGenParams(testAllGenAvgDiffsTable,
#                           "AvgAbsDiffAvgSurEff","AvgAbsDiffAvgRepEff",
#                           colorvar="Generation")
# 
# violinPanelPlot(filter(testAllGenDiffsTable,Generation %in% seq(0,20000,by=1000)),
#                 "Generation",
#                 c("AbsDiffAvgSurEff","AbsDiffAvgRepEff"),
#                 width = 11,
#                 height = 3,
#                 outPrefix = "./m3c2N3")



if (!exists("allGenInvFreqTable")){
  if (file.exists(paste0(resultsDir,"allGenInvFreqTable.Rdata"))){
    load(paste0(resultsDir,"allGenInvFreqTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenInvFreqTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

freqBinPanelPlot(allGenInvFreqTable,
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 8,
                 height = 2,
                 outPrefix = "./MutRateControl.")

freqBinPanelPlot(filter(allGenInvFreqTable,
                        Generation<=4000),
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 9,
                 height = 3,
                 outPrefix = "./MutRateControl.Short.")


histogramPanelPlot(filter(allGenInvFreqTable,
                          Generation==20000),
                   "InvFreq",
                   facet2vars = "MutRate",
                   weightvar = "OccurrenceCount",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./MutRateControl.")

histogramPanelPlot(filter(allGenInvFreqTable,
                          Generation==20000,
                          NoMaleCost==FALSE,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate%in%c(0,0.0001042),
                          InvRate==1.042e-4),
                   "Freq",
                   facet2vars = "MutRate",
                   weightvar = "OccurrenceCount",
                   colorvar = "FreqRankIncAnc",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./Fig2AMutRateControl.")

histogramPanelPlot(filter(allGenInvFreqTable,
                          Generation==20000,
                          NoMaleCost==FALSE,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate%in%c(0,0.0001042),
                          InvRate==1.042e-4,
                          Ancestral==FALSE),
                   "Freq",
                   facet2vars = "MutRate",
                   weightvar = "OccurrenceCount",
                   # colorvar = "FreqRank",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./Fig2AMutRateControl.")


if (!exists("allGenPopStatTable")){
  if (file.exists(paste0(resultsDir,"allGenPopStatTable.Rdata"))){
    load(paste0(resultsDir,"allGenPopStatTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenPopStatTable.Rdata', needs to be generated with Extract_Full.R")
  }
}


scatterPanelPlot(allGenPopStatTable,
                 "SurVar",
                 "RepVar",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")


scatterPanelPlot(filter(allGenPopStatTable, Generation==20000),
                         "SurVar",
                         "RepVar",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")



scatterPanelPlot(allGenPopStatTable,
                 "SurMean",
                 "RepMean",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")

scatterPanelPlot(filter(allGenPopStatTable,
                        Generation %in% seq(0,20000,by=500),
                        ReplicateNum %in% seq(0,100)),
                 "SurMean",
                 "RepMean",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")

scatterPanelPlot(filter(allGenPopStatTable, Generation==20000),
                         "SurMean",
                         "RepMean",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")

violinPanelPlot(filter(allGenPopStatTable, Generation %in% seq(0,4000,by=200)),
                "Generation",
                c("SurMean","RepMean"),
                facet2vars = "InvRate",
                width = 11,
                height = 4,
                outPrefix = "./InvRateControl.")

violinPanelPlot(filter(allGenPopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==0),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate0e+0")

violinPanelPlot(filter(allGenPopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==1e-4),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate1e-4")


if (!exists("allGenAvgPopStatTable")){
  if (file.exists(paste0(resultsDir,"allGenAvgPopStatTable.Rdata"))){
    load(paste0(resultsDir,"allGenAvgPopStatTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenAvgPopStatTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

scatterPanelPlot(filter(allGenAvgPopStatTable,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE,
                        EncounterNum==100,
                        MutRate==0.0001042,
                        CrossoverMapLength==0.436,
                        LifeHistoryStage=="Adult",
                        InvRate %in% c(0,0.0001042)),
                 "AvgSurMean",
                 "AvgRepMean",
                 # facet1vars = c("CrossoverMapLength"),
                 facet2vars = c("InvRate"),
                 # facet2vars = "InvRate",
                 colorvar = "Generation",
                 # x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 6,
                 height = 3*6/8,
                 outPrefix = "./Fig1AInvRateControl.")

# repAvgm3e100c2N3PopStatTable <- averagePopStatsByRep(m3e100c2N3PopStatTable)

scatterPanelPlot(repAvgm3e100c2N3PopStatTable,
                 "AvgSurMean",
                 "AvgRepMean",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")


if (!exists("finalGenAllInvPopAvgDiffsTable")){
  if (file.exists(paste0(resultsDir,"finalGenAllInvPopAvgDiffsTable.Rdata"))){
    load(paste0(resultsDir,"finalGenAllInvPopAvgDiffsTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'finalGenAllInvPopAvgDiffsTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

# finalGenAllInvPopAvgDiffsTable <- filter(allGenAllInvPopAvgDiffsTable,
#                                          Generation==max(allGenAllInvPopAvgDiffsTable$Generation))


subset1FinalGenAllInvPopAvgDiffsTable <- filter(finalGenAllInvPopAvgDiffsTable,
                                               MutRate==1.042e-4,
                                               InvRate==1.042e-4,
                                               CrossoverMapLength==4.36e-1,
                                               NoMaleCost==FALSE,
                                               NoFemaleCost==FALSE,
                                               EncounterNum==100)

subset2FinalGenAllInvPopAvgDiffsTable <- filter(finalGenAllInvPopAvgDiffsTable,
                                                MutRate==1.042e-2,
                                                InvRate==1.042e-2,
                                                CrossoverMapLength==4.36e1,
                                                NoMaleCost==FALSE,
                                                NoFemaleCost==FALSE,
                                                EncounterNum==100)

subsetFinalGenAllInvPopAvgDiffsTable <- rbind(filter(finalGenAllInvPopAvgDiffsTable,
                                                     Generation==max(finalGenAllInvPopAvgDiffsTable$Generation),
                                                     MutRate==1.042e-4,
                                                     InvRate==1.042e-4,
                                                     EncounterNum==100),
                                              filter(finalGenAllInvPopAvgDiffsTable,
                                                     Generation==max(finalGenAllInvPopAvgDiffsTable$Generation)-1,
                                                     LifeHistoryStage=="Adult",
                                                     MutRate==1.042e-4,
                                                     InvRate==1.042e-4,
                                                     EncounterNum==100))

# subsetFinalGenAllInvPopAvgDiffsTable <- classifyArrangementState(subsetFinalGenAllInvPopAvgDiffsTable)
# save(subsetFinalGenAllInvPopAvgDiffsTable,file=paste0(resultsDir,"subsetFinalGenAllInvPopAvgDiffsTable.Rdata"))
# load(paste0(resultsDir,"subsetFinalGenAllInvPopAvgDiffsTable.Rdata"))

# Trying to assess if HW-equilibrium is a way to assess if a population follows this?
# cutoff for common inversion 
# Expectation should be x axis, observation should be the Y axis
# Plot for the heterozygote instead? Triangle plots?
# Plot for the alt homozygote? (zero in on that space? don't use a solid line for the x axis so that points are visible, don't need log scale)
# distill to two most common karyotypes in the pops, normalize the frequency to the new total, separate by hom. for the most freq karyotype, the less freq karyotype
# Only show females b/c a discrete pool of survivors, unlike males

scatterPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                        MutRate==1.042e-4,
                        InvRate==1.042e-4,
                        EncounterNum==100,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE),
                 "HWExpectedHomozygoteFreq",
                 "HomozygoteFreq",
                 facet1vars = c("CrossoverMapLength"),
                 facet2vars = c("LifeHistoryStage","NoMaleCost"),
                 seg.coords = c(0,0,1,1),
                 width = 20,
                 height = 5,
                 outPrefix = "./Fig2BHomozygoteFrequency.")

scatterPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                        MutRate==1.042e-4,
                        InvRate==1.042e-4,
                        EncounterNum==100,
                        NoMaleCost==FALSE,
                        NoFemaleCost==FALSE),
                 "HWExpectedHeterozygoteFreq",
                 "HeterozygoteFreq",
                 facet1vars = c("CrossoverMapLength"),
                 facet2vars = c("LifeHistoryStage"),
                 seg.coords = c(0,0,.5,.5),
                 width = 6,
                 height = 6,
                 outPrefix = "./Fig2BHeterozygoteFrequency.")

scatterPanelMultiPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                             EncounterNum==100,
                             NoFemaleCost==FALSE),
                      c("HomozygoteFreq","FemaleHomozygoteFreq","MaleHomozygoteFreq"),
                      c("HWExpectedHomozygoteFreq","HWExpectedFemaleHomozygoteFreq","HWExpectedMaleHomozygoteFreq"),
                      facet1vars = c("MutRate","InvRate","CrossoverMapLength"),
                      facet2vars = c("LifeHistoryStage","NoMaleCost"),
                      colors = c("black","red","blue"),
                      seg.coords = c(0,0,1,1),
                      width = 20,
                      height = 10,
                      outPrefix = "./Fig2BHomozygoteFrequency.")


scatterPanelMultiPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                             MutRate==1.042e-4,
                             InvRate==1.042e-4,
                             EncounterNum==100,
                             CrossoverMapLength==0.436,
                             FrequencyRank==2,
                             NoMaleCost==FALSE,
                             NoFemaleCost==FALSE),
                      c("HWExpectedFemaleHeterozygoteFreq","HWExpectedFemaleHomozygoteFreq"),
                      c("FemaleHeterozygoteFreq","FemaleHomozygoteFreq"),
                      axis1lab="Hardy-Weinberg Expectation",
                      axis2lab="Observed Frequency of the Heterozygote", 
                      # facet1vars = c("InvRate","EncounterNum"),
                      facet2vars = c("LifeHistoryStage"),
                      facet2rev = TRUE,
                      colors = c("yellow3","purple2"),
                      # colors = hcl.colors(2, palette = "Plasma"),
                      seg.coords = c(0,0,.5,.5),
                      width = 5,
                      height = 4,
                      outPrefix = "./Fig2CHeterozygoteFrequency.")

scatterPanelMultiPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                             MutRate==1.042e-4,
                             InvRate==1.042e-4,
                             EncounterNum==100,
                             CrossoverMapLength==0.436,
                             FrequencyRank==2,
                             NoMaleCost==FALSE,
                             NoFemaleCost==FALSE),
                      c("HWExpectedFemaleHeterozygoteFreq","HWExpectedFemaleHomozygoteFreq"),
                      c("FemaleHeterozygoteFreq","FemaleHomozygoteFreq"),
                      axis1lab="Hardy-Weinberg Expectation",
                      axis2lab="Observed Genotype Frequency", 
                      # facet1vars = c("InvRate","EncounterNum"),
                      facet2vars = c("LifeHistoryStage"),
                      facet2rev = TRUE,
                      colors = c("orange3","purple2"),
                      # colors = hcl.colors(2, palette = "Plasma"),
                      seg.coords = c(0,0,.5,.5),
                      width = 8,
                      height = 3.5,
                      outPrefix = "./Fig2CHeterozygoteFrequency.")


histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult",
                          (Frequency > 0.01) & (Frequency < 0.99)),
                   "Frequency",
                   facet2vars = "NoMaleCost",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./Fig2DSurvivalRanked.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "MaleFreq",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   width = 5,
                   height = 2.5,
                   outPrefix = "./Fig4ASharedCosts.")
histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "FemaleFreq",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   width = 5,
                   height = 2.5,
                   outPrefix = "./Fig4ASharedCosts.")


histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "Frequency",
                   facet1vars = "ArrangementScenario",
                   facet2vars = "NoMaleCost",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   width = 8,
                   height = 7,
                   outPrefix = "./Fig2DScenarioSepSurvivalRanked.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==TRUE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "FemaleFreq",
                   facet1vars = "ArrangementScenario",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,1000),
                   width = 5,
                   height = 5.5,
                   outPrefix = "./Fig2DScenarioSepSurvivalRankedFemale.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          NoMaleCost==TRUE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "MaleFreq",
                   facet1vars = "ArrangementScenario",
                   colorvar = "NormalizedSurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,1000),
                   width = 5,
                   height = 5.5,
                   outPrefix = "./Fig2DScenarioSepSurvivalRankedMale.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "FemaleFreq",
                   facet2vars = "NoMaleCost",
                   colorvar = "SurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,800),
                   width = 10,
                   height = 3,
                   outPrefix = "./Fig2DSurvivalRankedFemale.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "MaleFreq",
                   facet2vars = "NoMaleCost",
                   colorvar = "SurvivalRank",
                   binwidths = c(.01),
                   y.limits = c(0,800),
                   width = 10,
                   height = 3,
                   outPrefix = "./Fig2DSurvivalRankedMale.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "MaleFreq",
                   facet2vars = "NoMaleCost",
                   colorvar = "SurvivalRank",
                   weightvar = "MaleHomozygoteFreq",
                   binwidths = c(.01),
                   y.limits = c(0,800),
                   width = 10,
                   height = 3,
                   outPrefix = "./Fig2DSurvivalRankedMaleHom.")

histogramPanelPlot(filter(subsetFinalGenAllInvPopAvgDiffsTable,
                          NoFemaleCost==FALSE,
                          EncounterNum==100,
                          MutRate==1.042e-4,
                          InvRate==1.042e-4,
                          LifeHistoryStage=="Adult"),
                   "MaleFreq",
                   facet2vars = "NoMaleCost",
                   colorvar = "SurvivalRank",
                   weightvar = "MaleHeterozygoteFreq",
                   binwidths = c(.01),
                   y.limits = c(0,800),
                   width = 10,
                   height = 3,
                   outPrefix = "./Fig2DSurvivalRankedMaleHet.")

