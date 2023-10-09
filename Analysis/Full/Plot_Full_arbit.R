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

resultsDir = '/Users/cmcallester/cMcAl/UWM/Pool/SAIsim/Results/Full/arbit/'
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





if (!exists("fullAllGenStatTable")){
  if (file.exists(paste0(resultsDir,"fullAllGenStatTable.Rdata"))){
    load(paste0(resultsDir,"fullAllGenStatTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'fullAllGenStatTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

if (!exists("allGenDiffsTable")){
  if (file.exists(paste0(resultsDir,"allGenDiffsTable.Rdata"))){
    load(paste0(resultsDir,"allGenDiffsTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'allGenDiffsTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

scatterPanelPlot(filter(allGenDiffsTable,
                        MutRate==1e-3,
                        ConversionRate==1e-2,
                        EncounterNum==100,
                        Generation %in% seq(0,20000,by=100),
                        ReplicateNum %in% seq(0,2)),
                 "AbsDiffAvgSurEff",
                 "AbsDiffAvgRepEff",
                 colorvar = "Generation",
                 # x.trans = "log10",
                 # y.trans = "log10",
                 x.limits = c(0,1),
                 width = 5,
                 height = 3,
                 outPrefix = "./FullGenScatter.")


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

# testAllGenDiffsTable <- calcDiffsTable(filter(fullAllGenStatTable,
#                                               MutRate==1e-3,
#                                               ConversionRate==1e-2,
#                                               EncounterNum==100))
#                                               # Generation %in% seq(0,20000,by=2000)))
# 
# testAllGenDiffsTable <- filter(allGenAvgDiffsTable,
#                                MutRate==1e-3,
#                                ConversionRate==1e-2,
#                                EncounterNum==100)
# 
# testAllGenAvgDiffsTable <- calcAvgDiffsTable(filter(fullAllGenStatTable,
#                                               MutRate==1e-3,
#                                               ConversionRate==1e-2,
#                                               EncounterNum==100))
#                                               # Generation %in% seq(0,20000,by=2000)))

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



if (!exists("i4e100c2N3InvFreqTable")){
  if (file.exists(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))){
    load(paste0(resultsDir,"i4e100c2N3InvFreqTable.Rdata"))
  } else {
    print("COULD NOT LOAD 'i4e100c2N3InvFreqTable.Rdata', needs to be generated with Extract_Full.R")
  }
}

freqBinPanelPlot(i4e100c2N3InvFreqTable,
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 8,
                 height = 2,
                 outPrefix = "./MutRateControl.")

freqBinPanelPlot(filter(i4e100c2N3InvFreqTable,
                        Generation<=4000),
                 "Generation",
                 "InvFreq",
                 facet2vars = "MutRate",
                 weightvar = "OccurrenceCount",
                 binwidths = c(100,.01),
                 width = 9,
                 height = 3,
                 outPrefix = "./MutRateControl.Short.")


histogramPanelPlot(filter(i4e100c2N3InvFreqTable,
                          Generation==20000),
                   "InvFreq",
                   facet2vars = "MutRate",
                   weightvar = "OccurrenceCount",
                   binwidths = c(.01),
                   width = 8,
                   height = 4,
                   outPrefix = "./MutRateControl.")


if (!exists("m3e100c2N3PopStatTable")){
  if (file.exists(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))){
    load(paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
  } else {
    m3e100c2N3PopStatTable <- extractPopStatTable(N3dir,c(1e-3),c(0,1e-4),c(100),1,
                                            c(1e-2),1e3,2e4,1,1000)
    save(m3e100c2N3PopStatTable,file=paste0(resultsDir,"m3e100c2N3PopStatTable.Rdata"))
  }
}


scatterPanelPlot(m3e100c2N3PopStatTable,
                 "SurVar",
                 "RepVar",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")


scatterPanelPlot(filter(m3e100c2N3PopStatTable, Generation==20000),
                         "SurVar",
                         "RepVar",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")



scatterPanelPlot(m3e100c2N3PopStatTable,
                 "SurMean",
                 "RepMean",
                 facet2vars = "InvRate",
                 colorvar = "Generation",
                 x.trans = "log10",
                 x.coord.trans = "reverse",
                 width = 8,
                 height = 3,
                 outPrefix = "./InvRateControl.")

scatterPanelPlot(filter(m3e100c2N3PopStatTable,
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

scatterPanelPlot(filter(m3e100c2N3PopStatTable, Generation==20000),
                         "SurMean",
                         "RepMean",
                         facet2vars = "InvRate",
                         width = 8,
                         height = 3,
                         outPrefix = "./InvRateControl.")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,4000,by=200)),
                "Generation",
                c("SurMean","RepMean"),
                facet2vars = "InvRate",
                width = 11,
                height = 4,
                outPrefix = "./InvRateControl.")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==0),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate0e+0")

violinPanelPlot(filter(m3e100c2N3PopStatTable, Generation %in% seq(0,10000,by=1000), InvRate==1e-4),
                "Generation",
                c("RepMean","SurMean"),
                width = 16,
                height = 6,
                outPrefix = "./InvRateControl.InvRate1e-4")


repAvgm3e100c2N3PopStatTable <- averagePopStatsByRep(m3e100c2N3PopStatTable)

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
