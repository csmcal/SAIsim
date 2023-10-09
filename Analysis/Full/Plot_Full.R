# For analyzing and plotting full SAIsim model data 
#

library(tidyverse)
library(svglite)
library(gganimate)



axisLabels <- c(
  "AvgAbsDiffAvgNumMut"="Abs Val Difference in Number of Mutants",
  # "AvgAbsDiffAvgSurEff"="Abs Val Difference in Survival Effect",
  # "AvgAbsDiffAvgRepEff"="Abs Val Difference in Reproductive Effect",
  "AvgAbsDiffAvgSurEff"="Difference in Survival Effect",
  "AvgAbsDiffAvgRepEff"="Difference in Reproductive Effect",
  "AbsDiffAvgNumMut"="Abs Val Difference in Number of Mutants",
  "AbsDiffAvgSurEff"="Abs Val Difference in Survival Effect",
  "AbsDiffAvgRepEff"="Abs Val Difference in Reproductive Effect",
  # "AvgAbsDiffAvgNumMut"="Magnitude of the Difference in Number of Mutants",
  # "AvgAbsDiffAvgSurEff"="Magnitude of the Difference in Total Survival Effect",
  # "AvgAbsDiffAvgRepEff"="Magnitude of the Difference in Total Reproductive Effect",
  # "AbsDiffAvgNumMut"="Magnitude of the Difference in Number of Mutants",
  # "AbsDiffAvgSurEff"="Magnitude of the Difference in Total Survival Effect",
  # "AbsDiffAvgRepEff"="Magnitude of the Difference in Total Reproductive Effect",
  "AvgNumMut"="Average Number of Mutations in the Arrangement",
  "AvgSurEff"="Average Survival Effect of the Arrangement",
  "AvgRepEff"="Average Reproductive Effect of the Arrangement",
  "AbsDiffCount"="Difference in Count between the Most Dominant Arrangement and all Others",
  "AvgAbsDiffCount"="Average Difference in Count between the Most Dominant Arrangement and all Others",
  "MaxArrangementCount"="Population Count of the Most Common Arrangement",
  "AvgMaxArrangementCount"="Average Population Count of the Most Common Arrangement",
  "AbsDiffFrequency"="Difference in Frequency between the Most Dominant Arrangement and all Others",
  "AvgAbsDiffFrequency"="Average Difference in Frequency between the Most Dominant Arrangement and all Others",
  "MaxArrangementFrequency"="Population Frequency of the Most Common Arrangement",
  "AvgMaxArrangementFrequency"="Average Population Frequency of the Most Common Arrangement",
  "Generation"="Simulated Generation",
  "InvCount"="Inversion Arrangement Count within the Population",
  "InvFreq"="Frequency of the Inversion Arrangement within the Population",
  "Freq"="Frequency of the Arrangement within its Population",
  "FreqRankIncAnc"="Frequency Rank of each Arrangement within its Population",
  "FreqRank"="Frequency Rank of each Arrangement within its Population",
  # "InvCount-Generation"="Number of Occurrences across all Simulations",
  # "Generation-InvCount"="Number of Occurrences across all Simulations",
  # "InvFreq-Generation"="Number of Occurrences across all Simulations",
  # "Generation-InvFreq"="Number of Occurrences across all Simulations",
  "InvCount-Generation"="Count",
  "Generation-InvCount"="Count",
  "InvFreq-Generation"="Count",
  "Generation-InvFreq"="Count",
  "SurMean"="Mean Survival Value",
  "AvgSurMean"="Mean Survival Value",
  # "AvgSurMean"="Mean Survival Value across Simulations",
  "SurVar"="Variance in Survival Value",
  "RepMean"="Mean Display Value",
  "AvgRepMean"="Mean Display Value",
  # "AvgRepMean"="Mean Display Value across Simulations",
  "RepVar"="Variance in Display Value",
  "HomozygoteFreq"="Observed Frequency of the Homozygote",
  "HWExpectedHomozygoteFreq"="Hardy-Weinberg Expectation",
  "HeterozygoteFreq"="Observed Frequency of the Heterozygote",
  "HWExpectedHeterozygoteFreq"="Hardy-Weinberg Expectation"
  )



savePlot <- function(p,saveFile,width = 6,height = 5) {
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}

genSaveFileName <- function(plotTypeStr,
                            axis1var, 
                            axis2var,
                            facet1vars=NULL,
                            facet2vars=NULL,
                            colorvar=NULL,
                            shapevar=NULL,
                            x.trans="identity",
                            y.trans="identity",
                            x.coord.trans="identity",
                            y.coord.trans="identity",
                            x.limits=NULL,
                            y.limits=NULL,
                            seg.coords=NULL,
                            outPrefix = "./"){
  saveFile <- paste0(outPrefix,plotTypeStr,'.',axis1var,'.',axis2var)
  if (!is.null(colorvar)){
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  if (!is.null(shapevar)){
    saveFile <- paste0(saveFile,'.Shaped.',shapevar)
  }
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    saveFile <- paste0(saveFile,'.Panelled.')
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
      saveFile <- paste0(saveFile,facet1str)
      if(!is.null(facet2vars)){
        saveFile <- paste0(saveFile,'~')
      }
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
      saveFile <- paste0(saveFile,facet2str)
    }
  }
  return(saveFile)
}


plotCanvas <- function(simData,
                       axis1lab,
                       axis2lab,
                       facet1vars=NULL,
                       facet2vars=NULL,
                       facet1rev=FALSE,
                       facet2rev=FALSE,
                       facet1.labs=NULL,
                       facet2.labs=NULL,
                       colorvar=NULL,
                       shapevar=NULL,
                       x.trans="identity",
                       y.trans="identity",
                       x.coord.trans="identity",
                       y.coord.trans="identity",
                       x.limits=NULL,
                       y.limits=NULL,
                       seg.coords=NULL,
                       height=NULL,
                       width=NULL) {
  
  p <- ggplot(simData) +
    xlab(axis1lab) +
    ylab(axis2lab) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      if (facet1rev) {
        facet1str <- paste0('fct_rev(',paste(facet1vars,collapse = ')+fct_rev('),')')
      }
      else{
        facet1str <- paste(facet1vars,collapse = '+')
      }
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      if (facet2rev) {
        facet2str <- paste0('fct_rev(',paste(facet2vars,collapse = ')+fct_rev('),')')
      }
      else{
        facet2str <- paste(facet2vars,collapse = '+')
      }
      if (facet2rev) {
        facet2str
      }
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    print(facetForm)
    
    if (!is.null(facet1.labs) | !is.null(facet2.labs)){
      labeller.str <- "labeller("
      if (!is.null(facet1.labs)){
        for (i in seq(length(facet1.labs))){
          labeller.str <- paste0(labeller.str,facet1vars[[i]],"=facet1.labs[[i]],")
        }
      }
      if (!is.null(facet2.labs)){
        for (i in seq(length(facet2.labs))){
          labeller.str <- paste0(labeller.str,facet2vars[[i]],"=facet2.labs[[i]],")
        }
      }
      labeller.str <- paste0(substr(labeller.str,1,nchar(labeller.str)-1),")")
      print(labeller.str)
      print(parse(text=labeller.str))
      
      p <- p + facet_grid(as.formula(facetForm),labeller = parse(text=labeller.str))
    } else {
      p <- p + facet_grid(as.formula(facetForm))
    }
  }
  
  if (!is.null(colorvar)){
    print("colorvar is not null")
    print(paste0(axisLabels[colorvar]))
    maxColorVarVal <- max(select(simData,colorvar))
    minColorVarVal <- min(select(simData,colorvar))
    # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
    
    p <- p +
      scale_colour_gradientn(
        name=axisLabels[colorvar],
        # name=NULL,
        limits = c(minColorVarVal,maxColorVarVal),
        colors = hcl.colors(5, palette = "Plasma"),
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      # scale_colour_gradient(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "#4B0055",
      #   high = "#FDE333",
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) +
      # scale_colour_gradient2(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "red",
      #   mid = "white",
      #   high = "blue",
      #   midpoint = midColorVarVal,
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) +
      guides(
        # size = "none",
        colour = guide_colourbar(title.position = "right")
      )+
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            # legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  
  if (!is.null(seg.coords)) {
    # p <- p + geom_segment(0,0,1,1)
    # p <- p + geom_segment(seg.coords[1],seg.coords[2],seg.coords[3],seg.coords[4])
    # p <- p + geom_segment(seg.coords[1],seg.coords[2],seg.coords[3],seg.coords[4],
    #                       mapping=aes(linetype="dashed"))
    p <- p + geom_segment(aes(x = seg.coords[1],
                              y = seg.coords[2],
                              xend = seg.coords[3],
                              yend = seg.coords[4]),
                          alpha=.01,
                          color="gray25",
                          # color="#40404009",
                          # scale_color_manual("#404040"),
                          linetype="dashed")
  }
  
  if (!is.null(x.limits)) {
    p <- p + scale_x_continuous(trans = x.trans,
                                limits=x.limits)
  }
  if (!is.null(y.limits)) {
    p <- p + scale_y_continuous(trans = y.trans,
                                limits=y.limits)
  }
  # p <- p + scale_x_continuous(trans = x.trans)
  # p <- p + scale_y_continuous(trans = y.trans)
  p <- p + coord_trans(x = x.coord.trans, y = y.coord.trans)
  
  if (x.trans == "log10" || y.trans == "log10") {
    sidesStr <- paste0(rep('l',as.numeric(y.trans == "log10")),
                       rep('b',as.numeric(x.trans == "log10")))
    p <- p + annotation_logticks(sides=sidesStr)
  }
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  return(p)
}

scatterPanelPlot <- function(simData, 
                             axis1var, 
                             axis2var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             facet1rev=FALSE,
                             facet2rev=FALSE,
                             colorvar=NULL,
                             shapevar=NULL,
                             x.trans="identity",
                             y.trans="identity",
                             x.coord.trans="identity",
                             y.coord.trans="identity",
                             x.limits=NULL,
                             y.limits=NULL,
                             seg.coords=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  # x.coord.trans = x.trans
  # x.coord.trans = y.trans
  saveFile <- genSaveFileName("Scatter",
                              axis1var, 
                              axis2var,
                              facet1vars=facet1vars,
                              facet2vars=facet2vars,
                              colorvar=colorvar,
                              shapevar=shapevar,
                              x.trans=x.trans,
                              y.trans=y.trans,
                              x.coord.trans=x.coord.trans,
                              y.coord.trans=y.coord.trans,
                              x.limits=x.limits,
                              y.limits=y.limits,
                              seg.coords=seg.coords,
                              outPrefix = outPrefix)
  
  axis1lab <- axis1var
  if(!is.na(axisLabels[axis1var])){
    axis1lab <- axisLabels[axis1var]
  }
  axis2lab <- axis2var
  if(!is.na(axisLabels[axis2var])){
    axis2lab <- axisLabels[axis2var]
  }
  
  p <- plotCanvas(simData,
                  axis1lab, 
                  axis2lab,
                  facet1vars=facet1vars,
                  facet2vars=facet2vars,
                  facet1rev=facet1rev,
                  facet2rev=facet2rev,
                  colorvar=colorvar,
                  shapevar=shapevar,
                  x.trans=x.trans,
                  y.trans=y.trans,
                  x.coord.trans=x.coord.trans,
                  y.coord.trans=y.coord.trans,
                  x.limits=x.limits,
                  y.limits=y.limits,
                  seg.coords=seg.coords)
  
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  
  # if (!is.null(shapevar)) {
  #   p <- p + geom_point(aes(x=!!a1,
  #                           y=!!a2),
  #                       # color=colorvar,
  #                       size=.25) #shape=21
  # } else {
  #   p <- p + geom_point(aes(x=!!a1,
  #                           y=!!a2),
  #                       # color=colorvar,
  #                       shape=21,
  #                       size=.25)
  # }
  
  if (!is.null(shapevar)) {
    s <- sym(shapevar)
    if (!is.null(colorvar)) {
      c <- sym(colorvar)
      p <- p + geom_point(aes(x=!!a1,
                              y=!!a2,
                              shape=!!s,
                              color=!!c),
                          size=.25) #shape=21
    } else {
      p <- p + geom_point(aes(x=!!a1,
                              y=!!a2,
                              shape=!!s),
                          size=.25) #shape=21
    }
  } else {
    
    if (!is.null(colorvar)) {
      c <- sym(colorvar)
      p <- p + geom_point(aes(x=!!a1,
                              y=!!a2,
                              color=!!c),
                          shape=21,
                          size=.25)
    } else {
      p <- p + geom_point(aes(x=!!a1,
                              y=!!a2),
                          shape=21,
                          size=.25)
    }
  }
  
  # p <- p + geom_point(aes(x=!!a1,
  #                         y=!!a2),
  #                     color=colorvar,
  #                     shape=shapevar,
  #                     size=.25) #shape=21
  
  p %>% savePlot(saveFile,width = width,height = height)
}

scatterPanelPlotB <- function(simData, 
                             axis1var, 
                             axis2var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             facet1rev=FALSE,
                             facet2rev=FALSE,
                             colorvar=NULL,
                             shapevar=NULL,
                             x.trans="identity",
                             y.trans="identity",
                             x.coord.trans="identity",
                             y.coord.trans="identity",
                             x.limits=NULL,
                             y.limits=NULL,
                             seg.coords=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  saveFile <- paste0(outPrefix,'Scatter.',axis1var,'.',axis2var)
  if (is.null(colorvar)){
    p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else {
    c <- sym(colorvar)
    p <- ggplot(simData, aes(x=!!a1, y=!!a2, color=!!c))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  
  p <- p +
    # geom_point(aes(x=log10(!!a1), y=log10(!!a2)),size=.25, shape=21) +
    geom_point(size=.25, shape=21) +
    xlab(axisLabels[axis1var]) +
    ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  if (!is.null(colorvar)){
    maxColorVarVal <- max(select(simData,colorvar))
    minColorVarVal <- min(select(simData,colorvar))
    # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
    
    p <- p +
      scale_colour_gradientn(
        name=axisLabels[colorvar],
        # name=NULL,
        limits = c(minColorVarVal,maxColorVarVal),
        colors = hcl.colors(5, palette = "Plasma"),
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      # scale_colour_gradient(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "#4B0055",
      #   high = "#FDE333",
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) +
      # scale_colour_gradient2(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "red",
      #   mid = "white",
      #   high = "blue",
      #   midpoint = midColorVarVal,
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) + 
      guides(
        # size = "none",
        colour = guide_colourbar(title.position = "right")
      )+
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  
  if (!is.null(seg.coords)) {
    # p <- p + geom_segment(0,0,1,1)
    # p <- p + geom_segment(seg.coords[1],seg.coords[2],seg.coords[3],seg.coords[4])
    # p <- p + geom_segment(seg.coords[1],seg.coords[2],seg.coords[3],seg.coords[4],
    #                       mapping=aes(linetype="dashed"))
    p <- p + geom_segment(aes(x = seg.coords[1],
                              y = seg.coords[2],
                              xend = seg.coords[3],
                              yend = seg.coords[4]),
                          alpha=.01,
                          color="gray25",
                          # color="#40404009",
                          # scale_color_manual("#404040"),
                          linetype="dashed")
  }
  
  if (!is.null(x.limits)) {
    p <- p + scale_x_continuous(trans = x.trans,
                                limits=x.limits)
  }
  if (!is.null(y.limits)) {
    p <- p + scale_y_continuous(trans = y.trans,
                                limits=y.limits)
  }
  # p <- p + scale_x_continuous(trans = x.trans)
  # p <- p + scale_y_continuous(trans = y.trans)
  p <- p + coord_trans(x = x.coord.trans, y = y.coord.trans)
  
  if (x.trans == "log10" || y.trans == "log10") {
    sidesStr <- paste0(rep('l',as.numeric(y.trans == "log10")),
                       rep('b',as.numeric(x.trans == "log10")))
    p <- p + annotation_logticks(sides=sidesStr)
  }
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




scatterPlotFinalGenParams <- function(fullFinalGenStatTable, 
                                      axis1var, 
                                      axis2var,
                                      colorvar=NULL,
                                      outPrefix = "./",
                                      width = 6,
                                      height = 5){
  for (m in unique(fullFinalGenStatTable$MutRate)){
    for (i in unique(fullFinalGenStatTable$InvRate)){
      for (e in unique(fullFinalGenStatTable$EncounterNum)){
        for (r in unique(fullFinalGenStatTable$CrossoverMapLength)){
          for (c in unique(fullFinalGenStatTable$ConversionRate)){
            for (N in unique(fullFinalGenStatTable$PopSize)){
              specificParamData <- filter(fullFinalGenStatTable,
                                          MutRate == m,
                                          InvRate == i,
                                          EncounterNum == e,
                                          CrossoverMapLength == r,
                                          ConversionRate == c,
                                          PopSize == N)
              specificParams <- select(specificParamData[1,],
                                       c(MutRate:PopSize))
              specificOutPrefix <- paste0(outPrefix,'m',m,'i',i,'e',e,'r',r,'c',c,'N',N,'.')
              print(specificOutPrefix)
              scatterPanelPlot(specificParamData, 
                                       axis1var, 
                                       axis2var,
                                       colorvar = colorvar,
                                       outPrefix = specificOutPrefix,
                                       width = width,
                                       height = height)
            }
          }
        }
      }
    }
  }
  return()
}


scatterPanelMultiPlot <- function(simData, 
                             axis1vars, 
                             axis2vars,
                             axis1lab=NULL,
                             axis2lab=NULL,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             facet1rev=FALSE,
                             facet2rev=FALSE,
                             colorvar=NULL,
                             colors=NULL,
                             shapevar=NULL,
                             shapes=NULL,
                             x.trans="identity",
                             y.trans="identity",
                             x.coord.trans="identity",
                             y.coord.trans="identity",
                             x.limits=NULL,
                             y.limits=NULL,
                             seg.coords=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  if(is.null(axis1lab)){
    axis1lab <- axis1vars[[1]]
    if(!is.na(axisLabels[axis1vars[[1]]])){
      axis1lab <- axisLabels[axis1vars[[1]]]
    }
  }
  if(is.null(axis2lab)){
    axis2lab <- axis2vars[[1]]
    if(!is.na(axisLabels[axis2vars[[1]]])){
      axis2lab <- axisLabels[axis2vars[[1]]]
    }
  }
  saveFile <- genSaveFileName("Scatter",
                              axis1lab, 
                              axis2lab,
                              facet1vars=facet1vars,
                              facet2vars=facet2vars,
                              colorvar=colorvar,
                              shapevar=shapevar,
                              x.trans=x.trans,
                              y.trans=y.trans,
                              x.coord.trans=x.coord.trans,
                              y.coord.trans=y.coord.trans,
                              x.limits=x.limits,
                              y.limits=y.limits,
                              seg.coords=seg.coords,
                              outPrefix = outPrefix)
  # p <- scatterPanelPlotter(simData,
  #                          axis1vars[[1]], 
  #                          axis2vars[[1]],
  #                          facet1vars=facet1vars,
  #                          facet2vars=facet2vars,
  #                          colorvar=colorvar,
  #                          x.trans=x.trans,
  #                          y.trans=y.trans,
  #                          x.coord.trans=x.coord.trans,
  #                          y.coord.trans=y.coord.trans,
  #                          x.limits=x.limits,
  #                          y.limits=y.limits,
  #                          seg.coords=seg.coords)
  p <- plotCanvas(simData,
                  axis1lab, 
                  axis2lab,
                  facet1vars=facet1vars,
                  facet2vars=facet2vars,
                  facet1rev=facet1rev,
                  facet2rev=facet2rev,
                  colorvar=colorvar,
                  shapevar=shapevar,
                  x.trans=x.trans,
                  y.trans=y.trans,
                  x.coord.trans=x.coord.trans,
                  y.coord.trans=y.coord.trans,
                  x.limits=x.limits,
                  y.limits=y.limits,
                  seg.coords=seg.coords)
  if (!is.null(shapes)) {
    for (i in seq(length(axis1vars))){
      p <- p +
        geom_point(aes(x=!!(sym(axis1vars[[i]])),
                       y=!!(sym(axis2vars[[i]]))),
                   color=colors[[i]],
                   shape=shapes[[i]],
                   size=.25) #shape=21
    }
  } else {
    for (i in seq(length(axis1vars))){
      p <- p +
        geom_point(aes(x=!!(sym(axis1vars[[i]])),
                       y=!!(sym(axis2vars[[i]]))),
                   color=colors[[i]],
                   shape=21,
                   size=.25)
    }
  }
  
  p %>% savePlot(saveFile,width = width,height = height)
}



scatterAnimation <- function(simData, 
                             axis1var, 
                             axis2var,
                             animvar,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             facet1rev=FALSE,
                             facet2rev=FALSE,
                             facet1.labs=NULL,
                             facet2.labs=NULL,
                             colorvar=NULL,
                             shapevar=NULL,
                             x.trans="identity",
                             y.trans="identity",
                             x.coord.trans="identity",
                             y.coord.trans="identity",
                             x.limits=NULL,
                             y.limits=NULL,
                             seg.coords=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  saveFile <- genSaveFileName("Animation.Scatter",
                              axis1var, 
                              axis2var,
                              facet1vars=facet1vars,
                              facet2vars=facet2vars,
                              colorvar=colorvar,
                              shapevar=shapevar,
                              x.trans=x.trans,
                              y.trans=y.trans,
                              x.coord.trans=x.coord.trans,
                              y.coord.trans=y.coord.trans,
                              x.limits=x.limits,
                              y.limits=y.limits,
                              seg.coords=seg.coords,
                              outPrefix = outPrefix)
  
  axis1lab <- axis1var
  if(!is.na(axisLabels[axis1var])){
    axis1lab <- axisLabels[axis1var]
  }
  axis2lab <- axis2var
  if(!is.na(axisLabels[axis2var])){
    axis2lab <- axisLabels[axis2var]
  }
  
  p <- plotCanvas(simData,
                  axis1lab, 
                  axis2lab,
                  facet1vars=facet1vars,
                  facet2vars=facet2vars,
                  facet1rev=facet1rev,
                  facet2rev=facet2rev,
                  facet1.labs=facet1.labs,
                  facet2.labs=facet2.labs,
                  colorvar=colorvar,
                  shapevar=shapevar,
                  x.trans=x.trans,
                  y.trans=y.trans,
                  x.coord.trans=x.coord.trans,
                  y.coord.trans=y.coord.trans,
                  x.limits=x.limits,
                  y.limits=y.limits,
                  seg.coords=seg.coords)
  
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  
  if (!is.null(colorvar)) {
    c <- sym(colorvar)
    p <- p + geom_point(aes(x=!!a1,
                            y=!!a2,
                            color=!!c),
                        size=.25) #shape=21
  } else {
    p <- p + geom_point(aes(x=!!a1,
                            y=!!a2),
                        shape=21,
                        size=.25)
  }
  
  t <- sym(animvar)
  
  # gganimate
  p <- p +
    labs(title = 'Generation: {floor(frame_time)}', x = axis1lab, y = axis2lab) +
    transition_time(!!t) +
    ease_aes('linear')
  
  # p %>% anim_save(paste0(saveFile,".gif"),width = width,height = height)
  # anim_save(paste0(saveFile,".gif"),animation = p)
  
  p <- animate(p,
               nframes=111,
               fps=20,
               # detail=3,
               # start_pause=10,
               end_pause=10,
               height = height,
               width = width,
               res = 150,
               units = "in")
  
  anim_save(paste0(saveFile,".gif"),animation = p)
}





violinPanelPlot <- function(simData, 
                            axis1var,
                            axis2vars,
                            facet1vars=NULL,
                            facet2vars=NULL,
                            color=FALSE,
                            outPrefix = "./",
                            width = 6,
                            height = 5) {
  a1 <- sym(axis1var)
  a21 <- sym(axis2vars[1])
  saveFile <- paste0(outPrefix,'Violin.',axis1var,'.',axis2vars[1])
  
  dualAxis <- (length(axis2vars) > 1)
  reqVars <- c(axis1var,facet1vars,facet2vars)
  tempSimData <- simData
  scale_shift_f <- function(x, scaleCoeff, shiftVal) {
    return(x * scaleCoeff - shiftVal)
  }
  inv_scale_shift_f <- function(x, scaleCoeff, shiftVal) {
    return( (x + shiftVal)/scaleCoeff )
  }
    
  if (dualAxis) {
    saveFile <- paste0(saveFile,'.',axis2vars[2])
    # a22 <- sym(axis2vars[2])
    maxA21 <- max(simData[axis2vars[1]])
    maxA22 <- max(simData[axis2vars[2]])
    minA21 <- min(simData[axis2vars[1]])
    minA22 <- min(simData[axis2vars[2]])
    scaleCoeff <- (maxA22-minA22)/(maxA21-minA21)
    shiftVal <- minA21-minA22
    # scaledVar <- paste0(axis2vars[2],"Scaled")
    scaledValVar <- "TempScaledVals"
    sourceVarCol <- "ValSource"
    sV <- sym(sourceVarCol)
    
    scaledSimData1 <- select(simData,all_of(c(reqVars,axis2vars[1])))
    scaledSimData1 <- cbind(scaledSimData1,axisLabels[axis2vars[1]])
    colnames(scaledSimData1) <- c(reqVars,scaledValVar,sourceVarCol)
    
    scaledSimData2 <- select(simData,all_of(c(reqVars,axis2vars[2])))
    scaledSimData2 <- cbind(scaledSimData2,axisLabels[axis2vars[2]])
    colnames(scaledSimData2) <- c(reqVars,scaledValVar,sourceVarCol)
    scaledSimData2[scaledValVar] <- inv_scale_shift_f(scaledSimData2[scaledValVar], scaleCoeff, shiftVal)
    
    tempSimData <- rbind(scaledSimData1,scaledSimData2)
    a22 <- sym(scaledValVar)
  }
  
  
  # saveFile <- paste0(outPrefix,'Histogram.',axis1var,'.',axis2var)
  
  # print(head(tempSimData,20))
  # print(tail(tempSimData,20))
  
  p <- ggplot(tempSimData, aes(x=!!a1))
  
  if (dualAxis) {
    p <- p +
      geom_violin(aes(y=!!a22, 
                      group=interaction(!!a1,!!sV), 
                      color=!!sV),
                  scale="count",
                  position="dodge") +
      # geom_boxplot(aes(y=!!a22, group=interaction(!!a1,!!sV), color=!!sV)) +
      geom_boxplot(aes(y=!!a22, 
                       group=interaction(!!a1,!!sV),
                       color=!!sV),
                   outlier.shape=1) +
                   # outlier.size=0.1,
                   # width=500,
                   # position=position_dodge(width=1000)) +
      # geom_boxplot(aes(y=!!a22, group=interaction(!!a1,!!sV), color=!!sV),width=~.*0.5,position=position_dodge(width=~.*1.1)) +
      scale_y_continuous(
        name = axisLabels[axis2vars[1]],
        sec.axis = sec_axis(~scale_shift_f(., scaleCoeff, shiftVal), name = axisLabels[axis2vars[2]])
      )
    # print(ggplot_build(p))
  } else {
    p <- p +
      geom_violin(aes(y=!!a21, group=!!a1)) +
      geom_boxplot(aes(y=!!a22, group=!!a1)) +
      ylab(axisLabels[axis2vars[1]])
  }
  
  p <- p +
    xlab(axisLabels[axis1var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  
  
  if (color) {
    saveFile <- paste0(saveFile,'.Colored')
    p <- p +
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




#
# Generating inversion SFS plots for neutral inversion dynamic comparisons
#


freqBinPanelPlot <- function(simData, 
                             axis1var, 
                             axis2var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             weightvar=NULL,
                             binwidths=NULL,
                             # colorvar=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  saveFile <- paste0(outPrefix,'FreqPlot.',axis1var,'.',axis2var)
  
  if (is.null(weightvar)){
    p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else {
    w <- sym(weightvar)
    p <- ggplot(simData, aes(x=!!a1, y=!!a2, weight=!!w))
  }
  
  if (is.null(binwidths)) {
    p <- p +
      geom_bin2d()
  } else {
    p <- p +
      geom_bin2d(binwidth=binwidths)
  }
  p <- p +
    xlab(axisLabels[axis1var]) +
    ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  p <- p +
    scale_fill_gradientn(
      name=axisLabels[paste0(axis1var,'-',axis2var)],
      # name=NULL,
      # limits = c(minColorVarVal,maxColorVarVal),
      colors = c("#FFFFFF",rev(hcl.colors(5, palette = "Blues")),rep(hcl.colors(1, palette = "Blues"),8)),
      # colors = c("#FFFFFF",rev(hcl.colors(4, palette = "Plasma")),rep(hcl.colors(1, palette = "Plasma"),3)),
      # colors = c("#FFFFFF",rev(hcl.colors(3, palette = "Plasma"))),
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    guides(
      # size = "none",
      colour = guide_colourbar(title.position = "right")
    )+
    theme(
      legend.direction = "vertical",
      legend.position = "right",
      # legend.key.height = unit(height/20, "in"),
      # legend.key.height = unit(2, "in"),
      legend.title = element_text(size = 10, angle = 90),
      legend.title.align = 0.5
    )
    # guides(
    #   colour = guide_colourbar(title.position = "right")
    # )+
    # theme(legend.title = element_text(size = 10, angle = 90),
    #   legend.title.align = 0.5,
    #   legend.key.height = unit(height/20, "in"),
    #   # legend.key.align = 0.5,
    #   legend.position = "right",
    #   legend.direction = "vertical"
    # )
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




histogramPanelPlot <- function(simData, 
                             axis1var,
                             facet1vars=NULL,
                             facet2vars=NULL,
                             weightvar=NULL,
                             colorvar=NULL,
                             binwidths=NULL,
                             x.trans="identity",
                             y.trans="identity",
                             x.coord.trans="identity",
                             y.coord.trans="identity",
                             x.limits=NULL,
                             y.limits=NULL,
                             outPrefix = "./",
                             width = 6,
                             height = 5) {
  a1 <- sym(axis1var)
  # a2 <- sym(axis2var)
  # saveFile <- paste0(outPrefix,'Histogram.',axis1var,'.',axis2var)
  saveFile <- paste0(outPrefix,'Histogram.',axis1var)
  
  if (is.null(weightvar) && is.null(colorvar)){
    p <- ggplot(simData, aes(x=!!a1))
    # p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else if (!is.null(weightvar) && is.null(colorvar)) {
    w <- sym(weightvar)
    p <- ggplot(simData, aes(x=!!a1, weight=!!w))
  } else if (is.null(weightvar) && !is.null(colorvar)) {
    c <- sym(colorvar)
    simData[[colorvar]] <- as.factor(simData[[colorvar]])
    # p <- ggplot(simData, aes(x=!!a1, fill=!!c, colour=!!c))
    p <- ggplot(simData, aes(x=!!a1, fill=!!c))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  } else {
    c <- sym(colorvar)
    # print(simData[[colorvar]])
    simData[[colorvar]] <- as.factor(simData[[colorvar]])
    # print(class(simData[[colorvar]]))
    # print(levels(simData[[colorvar]]))
    w <- sym(weightvar)
    # p <- ggplot(simData, aes(x=!!a1, fill=!!c, colour=!!c, weight=!!w))
    p <- ggplot(simData, aes(x=!!a1, fill=!!c, weight=!!w))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  
  if (is.null(binwidths)) {
    p <- p +
      geom_histogram()
  } else {
    p <- p +
      geom_histogram(binwidth = binwidths)
  }
  
  
  if (!is.null(x.limits)) {
    p <- p + scale_x_continuous(trans = x.trans)
  }
  if (!is.null(y.limits)) {
    p <- p + scale_y_continuous(trans = y.trans)
  }
  p <- p + coord_trans(xlim = x.limits, ylim = y.limits,
                       x = x.coord.trans, y = y.coord.trans)
  
  
  p <- p +
    xlab(axisLabels[axis1var]) +
    ylab("Count Across Simulations") +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  p <- p + theme(strip.text = element_text(face = "bold"))
  
  
  # if (!is.null(colorvar)){
  #   maxColorVarVal <- max(select(simData,colorvar))
  #   minColorVarVal <- min(select(simData,colorvar))
  #   # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
  #   p <- p +
  #     scale_colour_gradientn(
  #       name=axisLabels[paste0(colorvar)],
  #       # name=NULL,
  #       limits = c(minColorVarVal,maxColorVarVal),
  #       colors = c("#FFFFFF",rev(hcl.colors(5, palette = "Plasma"))),
  #       space = "Lab",
  #       na.value = "grey50",
  #       guide = "colourbar",
  #       aesthetics = "fill"
  #     ) +
  #     guides(
  #       colour = guide_colourbar(title.position = "right")
  #     )+
  #     theme(legend.title = element_text(size = 10, angle = 90),
  #           legend.title.align = 0.5,
  #           legend.key.height = unit(height/10, "in"),
  #           # legend.key.align = 0.5,
  #           legend.position = "right",
  #           legend.direction = "vertical"
  #     )
  # }
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}





# Designed for plotting the observed and HW expected homozygous individual frequencies of inv

xyFrequencyPlot <- function(simData, 
                            axis1var, 
                            axis2var,
                            facet1vars=NULL,
                            facet2vars=NULL,
                            colorvar=NULL,
                            x.trans="identity",
                            y.trans="identity",
                            x.coord.trans="identity",
                            y.coord.trans="identity",
                            x.limits=NULL,
                            y.limits=NULL,
                            outPrefix = "./",
                            width = 6,
                            height = 5) {
  a1 <- sym(axis1var)
  a2 <- sym(axis2var)
  saveFile <- paste0(outPrefix,'XY.',axis1var,'.',axis2var)
  if (is.null(colorvar)){
    p <- ggplot(simData, aes(x=!!a1, y=!!a2))
  } else {
    c <- sym(colorvar)
    p <- ggplot(simData, aes(x=!!a1, y=!!a2, color=!!c))
    saveFile <- paste0(saveFile,'.Colored.',colorvar)
  }
  
  p <- p +
    # geom_point(aes(x=log10(!!a1), y=log10(!!a2)),size=.25, shape=21) +
    geom_point(size=.25, shape=21) +
    xlab(axisLabels[axis1var]) +
    ylab(axisLabels[axis2var]) +
    theme(axis.text=element_text(size=14)) +
    theme_bw()
  
  if (!is.null(facet1vars) || !is.null(facet2vars)){
    if (!is.null(facet1vars)){
      facet1str <- paste(facet1vars,collapse = '+')
    } else {
      facet1str <- '.'
    }
    if (!is.null(facet2vars)){
      facet2str <- paste(facet2vars,collapse = '+')
    } else {
      facet2str <- '.'
    }
    facetForm <- paste(facet1str,facet2str,sep='~')
    p <- p + facet_grid(as.formula(facetForm))
    saveFile <- paste0(saveFile,'.Panelled.',facet1str,'.',facet2str)
  }
  
  if (!is.null(colorvar)){
    maxColorVarVal <- max(select(simData,colorvar))
    minColorVarVal <- min(select(simData,colorvar))
    # midColorVarVal <- (maxColorVarVal+minColorVarVal)/2
    
    p <- p +
      scale_colour_gradientn(
        name=axisLabels[colorvar],
        # name=NULL,
        limits = c(minColorVarVal,maxColorVarVal),
        colors = hcl.colors(5, palette = "Plasma"),
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
      ) +
      # scale_colour_gradient(
      #   name=axisLabels[colorvar],
      #   # name=NULL,
      #   limits = c(minColorVarVal,maxColorVarVal),
      #   low = "#4B0055",
      #   high = "#FDE333",
      #   space = "Lab",
      #   na.value = "grey50",
      #   guide = "colourbar",
      #   aesthetics = "colour"
      # ) +
    # scale_colour_gradient2(
    #   name=axisLabels[colorvar],
    #   # name=NULL,
    #   limits = c(minColorVarVal,maxColorVarVal),
    #   low = "red",
    #   mid = "white",
    #   high = "blue",
    #   midpoint = midColorVarVal,
    #   space = "Lab",
    #   na.value = "grey50",
    #   guide = "colourbar",
    #   aesthetics = "colour"
    # ) + 
    guides(
      # size = "none",
      colour = guide_colourbar(title.position = "right")
    )+
      theme(legend.title = element_text(size = 10, angle = 90),
            legend.title.align = 0.5,
            legend.key.height = unit(height/10, "in"),
            # legend.key.align = 0.5,
            legend.position = "right",
            legend.direction = "vertical"
      )
  }
  
  if (!is.null(x.limits)) {
    p <- p + scale_x_continuous(trans = x.trans,
                                limits=x.limits)
  }
  if (!is.null(y.limits)) {
    p <- p + scale_y_continuous(trans = y.trans,
                                limits=y.limits)
  }
  # p <- p + scale_x_continuous(trans = x.trans)
  # p <- p + scale_y_continuous(trans = y.trans)
  p <- p + coord_trans(x = x.coord.trans, y = y.coord.trans)
  
  if (x.trans == "log10" || y.trans == "log10") {
    sidesStr <- paste0(rep('l',as.numeric(y.trans == "log10")),
                       rep('b',as.numeric(x.trans == "log10")))
    p <- p + annotation_logticks(sides=sidesStr)
  }
  
  
  # if (max(select(simData,axis1var)) <= 1.0 && min(select(simData,axis1var)) >= 0.0){
  #   p <- p + scale_x_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  # if (max(select(simData,axis2var)) <= 1.0 && min(select(simData,axis2var)) >= 0.0){
  #   p <- p + scale_y_continuous(breaks=c(0,0.5,1),
  #                               limits=c(0,1))
  # }
  
  PNGsaveFile <- paste0(saveFile,".png")
  SVGsaveFile <- paste0(saveFile,".svg")
  ggsave(PNGsaveFile, plot = p, dpi=600, device = 'png',
         units="in", width=width, height=height)
  ggsave(SVGsaveFile, plot = p, dpi=600, device = 'svg',
         units="in", width=width, height=height)
}




