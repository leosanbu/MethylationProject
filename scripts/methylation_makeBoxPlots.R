args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape2))

## Script to make plots from extracted IPD values of predicted methylated motifs   ##
##                                                                                 ##
## July 2017                                                                       ##
##                                                                                 ##
## Leonor Sánchez-Busó (Infection Genomics)                                        ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire      ##

## Run this script using:
## Rscript methylation_makePlots.R /path/to/motif/folders query_motifs.txt
##
## IMPORTANT: order of arguments matters!
## /path/to/motif/folders: path leading to the folders containing *IPDs.txt files
## query_motifs.txt:       list of motifs to work with

pathtoipds <- args[1]
motifslist <- args[2]
ylimit <- args[3]
if (is.na(ylimit)){ ylimit <- 15 }  # for when you want to re-do one with another y-limit

## FUNCTIONS ##
readMotifs <- function(motifslist){
  queryFile <- scan(motifslist, what="\n")
  queryMat <- matrix(1:2, 1, 2)
  name <- ""
  seq <- ""
  for (i in 1:length(queryFile)){
    g <- grep("^>", queryFile[i])
    if (length(g)>0){
      if (seq != ""){
        queryMat <- rbind(queryMat, c(name, seq))
      }
      name = queryFile[i]
      seq = ""
    }else{
      if (name != ""){
        seq = queryFile[i]
      }
    }
  }
  queryMat <- rbind(queryMat, c(name, seq))
  queryMat <- queryMat[-1,, drop=F]
  queryMat[,1] <- gsub(">", "", queryMat[,1])
  return(queryMat)
}
createDataFrames <- function(files){
  tabList <- list()
  nams <- c()
  for (i in 1:length(files)){
    strain <- gsub("_IPDs.txt", "", files[i])
    tab <- read.table(paste0(pathtoipds, "/", motifName, "/", files[i]), header=T)

    # Melt for basic plots
    tmpmelt <- melt(tab[,c(2,9:ncol(tab))])
    tabList[[i]] <- tmpmelt
    nams <- c(nams, strain)
  }
  names(tabList) <- nams
  return(list(tabList))
}
makeBoxPlots <- function(outList, motif, titlesize=16, textsize=16, dotsize=0.75, ylimit=ylimit){
  
  queryMotif <- motif
  
  ## Process list of data frames ##
  
  # Create and save plots #
  plotList <- as.list(1:length(outList))
  nams <- c()
  for (i in 1:length(outList)){
    tmpdf <- outList[[i]]
    nam <- names(outList)[i]
    nams <- c(nams, nam)
    
    pl <- ggplot(tmpdf, aes(x = variable, y = value, fill = IPDstrand))+
      geom_boxplot(show.legend=FALSE, outlier.size=.75)+xlab("")+ylab("IPD ratio")+
      ggtitle(nam)+ylim(c(0,ylimit))+
      scale_x_discrete(labels=strsplit(queryMotif, "")[[1]])+
      theme_bw()+
      theme(text = element_text(size=textsize),
            plot.title = element_text(lineheight=.8, face="bold", size=titlesize))
    
    class(pl) <- "ggplot"
    plotList[[i]] <- pl
  }
  names(plotList) <- nams
  plotList <- plotList[order(names(plotList))]
  
  return(plotList)
}

## MAIN ##

# Read motifs
queryMat <- readMotifs(motifslist)

for (i in 1:nrow(queryMat)){
  motifName <- queryMat[i,1]
  motifSeq <- queryMat[i,2]
  
  outfilename <- paste0(pathtoipds, "/", motifName, "/", motifSeq, ".pdf")
  if (file.exists(outfilename)){
    cat(motifName, "\t", outfilename, "exists, skipping...\n")
  }else{
    
    # Read input files
    files <- dir(paste0(pathtoipds, "/", motifName), pattern="_IPDs.txt")
    
    # Create data frames
    suppressMessages(dataframes <- createDataFrames(files))
    
    # Plot 1
    # Boxplots of IPDs per base/strand
    tabList <- dataframes[[1]]
    
    plots <- makeBoxPlots(tabList, motifSeq, ylimit=ylimit)
    pdf(outfilename, width=18, height=15, useDingbats=F)
    grid.arrange(grobs=plots, ncol=5, nrow=5)
    dev.off()
    
    cat(motifName, "\t", outfilename, "created!\n")
  }
}


