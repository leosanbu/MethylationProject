args <- commandArgs(TRUE)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(data.table)
library(RColorBrewer)
library(ape)

## Script to study genome-wide IPD ratios in annotated 5kb windows              ##
##                                                                              ##
## August 2017                                                                  ##
##                                                                              ##
## Leonor Sánchez-Busó (Infection Genomics)                                     ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire   ##

get_GFF <- function(x){
  gff.tmp <- system(paste("less", x, "| grep \\\t"), intern=TRUE)
  gff.tmp <- strsplit(gff.tmp, "\t")
  gff.tmp <- do.call(rbind, gff.tmp)
}

# Give as input file with paths to modifications.csv #
paths <- read.table(args[1]) # path_to_modifications_file.txt
my_palette <- colorRampPalette(c("white", "red"))(n = 100)
fna.paths <- args[2] # Path to .fna files
fnas <- dir(fna.paths, pattern=".fna")
gff.paths <- args[3] # path to .gff files

ord <- c("WHO_Z", "WHO_X", "WHO_W", "WHO_K", "WHO_Y", "WHO_V", "WHO_P", "WHO_M",
         "WHO_U", "WHO_N", "WHO_G", "WHO_O", "WHO_L", "WHO_F", "FA1090")

pdf("gwide_methylation.pdf")
m <- matrix(1:15, 15, 1)
layout(m)

# plots <- list()
for (i in 1:length(ord)){
  cat(ord[i], "\n")
  g <- which(paths[,1]==ord[i])
  modif <- fread(as.vector(paths[g,2]), header=T, sep=",")  # Use 'fread' to load very big tables #
  strain <- as.vector(paths[g,1])
  IPD0 <- modif[which(modif$strand==0),]
  IPD0 <- IPD0$ipdRatio
  IPD1 <- modif[which(modif$strand==1),]
  IPD1 <- IPD1$ipdRatio

  fas <- read.dna(paste0(fna.paths, strain, ".fna"), "fasta")
  
  # Genome-wide IPD ratios by window #
  coords <- coords_table(window=5000, step=2500, max=ncol(fas))
  if (coords[nrow(coords),2]>ncol(fas)){ coords[nrow(coords),2] <- ncol(fas) }
  
  # Check metadata here: phages, GGI, GIs #
  if (ord[i]=="WHO_F"){
    gffname <- paste0(gff.paths, "/", ord[i], ".gff")
  }else if (ord[i]=="FA1090"){
    gffname <- paste0(gff.paths, "/Ngono_FA1090.gff")
  }else{
    gffname <- paste0(gff.paths, "/", ord[i], "_with_plasmids.gff")
  }
  gff <- get_GFF(gffname)
  if (ord[i]=="FA1090"){
    gff <- gff[which(gff[,1]=="NC_002946.2"),]
  }else{
    gff <- gff[which(gff[,1]==ord[i]),]
  }
  extrainfo <- read.table(paste0(gff.paths, "/", ord[i], "_chr_features.bed"), header=T)
  extrainfo <- extrainfo[which(extrainfo$value %in% c("GGI", "GI")),, drop=F]
  genes_ini <- as.numeric(gff[,4])
  genes_end <- as.numeric(gff[,5])
  col <- c()
  for (x in 1:nrow(coords)){
    ini <- as.numeric(coords[x,1])
    end <- as.numeric(coords[x,2])
    # check first in extrainfo
    if (nrow(extrainfo)>0){
      outco <- ""
      for (y in 1:nrow(extrainfo)){
        st <- as.numeric(as.vector(extrainfo[y,2]))
        fi <- as.numeric(as.vector(extrainfo[y,3]))
        wh <- length(which(ini:end %in% st:fi))
        if (wh>0){
          if (extrainfo[y,5]=="GGI"){
            co <- "goldenrod1"
          }else{
            co <- "darkorange2"
          }
          
        }else{
          # check GFF
          wh <- which(genes_ini>=ini & genes_end<=end)
          if (length(wh)>0){
            check_min <- as.numeric(gff[min(wh)-1,5])
            if (length(check_min)>0){
              if (check_min>ini){ wh <- c(min(wh)-1, wh) }
            }          
            check_max <- as.numeric(gff[max(wh)+1,4])
            if (length(check_max)>0){
              if (check_max<end){ wh <- c(wh, max(wh)+1) }
            }
            tmp_gff <- gff[wh,, drop=F]
            agv <- grep("opA54|opacity|piiC|pilE|pilS|mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
            if (length(agv)>0){
              co <- "yellowgreen"
            }else{
              phages <- grep("phage", tmp_gff[,9])
              if (length(phages)>0){
                co <- "darkmagenta"
              }else{
                rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
                if (length(rib)>0){
                  co <- "grey50"
                }else{
                  co <- ""
                }
              }        
            }
          }else{
            wh <- c(max(which(genes_end<=end)), min(which(genes_ini>=ini)))
            tmp_gff <- gff[wh,, drop=F]
            agv <- grep("opA54|opacity|piiC|pilE|pilS|mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
            if (length(agv)>0){
              co <- "yellowgreen"
            }else{
              phages <- grep("phage", tmp_gff[,9])
              if (length(phages)>0){
                co <- "darkmagenta"
              }else{
                rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
                if (length(rib)>0){
                  co <- "grey50"
                }else{
                  co <- ""
                }
              }        
            }
          }
        }
        if (co!=""){
          outco <- co
        }
      }
      col <- c(col, outco)
    }else{ # no extra info, check GFF
      # check GFF
      wh <- which(genes_ini>=ini & genes_end<=end)
      if (length(wh)>0){
        check_min <- as.numeric(gff[min(wh)-1,5])
        if (length(check_min)>0){
          if (check_min>ini){ wh <- c(min(wh)-1, wh) }
        }          
        check_max <- as.numeric(gff[max(wh)+1,4])
        if (length(check_max)>0){
          if (check_max<end){ wh <- c(wh, max(wh)+1) }
        }
        tmp_gff <- gff[wh,, drop=F]
        agv <- grep("opA54|opacity|piiC|pilE|pilS|mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
        if (length(agv)>0){
          col <- c(col, "yellowgreen")
        }else{
          phages <- grep("phage", tmp_gff[,9])
          if (length(phages)>0){
            col <- c(col, "darkmagenta")
          }else{
            rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
            if (length(rib)>0){
              col <- c(col, "grey50")
            }else{
              col <- c(col, "")
            }
          }        
        }
      }else{
        wh <- c(max(which(genes_end<=end)), min(which(genes_ini>=ini)))
        tmp_gff <- gff[wh,, drop=F]
        agv <- grep("opA54|opacity|piiC|pilE|pilS|mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
        if (length(agv)>0){
          col <- c(col, "yellowgreen")
        }else{
          phages <- grep("phage", tmp_gff[,9])
          if (length(phages)>0){
            col <- c(col, "darkmagenta")
          }else{
            rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
            if (length(rib)>0){
              col <- c(col, "grey50")
            }else{
              col <- c(col, "")
            }          
          }        
        }
      }
    }
  }
  
  getMean <- apply(coords, 1, function(x){
    zero <- mean(IPD0[x[1]:x[2]])
    one <- mean(IPD1[x[1]:x[2]])
    list(zero,one)
  })
  
  IPD0win <- unlist(lapply(getMean, function(x) x[[1]]))
  IPD1win <- unlist(lapply(getMean, function(x) x[[2]]))
  
  # Plot background #
  coords <- cbind(coords, col)
  coords.filt <- coords[which(col!=""),]
  
  par(mar=c(0,0,0,0))
  plot(x=as.numeric(coords[,1]), y=rep(0, nrow(coords)), type="n", axes=F, xlab="", ylab="", ylim=c(0.9,1.25))
  segments(x0=as.numeric(coords.filt[,1]), x1=as.numeric(coords.filt[,1]),
           y0=rep(0.98, nrow(coords.filt)), y1=rep(1.25, nrow(coords.filt)), col=adjustcolor(coords.filt[,3], alpha=.5))
  
  # Genome-wide IPD patterns #
  par(new=T)
  plot(IPD0win, type="l", ylim=c(0.9,1.25), col="gray30", axes=F, xlab="", ylab="")
  par(new=T)
  plot(IPD1win, type="l", ylim=c(0.9,1.25), col="orangered", axes=F, xlab="", ylab="")

}
dev.off()

