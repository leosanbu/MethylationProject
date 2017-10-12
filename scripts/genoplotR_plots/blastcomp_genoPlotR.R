args <- commandArgs(TRUE)
Sys.setlocale('LC_ALL','C') 
library(genoPlotR)

## Script to read crunch files from ACT and plot genome comparisons using genoplotR ##
##                                                                                  ##
## July 2017                                                                        ##
##                                                                                  ##
## Leonor Sánchez-Busó (Infection Genomics)                                         ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire       ##

files <- dir(".", pattern="crunch")
embl.path <- args[1]                          # Path to annotation files in EMBL format
extra.info <- read.table(args[2], header=T)   # Table with extra information to plot, i.e. GGI_extrainfo.tab
  
# Are the crunch files gzipped? #
wh <- grep(".crunch.gz", files)
if (length(wh)>0){
  for (i in wh){
    system(paste("gunzip", files[i]))
  }
}

# Read tree #
tree_file <- "RAxML_bestTree.ml_WHO_and_FA1090.midpoint.tre"
tree_str <- readLines(tree_file)
tree <- newick2phylog(tree_str)

# Load EMBL files #
# Keep tree order! #
select.embl <- paste(names(tree$leaves), ".embl", sep="")
embl <- dir(embl.path, pattern=".embl")

# Create and append list for all the embl files #
dna_segs <- list()
tags <- c("CDS", "rRNA", "tRNA", "tmRNA") 
for (i in 1:length(select.embl)){
  sel <- select.embl[i]
  cat("# Loading", sel, "...\n")
  cur <- read_dna_seg_from_embl(paste(embl.path, sel, sep=""), tagsToParse = tags)
  # Edit colours, set all to blue except for phages #
  cur$col <- "navy"
  phages <- grep("phage", cur$product)
  cur$col[phages] <- "pink"
  
  # Check if there is extra info and add it in that case #
  if(!is.null(extra.info)){
    sel <- gsub(".embl", "", sel)
    wh <- which(extra.info$strain==sel)
    if(length(wh)>0){
      # In case there are more than 1 extra info #
      for (x in wh){
        sub.info <- extra.info[x,]
        # Check which features overlap with the extra info and add 
        # the name to the product so then can be found with grep
        ini <- sub.info[1,2]
        end <- sub.info[1,3]
        wh.above <- which(cur$start>=ini)
        wh.below <- which(cur$end<=end)
        int <- intersect(wh.above, wh.below)
        cur[int,9] <- paste(cur[int,9], ",__", sub.info$name, "__", sep="")
        cur[int,13] <- as.character(sub.info$col)
        
        tmp <- rep("NA", length(cur[1,]))
        names(tmp) <- colnames(cur)
        rm <- which(!colnames(extra.info) %in% colnames(cur))
        if (length(rm)>0){ sub.info <- sub.info[,-rm] }
        tmp[names(sub.info)] <- as.matrix(sub.info)
        names(tmp) <- colnames(cur)
        tmp["strand"] <- 1
        tmp["feature"] <- "CDS"
        tmp["gene_type"] <- "bars"
        tmp["lty"] <- 1
        tmp["lwd"] <- 1
        tmp["pch"] <- 8
        tmp["cex"] <- 1
        tmp["length"] <- sub.info$end-(sub.info$start-1)
        cur <- rbind(cur, tmp)
      }
    }
  }
  cur$start <- as.numeric(cur$start)
  cur$end <- as.numeric(cur$end)
  cur$length <- as.numeric(cur$length)
  cur$lwd <- as.numeric(cur$lwd)
  cur$lty <- as.numeric(cur$lty)
  cur$pch <- as.numeric(cur$pch)
  cur$cex <- as.numeric(cur$cex)
  dna_segs[[i]] <- cur 
}
names(dna_segs) <- gsub(select.embl, pattern=".embl", replacement="")

# Adapt crunch output files and create comparison list #
comp.list <- list()
comp.tab <- cbind(names(dna_segs)[-length(names(dna_segs))], names(dna_segs)[-1])
files <- dir(".", pattern="crunch")
for (i in 1:nrow(comp.tab)){
  f1 <- comp.tab[i,1]
  f2 <- comp.tab[i,2]
  # Get the right comparison file #
  wh <- grep(f1, files)
  wh <- intersect(wh, grep(f2, files))
  f <- files[wh]
  cat(f1, f2, f, "\n")
  cur.tab <- read.table(f)
  # Remove hits <250 bp length <99% ID #
  aln.len <- cur.tab$V4-(cur.tab$V3-1)
  wh <- which(aln.len<1000)
  if (length(wh)>0){ cur.tab <- cur.tab[-wh,] }
  wh <- which(cur.tab$V2<99)
  if (length(wh)>0){ cur.tab <- cur.tab[-wh,] }
  
  # Set the right format #
  # Check query/subject order #
  g1 <- which(cur.tab[1,]==f1)
  g2 <- which(cur.tab[1,]==f2)
  if (g1<g2){
    cur.tab <- as.data.frame(cbind(cur.tab$V3, cur.tab$V4, cur.tab$V6, cur.tab$V7))
  }else{
    cur.tab <- as.data.frame(cbind(cur.tab$V6, cur.tab$V7, cur.tab$V3, cur.tab$V4))
  }
  
  colnames(cur.tab) <- c("start1", "end1", "start2", "end2")

  comp <- comparison(cur.tab)
  comp$col <- NA
  comp$col[which(comp$direction==1)] <- "red"
  comp$col[which(comp$direction==-1)] <- "royalblue"
  comp <- comp[,c(1:4, 6, 5)]
  comp.list[[i]] <- comp
}

pdf("WHO_refs_genoplotR.pdf", width=15, height=10)
plot_gene_map(dna_segs=dna_segs, comparisons=comp.list, tree=tree, tree_branch_labels_cex=0, 
              tree_scale=T, tree_width=3)

names(dna_segs) <- NULL
plot_gene_map(dna_segs=dna_segs, comparisons=comp.list)
dev.off()
