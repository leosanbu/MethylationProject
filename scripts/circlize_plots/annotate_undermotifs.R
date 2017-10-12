args <- commandArgs(TRUE)

## Script to annotate genomic regions with underrepresentation of certain motifs    ##
##                                                                                  ##
## July 2017                                                                        ##
##                                                                                  ##
## Leonor Sánchez-Busó (Infection Genomics)                                         ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire       ##


get_GFF <- function(x){
  gff.tmp <- system(paste("less", x, "| grep \\\t"), intern=TRUE)
  gff.tmp <- strsplit(gff.tmp, "\t")
  gff.tmp <- do.call(rbind, gff.tmp)
}

gffiles <- "../GFF_files/" # Location of GFF files

strain <- args[1]#"WHO_F"            # strain name
gsize <- as.numeric(args[2])#2292467 # genome size

bedfiles <- dir(".", pattern=strain)
bedfiles <- bedfiles[grep("bed", bedfiles)]
g <- grep("features", bedfiles)
features <- bedfiles[g]
bedfiles <- bedfiles[-g]
bedfiles <- bedfiles[grep("_chr_motif", bedfiles)]

motiflist <- list()
for (x in 1:length(bedfiles)){
  bed <- read.table(bedfiles[x], header=T)
  if (strain=="FA1090"){
    bed <- bed[which(bed$chr=="NC_002946"),] 
  }else{
    bed <- bed[which(bed$chr==strain),]
  }
  
  one.every <- ceiling(gsize/nrow(bed))
  cat(strain, bedfiles[x], one.every, "\n")
  start <- seq(1, gsize, by=one.every)
  end <- seq(one.every, gsize, by=one.every)
  end <- c(end, gsize)
  coords <- cbind(start, end)
  coords <- cbind(coords, rep(1, nrow(coords)))
  colnames(coords)[3] <- "expected"
  
  # calculate observed number of motifs
  observed <- c()
  for (i in 1:nrow(coords)){
    ini <- as.numeric(as.vector(coords[i,1]))
    end <- as.numeric(as.vector(coords[i,2]))
    hmany <- length(which(bed$start %in% ini:end))
    observed <- c(observed, hmany)
  }
  coords <- cbind(coords, observed)
  
  # overrepresentation of motifs
  # >3 observed motifs (mean(=1)+3*sd(=1))
  # over <- coords[which(coords[,4]>3),]
  
  # underrepresentation of motifs
  # >3 consecutive regions == 0
  count <- 0
  countlist <- 0
  regionlist <- c()
  for (i in 1:length(observed)){
    ob <- observed[i]
    if (ob==0){ 
      count <- count+1 
    }else{
      if (count>4){
        countlist <- countlist+1
        region <- sort((i-1):(i-count))
        tmp <- coords[region,]
        co <- range(tmp[,1:2])
        len <- co[2]-(co[1]-1)
        if (len>=8000){
          regionlist <- c(regionlist, paste0(paste(co, collapse="-"), ":", len))
        }
      }
      count <- 0
    }
  }
  motiflist[[x]] <- regionlist
  if (x==length(bedfiles) & is.null(regionlist)){
    motiflist[[x]] <- NA
  }
}
names(motiflist) <- bedfiles

# Final table
tab <- motiflist[[1]]
tab <- strsplit(tab, ":")
tab <- do.call(rbind, tab)
co <- strsplit(tab[,1], "-")
co <- do.call(rbind, co)
tab <- cbind(rep(names(motiflist)[1], nrow(tab)), co, tab[,2])
for (i in 2:length(motiflist)){
  if (!is.null(motiflist[[i]])){
    if (!is.na(motiflist[[i]][1])){
      tmp <- motiflist[[i]]
      tmp <- strsplit(tmp, ":")
      tmp <- do.call(rbind, tmp)
      co <- strsplit(tmp[,1], "-")
      co <- do.call(rbind, co)
      tmp <- cbind(rep(names(motiflist)[i], nrow(tmp)), co, tmp[,2])
      tab <- rbind(tab, tmp)
    }
  }
}
tab <- tab[order(as.numeric(tab[,2])),]

underbed <- cbind(rep(strain, nrow(tab)), tab[,2:3], rep(strain, nrow(tab)), tab[,1], rep("black", nrow(tab)))
colnames(underbed) <- c("strain", "start", "end", "features", "value", "clr")
write.table(underbed, paste0(strain, "_unmeth8k.bed"), sep="\t", quote = F, row.names=F)

# Annotate regions 

if (strain=="WHO_F"){
  gffname <- paste0(gffiles, strain, ".gff")
}else if (strain=="FA1090"){
  gffname <- paste0(gffiles, "Ngono_FA1090.gff")
}else{
  gffname <- paste0(gffiles, strain, "_with_plasmids.gff")
}

gff <- get_GFF(gffname)
g <- grep("pCryptic", gff[,9])
g <- c(g, grep("pConjugative", gff[,9]))
g <- c(g, grep("pBlaTEM", gff[,9]))
g <- sort(unique(g))
if (length(g)>0){
  gff <- gff[-g,]
}
genes_ini <- as.numeric(gff[,4])
genes_end <- as.numeric(gff[,5])
for (i in 1:nrow(underbed)){
  start <- as.numeric(underbed[i,2])
  end <- as.numeric(underbed[i,3])
  wh <- which(genes_ini>=start & genes_end<=end)
  check_min <- as.numeric(gff[min(wh)-1,5])
  if (check_min>start){ wh <- c(min(wh)-1, wh) }
  check_max <- as.numeric(gff[max(wh)+1,4])
  if (check_max<end){ wh <- c(wh, max(wh)+1) }
  tmp_gff <- gff[wh,]
  cat(paste0("##", paste(underbed[i,c(1:3,5)], collapse=" "), "\n"), file=paste0(strain, "_unmeth8k.gff"), append=TRUE)
  write.table(tmp_gff, file=paste0(strain, "_unmeth8k.gff"), sep="\t", quote = F, col.names=F, row.names=F, append=TRUE)
}

