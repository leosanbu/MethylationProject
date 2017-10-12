library(gplots)
library(reshape2)
library(RColorBrewer)
library(dendextend)

## Script to detect differentially methylated genes in individual motifs       ##
##                                                                             ##
## October 2017                                                                ##
##                                                                             ##
## Leonor Sánchez-Busó (Infection Genomics)                                    ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire  ##


# Must be run for each motif individually #
setwd("../motif24/")

mot_enrich <- 'motif_enrichment.txt'
meth_enrich <- 'methylation_enrichment.txt'

mot <- read.csv(mot_enrich, header=T)
met <- read.csv(meth_enrich, header=T)

mot.df <- melt(mot)
met.df <- melt(met)
df <- cbind(mot.df, met.df$value)
# df <- df[which(df$value>0),]
colnames(df)[3:4] <- c("motifs", "methylated")
df$prop <- (df$methylated/df$motifs)
df$prop[which(is.nan(df$prop))] <- NA

# get clusters methylated at 100% in at least one cluster
fullfreq <- unique(df$cluster[which(df$prop==1)])
full.df <- df[which(df$cluster %in% fullfreq),]

mot.fullfreq <- mot[which(mot$cluster %in% fullfreq),]
met.fullfreq <- met[which(met$cluster %in% fullfreq),]
mat.fullfreq <- (met.fullfreq[,2:ncol(met.fullfreq)]/mot.fullfreq[,2:ncol(mot.fullfreq)])*100
rownames(mat.fullfreq) <- met.fullfreq$cluster

nummat <- as.numeric(as.vector(mat.fullfreq[,1]))
for (i in 2:ncol(mat.fullfreq)){
  tmp <- as.numeric(as.vector(mat.fullfreq[,i]))
  wh <- which(is.nan(tmp))
  tmp[wh] <- NA
  nummat <- cbind(nummat, tmp)
}

rownames(nummat) <- rownames(mat.fullfreq)
colnames(nummat) <- colnames(mat.fullfreq)
  
hm.na <- apply(nummat, 1, function(x) length(which(is.na(x))))
nummat <- nummat[which(hm.na==0),] # get clusters with complete data
colnames(nummat) <- gsub("_with_plasmids", "", colnames(nummat))
colnames(nummat) <- gsub("Ngonorrhoeae_", "", colnames(nummat))
mybreaks <- c(0, 25, 50, 75, 100)
mycols <- brewer.pal(5, "YlOrRd")

pdf("motif24_heatmap.pdf")
h <- heatmap.2(nummat, na.rm=TRUE, tracecol=NA, labRow=FALSE, col=mycols, scale="none", na.color = "grey50")
dev.off()
write.table(nummat, "proportion_matrix_differential_methylation_allstrains.csv", quote = F)

## Select the strains with the target motif methylated ##

# sel_strains <- c("Ngono_FA1090", "WHO_F", "WHO_P", "WHO_Y", "NCTC13798", "NCTC12700") # CCACC
# sel_strains <- c("NCTC13801", "WHO_U", "WHO_G", "WHO_N") # GACN7TGC
# sel_strains <- c("NCTC10931", "WHO_M") # GACN6TGC
# sel_strains <- c("Ngono_FA1090", "WHO_W") # GCAN8TGC
# sel_strains <- c("Ngono_FA1090", "WHO_O", "NCTC12700") # motif31
# sel_strains <- c("Ngono_FA1090", "WHO_O", "NCTC12700", "NCTC10928", "NCTC13801", "NCTC13799", "NCTC10931")
# sel_strains <- colnames(nummat)[-c(10:12)]
# sel_strains <- c("NCTC12700", "NCTC10931")
# sel_strains <- c("NCTC13799", "NCTC13798", "WHO_N")

nummat_sel <- nummat[,which(colnames(nummat) %in% sel_strains)]
hm.na <- apply(nummat_sel, 1, function(x) length(which(is.na(x))))
nummat_sel <- nummat_sel[which(hm.na==0),]

pdf("motif24_heatmap_activeRMs.pdf")
heatmap.2(nummat_sel, na.rm=TRUE, tracecol=NA, labRow=FALSE, col=mycols, scale="none")
dev.off()

ap <- rowSums(nummat_sel)
nummat_sel2 <- nummat_sel[which(ap<1400),]

pdf("motif24_heatmap_activeRMs_diff_WHO_only.pdf")
h <- heatmap.2(nummat_sel2, na.rm=TRUE, tracecol=NA, col=mycols, scale="none")
dev.off()
write.csv(nummat_sel2, "proportion_matrix_differential_methylation_WHO_only.csv", quote = F)




