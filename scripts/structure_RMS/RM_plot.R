setwd("RM_plots/")
library(genoPlotR)

tab <- read.table("rm_plots_data.txt", header=T, sep="\t")

rms <- unique(tab$rm)

list.segs <- list()
for (i in 1:length(rms)){
  subtab <- tab[which(tab$rm==rms[i]),-c(1,2)]
  tmp_segs <- dna_seg(subtab)
  list.segs[[i]] <- tmp_segs
}

pdf("RM_structure.pdf")
plot_gene_map(dna_segs=list.segs)
dev.off()
