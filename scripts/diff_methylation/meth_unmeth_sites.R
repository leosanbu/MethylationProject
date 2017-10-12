args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(data.table))

## Script to study methylated and unmethylated sites                           ##
##                                                                             ##
## July 2017                                                                   ##
##                                                                             ##
## Leonor Sánchez-Busó (Infection Genomics)                                    ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire  ##

## This script contains several parts, not intended to run all together ##

fuzzpath <- args[1] # path containing fuzznuc files
pathfile <- args[2] # paths_to_modifications_file.txt
sample_size <- 10000
  
## FUNCTIONS ##

parseFuzznuc <- function(fuzznuc.file){
  
  pathtofile <- paste0(fuzzpath, "/", fuzznuc.file)
  f <- file(pathtofile, "rb")
  a <- readChar(f, file.info(pathtofile)$size, useBytes=T)
  fuzz <- strsplit(a,"\n",fixed=T,useBytes=T)[[1]]
  close(f)
  
  gseq <- grep("# Sequence", fuzz)
  gtab <- matrix(1:2, 1, 2)
  for (i in 1:length(gseq)){
    if (i==length(gseq)){
      gtab <- rbind(gtab, c(gseq[i], length(fuzz)))
    }else{
      gtab <- rbind(gtab, c(gseq[i], gseq[i+1]-1))
    }
  }
  gtab <- gtab[-1,]
  if (is.null(dim(gtab))){ gtab <- t(as.matrix(gtab)) }
  
  fuzzlist <- as.list(1:nrow(gtab))
  for (i in 1:nrow(gtab)){
    interval <- gtab[i,1]:gtab[i,2]
    fragment <- fuzz[interval]
    rm.hash <- grep("^#", fragment)
    if (length(rm.hash)>1){ fragment <- fragment[-rm.hash] }
    rm <- which(fragment=="")
    rm <- sort(c(rm, grep("Start", fragment)))
    if (length(rm)>0){ fragment <- fragment[-rm] }
    
    writeLines(fragment, "tmp.txt")
    split.fragment <- fread("tmp.txt")
    system("rm tmp.txt")
    
    fuzzlist[[i]] <- split.fragment
  }
  names.gseq <- gsub("# Sequence: ", "", fuzz[gseq])
  names.gseq <- gsub(".*contig", "contig", names.gseq)
  names.gseq <- gsub(" .*", "", names.gseq)
  names.gseq <- gsub("\\..*", "", names.gseq)
  names(fuzzlist) <- names.gseq
  return(fuzzlist)
}
buildDataFrame <- function(tab){
  gs <- as.vector(as.matrix(tab[which(tab$base=="G"),9]))
  cs <- as.vector(as.matrix(tab[which(tab$base=="C"),9]))
  ts <- as.vector(as.matrix(tab[which(tab$base=="T"),9]))
  as <- as.vector(as.matrix(tab[which(tab$base=="A"),9]))
  if (!is.na(sample_size)){
    gs <- sample(gs, sample_size, replace = F)
    cs <- sample(cs, sample_size, replace = F)
    ts <- sample(ts, sample_size, replace = F)
    as <- sample(as, sample_size, replace = F)
  }
  
  df <- data.frame(strain=rep(strain, length(gs)+length(cs)+length(ts)+length(as)), 
                   base=c(rep("G", length(gs)), rep("C", length(cs)),
                          rep("T", length(ts)), rep("A", length(as))),
                   ipd=c(gs, cs, ts, as))
  return(df)
}

## MAIN ##

paths <- read.table(pathfile)
fuzznucfiles <- dir(fuzzpath, pattern=".fuzznuc$")

# Get IPD values for unmodified bases (outside predicted motifs) #

for (i in 1:length(fuzznucfiles)){
  cat(fuzznucfiles[i], "\n")
  strain <- gsub(fuzznucfiles[i], pattern="\\..*", replacement="")
  g <- which(paths[,1]==strain)
  file <- as.character(paths[g,2])
  
  # Read modifications files from PacBio #    
  modifications <- fread(file, header=T, sep=",")  # Use 'fread' to load very big tables #
  modifications[[1]] <- gsub(" .*", "", modifications[[1]])
  modifications[[1]] <- gsub(".*\\|", "", modifications[[1]])
  
  # Read fuzznuc file #
  fuzznuc <- parseFuzznuc(fuzznucfiles[i])
  contigs <- names(fuzznuc)
  
  # From modifications file, exclude all bases in the detected motifs #
  # Get all bases in the motifs #
  for (con in 1:length(contigs)){
    contig <- fuzznuc[[con]]
    mod <- modifications[which(modifications[[1]]==contigs[con]),]
    all_coordinates <- sort(unique(unlist(apply(contig, 1, function(y) y[1]:y[2]))))
    if (con == 1){
      modif_outside <- mod[which(!mod$tpl %in% all_coordinates),]
    }else{
      modif_outside <- rbind(modif_outside, mod[which(!mod$tpl %in% all_coordinates),])
    }
  }
  if (i == 1){
    outside_df <- buildDataFrame(modif_outside)
  }else{
    outside_df <- rbind(outside_df, buildDataFrame(modif_outside))
  }
}

write.table(outside_df, "IPDs_unmethylated_sites_1M_sampled.txt", sep="\t", quote = F, row.names=F)


# Get IPD values for modified bases (within the predicted motifs) #

selection <- c("motif1","motif6","motif10","motif11","motif13","motif14","motif16","motif18","motif19","motif21","motif23","motif24","motif31")

meth.df <- data.frame()
for (i in 1:length(selection)){
  motif <- selection[i]
  ipdfiles <- dir(paste0("../IPDs_per_motif/", motif), pattern=".txt")
  for (f in 1:length(ipdfiles)){
    file <- read.table(paste0("../IPDs_per_motif/", motif, "/", ipdfiles[f]), header=T)
    strain <- gsub(paste0("_", motif, "_IPDs.txt"), "", ipdfiles[f])
    cat(strain, "\t", motif, "\t", floor(nrow(file)/2), "\n")
    if (motif=="motif1"){
      methbase1 <- file$C[which(file$IPDstrand=="IPD0")]
      methbase2 <- file$G.1[which(file$IPDstrand=="IPD1")]
      len <- length(methbase1)+length(methbase2)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("4mC", len), ipd=c(methbase1, methbase2))
    }else if (motif=="motif6"){
      methbase1 <- file$A[which(file$IPDstrand=="IPD0")]
      len <- length(methbase1)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("6mA", len), ipd=c(methbase1))
    }else if (motif %in% c("motif10", "motif11", "motif13", "motif14", "motif16", "motif19")){
      methbase1 <- file$A[which(file$IPDstrand=="IPD0")]
      methbase2 <- file$T[which(file$IPDstrand=="IPD1")]
      len <- length(methbase1)+length(methbase2)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("6mA", len), ipd=c(methbase1, methbase2))
    }else if (motif=="motif18"){
      methbase1 <- file$A.1[which(file$IPDstrand=="IPD0")]
      len <- length(methbase1)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("6mA", len), ipd=c(methbase1))
    }else if (motif=="motif21"){
      methbase1 <- file$C[which(file$IPDstrand=="IPD0")]
      methbase2 <- file$G.2[which(file$IPDstrand=="IPD1")]
      len <- length(methbase1)+length(methbase2)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("5mC", len), ipd=c(methbase1, methbase2))
    }else if (motif=="motif23"){
      methbase1 <- file$C[which(file$IPDstrand=="IPD0")]
      len <- length(methbase1)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("4mC", len), ipd=c(methbase1))
    }else if (motif %in% c("motif24", "motif31")){
      methbase1 <- file$A[which(file$IPDstrand=="IPD0")]
      len <- length(methbase1)
      cur.df <- data.frame(strain=rep(strain, len), motif=rep(motif, len), base=rep("6mA", len), ipd=c(methbase1))
    }
    # Add to general df
    if (i==1 & f==1){
      meth.df <- cur.df
    }else{
      meth.df <- rbind(meth.df, cur.df)
    }
  }
}

write.table(meth.df, "IPDs_target_methylated_sites.txt", sep="\t", quote = F, row.names=F)

meth.filt <- meth.df[which(meth.df$ipd>2 & meth.df$ipd<15),]
meth.filt <- meth.filt[,-2]

outside.G <- outside_df[which(outside_df$base=="G"),]
outside.C <- outside_df[which(outside_df$base=="C"),]
outside.A <- outside_df[which(outside_df$base=="A"),]
outside.T <- outside_df[which(outside_df$base=="T"),]

meth.filt <- rbind(meth.filt, outside.G, outside.C, outside.A, outside.T)

ggplot(meth.filt, aes(x=base, y=ipd, fill=base))+
  geom_boxplot(outlier.size=.75, outlier.colour="grey70", alpha=.5)+
  ylim(c(0,10))+
  xlim("6mA", "4mC", "5mC", "A", "C", "G", "T")+
  xlab("")+ylab("IPD ratio")+
  geom_hline(yintercept = 2, linetype="dashed", col="grey70")
ggsave("boxplots.pdf", useDingbats=F, width=6, height=5)

methylated <- rep("unmethylated", nrow(meth.filt))
methylated[which(meth.filt$base=="4mC")] <- "methylated"
methylated[which(meth.filt$base=="5mC")] <- "methylated"
methylated[which(meth.filt$base=="6mA")] <- "methylated"
meth.filt$meth <- methylated
meth.filt$base <- gsub("6m", "", meth.filt$base)

#A
meth.filt.a <- meth.filt[which(meth.filt$base=="A"),]

dummy <- data.frame(X = c("methylated", "unmethylated"), 
                    Z = c(mean(meth.filt.a$ipd[which(meth.filt.a$meth=="methylated")]),
                          mean(meth.filt.a$ipd[which(meth.filt.a$meth=="unmethylated")])))
dummy$X <- as.factor(dummy$X)
# 1   methylated 5.6299347
# 2 unmethylated 0.9766355

# > quantile(meth.filt.a$ipd[which(meth.filt.a$meth=="methylated")])
# 0%    25%    50%    75%   100% 
# 2.001  4.318  5.438  6.699 14.993 

p1<-ggplot(meth.filt.a, aes(x=ipd, fill=meth))+
  geom_density(alpha=.5)+
  facet_grid(meth~base)+
  xlim(c(0,10))+
  geom_vline(data=dummy, aes(xintercept=Z), linetype="dashed")

#C
meth.filt.c <- meth.filt[which(meth.filt$base %in% c("4mC", "5mC", "C")),]
meth.filt.c$meth[which(meth.filt.c$base=="4mC")] <- "4mC"
meth.filt.c$meth[which(meth.filt.c$base=="5mC")] <- "5mC"
meth.filt.c$base <- "C"
meth.filt.c$base <- as.factor(meth.filt.c$base)
meth.filt.c$meth <- as.factor(meth.filt.c$meth)

dummy <- data.frame(X = c("4mC", "5mC", "unmethylated"),
                    Z = c(mean(meth.filt.c$ipd[which(meth.filt.c$meth=="4mC")]),
                          mean(meth.filt.c$ipd[which(meth.filt.c$meth=="5mC")]),
                          mean(meth.filt.c$ipd[which(meth.filt.c$meth=="unmethylated")])))
dummy$X <- as.factor(dummy$X)
# 1          4mC 3.244325
# 2          5mC 2.814361
# 3 unmethylated 1.197051

# > quantile(meth.filt.c$ipd[which(meth.filt.c$meth=="4mC")])
# 0%     25%     50%     75%    100% 
# 2.0010  2.3140  2.7780  3.5845 14.6460 
# 
# > quantile(meth.filt.c$ipd[which(meth.filt.c$meth=="5mC")])
# 0%    25%    50%    75%   100% 
# 2.001  2.178  2.444  2.978 14.998 

p2<-ggplot(meth.filt.c, aes(x=ipd, fill=meth))+
  geom_density(alpha=.5)+
  facet_grid(meth~base)+
  xlim(c(0,10))+
  geom_vline(data=dummy, aes(xintercept=Z), linetype="dashed")

#G-T
meth.filt.gt <- meth.filt[which(meth.filt$base %in% c("G", "T")),]

dummy <- data.frame(X = c("G", "T"),
                    Z = c(median(meth.filt.gt$ipd[which(meth.filt.gt$base=="G")]),
                          median(meth.filt.gt$ipd[which(meth.filt.gt$base=="T")])))
dummy$X <- as.factor(dummy$X)
# 1 G 1.374061
# 2 T 1.135675
p3<-ggplot(meth.filt.gt, aes(x=ipd, fill=meth))+
  geom_density(alpha=.5)+
  facet_grid(base~.)+
  xlim(c(0,10))+
  geom_vline(data=dummy, aes(xintercept=Z), linetype="dashed")
  
pdf("density_plots.pdf", width=5, height=8)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()

unmeth <- read.table("../../methyl_vs_unmethyl/IPDs_unmethylated_sites_1M_sampled.txt", header=T)
outside.C <- unmeth[which(unmeth$base=="C"),]


### Individual cases ###

## motif1

setwd("motif1/")
motif1files <- dir(".", pattern="_motif1_IPDs.txt")
strains <- gsub("_motif1_IPDs.txt", "", motif1files)

for (i in 1:length(motif1files)){
  file <- read.table(motif1files[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file)),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])), 
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.1[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m4C methylation in AAANCGGTTNNC")
ggsave("motif1_m4C_IPDs_vs_unmethylated.pdf", width=8, height=5)

# c0 <- tab$ipd[which(tab$base=="C0")]
# c1 <- tab$ipd[which(tab$base=="C1")]
# cu <- tab$ipd[which(tab$base=="C")]
# 
# wilcox.test(c0, cu)$p.value
# # 1.520415e-58
# 
# wilcox.test(c1, cu)$p.value
# # 2.815364e-225

strains <- unique(tab$strain)
c0test <- c()
c1test <- c()
for (i in 1:length(strains)){
  subtab <- tab[which(tab$strain==strains[i]),]
  c0 <- subtab$ipd[which(subtab$base=="C0")]
  c1 <- subtab$ipd[which(subtab$base=="C1")]
  cu <- subtab$ipd[which(subtab$base=="C")]
  test1 <- wilcox.test(c0, cu)$p.value
  c0test <- c(c0test, test1)
  test2 <- wilcox.test(c1, cu)$p.value
  c1test <- c(c1test, test2)
}

tests <- cbind(as.vector(strains), c0test, c1test, p.adjust(as.numeric(c0test), "bonferroni"),
               p.adjust(as.numeric(c1test), "bonferroni"))
colnames(tests) <- c("strain", "C0_pval", "C1_pval", "C0_bonf", "C1_bonf")
write.table(tests, "pvalues_mC_motif1.txt", sep="\t", quote = F, row.names=F)

## motif23

setwd("../motif23")

motiffiles <- dir(".", pattern=".txt")
strains <- gsub("_motif23_IPDs.txt", "", motiffiles)

for (i in 1:length(files)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file)),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m4C methylation in GGCSCDND")
ggsave("motif23_m4C_IPDs_vs_unmethylated.pdf", width=8, height=5)

## motif21

setwd("../motif21")

motiffiles <- dir(".", pattern=".txt")
strains <- gsub("_motif21_IPDs.txt", "", motiffiles)

for (i in 1:length(files)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file)),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.2[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in GCGGND")
ggsave("motif21_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

## motif18

setwd("../motif18")

motiffiles <- dir(".", pattern=".txt")
strains <- gsub("_motif18_IPDs.txt", "", motiffiles)

for (i in 1:length(files)){
  file <- read.table(motiffiles[i], header=T)
  meth <- file[which(file$A.1>2 & file$IPDstrand=="IPD0"),]
  cat(strains[i], "\t", nrow(file[which(file$IPDstrand=="IPD0"),]), 
                  "\t", nrow(meth[which(meth$IPDstrand=="IPD0"),]), "\n")
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("A0", nrow(file[which(file$IPDstrand=="IPD0"),]))),
                       ipd=c(file$A.1[which(file$IPDstrand=="IPD0")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.A)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m6A methylation in GCAGA")
ggsave("motif18_m6A_IPDs_vs_unmethylated.pdf", width=8, height=5)


motiffiles <- dir("ngoAXII_3/", pattern=".txt")
strains <- gsub("_ngoAXII_3_IPDs.txt", "", motiffiles)

for (i in 1:length(files)){
  file <- read.table(paste0("ngoAXII_3/", motiffiles[i]), header=T)
  meth <- file[which(file$A.2>2 & file$IPDstrand=="IPD0"),]
  cat(strains[i], "\t", nrow(file[which(file$IPDstrand=="IPD0"),]), 
      "\t", nrow(meth[which(meth$IPDstrand=="IPD0"),]), "\n")
}

motiffiles <- dir("ngoAXII_4/", pattern=".txt")
strains <- gsub("_ngoAXII_4_IPDs.txt", "", motiffiles)

for (i in 1:length(files)){
  file <- read.table(paste0("ngoAXII_4/", motiffiles[i]), header=T)
  meth <- file[which(file$A.2>2 & file$IPDstrand=="IPD0"),]
  cat(strains[i], "\t", nrow(file[which(file$IPDstrand=="IPD0"),]), 
      "\t", nrow(meth[which(meth$IPDstrand=="IPD0"),]), "\n")
}

## ngoAI

setwd("../ngoAI/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAI_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.1[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAI_pvalues.txt"), sep="\t", quote = F, row.names=F)


## ngoAII

setwd("../ngoAII/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAII_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.1[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "fdr")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "fdr")
write.table(m, paste0("ngoAII_pvalues.txt"), sep="\t", quote = F, row.names=F)


## ngoAIII

setwd("../ngoAIII/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAIII_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C.1[which(file$IPDstrand=="IPD0")],
                             file$G.1[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAIII_pvalues.txt"), sep="\t", quote = F, row.names=F)

## ngoAIV

setwd("../ngoAIV/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAIV_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.2[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in GCCGGC")
ggsave("ngoAIV_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

ggplot(tab, aes(ipd, fill=base))+
  geom_density(alpha=.25)+
  facet_grid(strain~.)+xlim(c(0,5))+
  theme(axis.text.y=element_text(size=5))
ggsave("test.pdf", width=3, height=10)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAIV_pvalues.txt"), sep="\t", quote = F, row.names=F)


## ngoAVII

setwd("../ngoAVII/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAVII_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.2[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in GCGGC")
ggsave("ngoAVII_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

ggplot(tab, aes(ipd, fill=base))+
  geom_density(alpha=.25)+
  facet_grid(strain~.)+xlim(c(0,5))+
  theme(axis.text.y=element_text(size=5))
ggsave("test.pdf", width=3, height=10)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAVII_pvalues.txt"), sep="\t", quote = F, row.names=F)


# NgoAXV

setwd("../ngoAXV/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAXV_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C[which(file$IPDstrand=="IPD0")],
                             file$G.1[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in GGNNCC")
ggsave("ngoAXV_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

ggplot(tab, aes(ipd, fill=base))+
  geom_density(alpha=.25)+
  facet_grid(strain~.)+xlim(c(0,5))+
  theme(axis.text.y=element_text(size=5))
ggsave("test.pdf", width=3, height=10)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAXV_pvalues.txt"), sep="\t", quote = F, row.names=F)

# NgoAXIV

setwd("../ngoAXIV/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_ngoAXIII_IPDs.txt", "", motiffiles) # they are really CCGG-ngoAXIV

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD0"),])),
                       base=c(rep("C0", nrow(file[which(file$IPDstrand=="IPD0"),])),
                              rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$C.1[which(file$IPDstrand=="IPD0")],
                             file$G[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in CCGG")
ggsave("ngoAXIV_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

m <- matrix(1:3, 1, 3)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab0 <- df[which(df$strain==strain & df$base=="C0"),]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest0 <- wilcox.test(subtab0$ipd, unmet)$p.value
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest0, ttest1))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C0", "pvalue_C1")

m <- as.data.frame(m)

m$bonf_pvalue_C0 <- p.adjust(as.numeric(as.vector(m$pvalue_C0)), "bonferroni")
m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
write.table(m, paste0("ngoAXIV_pvalues.txt"), sep="\t", quote = F, row.names=F)

## NgoAXVI

setwd("../motif24/")
motiffiles <- dir(".", pattern="IPDs.txt")
strains <- gsub("_motif24_IPDs.txt", "", motiffiles)

for (i in 1:length(motiffiles)){
  file <- read.table(motiffiles[i], header=T)
  tmp.df <- data.frame(strain=rep(strains[i], nrow(file[which(file$IPDstrand=="IPD1"),])),
                       base=c(rep("C1", nrow(file[which(file$IPDstrand=="IPD1"),])),
                              rep("C2", nrow(file[which(file$IPDstrand=="IPD1"),])),
                              rep("C3", nrow(file[which(file$IPDstrand=="IPD1"),]))),
                       ipd=c(file$G[which(file$IPDstrand=="IPD1")],
                             file$G.1[which(file$IPDstrand=="IPD1")],
                             file$G.2[which(file$IPDstrand=="IPD1")]))
  if (i==1){
    df <- tmp.df
  }else{
    df <- rbind(df, tmp.df)
  }
}

tab <- rbind(df, outside.C)

ggplot(tab, aes(x=strain, y=ipd, fill=base))+
  geom_boxplot(alpha=.5, outlier.shape = NA)+ylim(c(0,3))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  xlab("Strain")+ylab("IPD ratio")+
  ggtitle("m5C methylation in GGTGA")
ggsave("ngoAXVI_m5C_IPDs_vs_unmethylated.pdf", width=8, height=5)

m <- matrix(1:4, 1, 4)
for (i in 1:length(strains)){
  strain <- strains[i]
  subtab1 <- df[which(df$strain==strain & df$base=="C1"),]
  subtab2 <- df[which(df$strain==strain & df$base=="C2"),]
  subtab3 <- df[which(df$strain==strain & df$base=="C3"),]
  unmet <- outside.C$ipd[which(outside.C$strain==strain)]
  ttest1 <- wilcox.test(subtab1$ipd, unmet)$p.value
  ttest2 <- wilcox.test(subtab2$ipd, unmet)$p.value
  ttest3 <- wilcox.test(subtab3$ipd, unmet)$p.value
  m <- rbind(m, c(strain, ttest1, ttest2, ttest3))
}
m <- m[-1,]
colnames(m) <- c("strain", "pvalue_C1", "pvalue_C2", "pvalue_C3")

m <- as.data.frame(m)

m$bonf_pvalue_C1 <- p.adjust(as.numeric(as.vector(m$pvalue_C1)), "bonferroni")
m$bonf_pvalue_C2 <- p.adjust(as.numeric(as.vector(m$pvalue_C2)), "bonferroni")
m$bonf_pvalue_C3 <- p.adjust(as.numeric(as.vector(m$pvalue_C3)), "bonferroni")

write.table(m, paste0("ngoAXVI_pvalues.txt"), sep="\t", quote = F, row.names=F)

m2 <- m[,c(1,5:7)]
m2 <- melt(m2)

ggplot(m2, aes(x=variable, y=-log10(value), col=variable))+
  geom_boxplot()
ggsave("motif24_mC_summary.pdf", useDingbats=F, width=6, height=5)



