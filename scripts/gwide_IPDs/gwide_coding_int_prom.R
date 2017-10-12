suppressMessages(library(data.table))

## Script to study IPD ratios in coding/non-coding/promoger regions             ##
##                                                                              ##
## September 2017                                                               ##
##                                                                              ##
## Leonor Sánchez-Busó (Infection Genomics)                                     ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire   ##

gffpath <- "GFF_files/"#args[1]
pathfile <- "paths_to_modifications_file.txt"#args[2]
unmeth <- read.table("methyl_vs_unmethyl/IPDs_unmethylated_sites_1M_sampled.txt", header=T)

## FUNCTIONS ##

get_GFF <- function(x){
  gff.tmp <- system(paste("less", x, "| grep \\\t"), intern=TRUE)
  gff.tmp <- strsplit(gff.tmp, "\t")
  gff.tmp <- do.call(rbind, gff.tmp)
}

## MAIN ##

# Get GFF files #
gfffiles <- dir(gffpath, pattern=".gff")

# Get paths to modifications files #
paths <- read.table(pathfile)

plotlist <- list()
for (i in 1:length(gfffiles)){
  strain <- gsub(".gff", "", gfffiles[i])
  strain <- gsub("Ngono_", "", strain)
  strain <- gsub("Ngonorrhoeae_", "", strain)
  strain <- gsub("_with_plasmids", "", strain)
  cat(strain, "\n")
  
  gff <- get_GFF(paste0(gffpath, "/", gfffiles[i]))
  
  modifications <- fread(as.vector(paths[which(paths[,1]==strain),2]), header=T, sep=",")
  coding_coords <- unique(unlist(apply(gff, 1, function(x) as.numeric(x[4]):as.numeric(x[5]))))
  promoter_coords <- unique(unlist(apply(gff, 1, function(x){
    strand <- x[7]
    if (strand=="+"){
      prom <- (as.numeric(x[4])-100):(as.numeric(x[4])-1)
    }else if (strand=="-"){
      prom <- (as.numeric(x[5])+1):(as.numeric(x[5])+100)
    }
  })))
  
  intergenic_coords <- 1:max(modifications$tpl)
  intergenic_coords <- intergenic_coords[which(!intergenic_coords %in% c(coding_coords, promoter_coords))]
  
  coding_modif0A <- modifications$ipdRatio[which(modifications$tpl %in% coding_coords & modifications$strand==0 & modifications$base=="A")]
  coding_modif0A <- coding_modif0A[which(coding_modif0A>3)]
  coding_modif0C <- modifications$ipdRatio[which(modifications$tpl %in% coding_coords & modifications$strand==0 & modifications$base=="C")]
  coding_modif0C <- coding_modif0C[which(coding_modif0C>3)]
  coding_modif1A <- modifications$ipdRatio[which(modifications$tpl %in% coding_coords & modifications$strand==1 & modifications$base=="A")]
  coding_modif1A <- coding_modif1A[which(coding_modif1A>3)]
  coding_modif1C <- modifications$ipdRatio[which(modifications$tpl %in% coding_coords & modifications$strand==0 & modifications$base=="C")]
  coding_modif1C <- coding_modif1C[which(coding_modif1C>3)]
  
  prom_modif0A <- modifications$ipdRatio[which(modifications$tpl %in% promoter_coords & modifications$strand==0 & modifications$base=="A")]
  prom_modif0A <- prom_modif0A[which(prom_modif0A>3)]
  prom_modif0C <- modifications$ipdRatio[which(modifications$tpl %in% promoter_coords & modifications$strand==0 & modifications$base=="C")]
  prom_modif0C <- prom_modif0C[which(prom_modif0C>3)]
  prom_modif1A <- modifications$ipdRatio[which(modifications$tpl %in% promoter_coords & modifications$strand==1 & modifications$base=="A")]
  prom_modif1A <- prom_modif1A[which(prom_modif1A>3)]
  prom_modif1C <- modifications$ipdRatio[which(modifications$tpl %in% promoter_coords & modifications$strand==0 & modifications$base=="C")]
  prom_modif1C <- prom_modif1C[which(prom_modif1C>3)]
  
  inter_modif0A <- modifications$ipdRatio[which(modifications$tpl %in% intergenic_coords & modifications$strand==0 & modifications$base=="A")]
  inter_modif0A <- inter_modif0A[which(inter_modif0A>3)]
  inter_modif0C <- modifications$ipdRatio[which(modifications$tpl %in% intergenic_coords & modifications$strand==0 & modifications$base=="C")]
  inter_modif0C <- inter_modif0C[which(inter_modif0C>3)]
  inter_modif1A <- modifications$ipdRatio[which(modifications$tpl %in% intergenic_coords & modifications$strand==1 & modifications$base=="A")]
  inter_modif1A <- inter_modif1A[which(inter_modif1A>3)]
  inter_modif1C <- modifications$ipdRatio[which(modifications$tpl %in% intergenic_coords & modifications$strand==0 & modifications$base=="C")]
  inter_modif1C <- inter_modif1C[which(inter_modif1C>3)]
  
  unmethA <- unmeth$ipd[which(unmeth$strain==strain & unmeth$base=="A")]
  unmethC <- unmeth$ipd[which(unmeth$strain==strain & unmeth$base=="C")]
  
  dfA <- data.frame(base="A",
                    type=c(rep("coding", length(coding_modif0A)+length(coding_modif1A)),
                           rep("promoter", length(prom_modif0A)+length(prom_modif1A)),
                           rep("intergenic", length(inter_modif0A)+length(inter_modif1A)),
                           rep("Unmeth", length(unmethA))),
                    strand=c(rep("IPD0", length(coding_modif0A)),
                             rep("IPD1", length(coding_modif1A)),
                             rep("IPD0", length(prom_modif0A)),
                             rep("IPD1", length(prom_modif1A)),
                             rep("IPD0", length(inter_modif0A)),
                             rep("IPD1", length(inter_modif1A)),
                             rep("Unmeth", length(unmethA))),
                   ipd=c(coding_modif0A, coding_modif1A, prom_modif0A, prom_modif1A, inter_modif0A, inter_modif1A, unmethA))
  
  dfC <- data.frame(base="C",
                    type=c(rep("coding", length(coding_modif0C)+length(coding_modif1C)),
                           rep("promoter", length(prom_modif0C)+length(prom_modif1C)),
                           rep("intergenic", length(inter_modif0C)+length(inter_modif1C)),
                           rep("Unmeth", length(unmethC))),
                    strand=c(rep("IPD0", length(coding_modif0C)),
                             rep("IPD1", length(coding_modif1C)),
                             rep("IPD0", length(prom_modif0C)),
                             rep("IPD1", length(prom_modif1C)),
                             rep("IPD0", length(inter_modif0C)), 
                             rep("IPD1", length(inter_modif1C)),
                             rep("Unmeth", length(unmethC))),
                   ipd=c(coding_modif0C, coding_modif1C, prom_modif0C, prom_modif1C, inter_modif0C, inter_modif1C, unmethC))
  df <- rbind(dfA, dfC)
  df$strain <- strain
  plotlist[[i]] <- df
}

alldf <- do.call(rbind, plotlist)

ggplot(alldf, aes(x=type, y=ipd, fill=strand))+geom_boxplot(outlier.size=.5)+
  xlim(c("promoter", "coding", "intergenic", "Unmeth"))+
  facet_grid(.~base)+
  ylim(c(0,25))+xlab("Region")+ylab("IPD ratio")


cod <- alldf[which(alldf$type=="coding"),]
s <- sample(1:nrow(cod), 10000, replace=F)
cod <- cod[s,]

prom <- alldf[which(alldf$type=="promoter"),]
s <- sample(1:nrow(prom), 10000, replace=F)
prom <- prom[s,]

inter <- alldf[which(alldf$type=="intergenic"),]
s <- sample(1:nrow(inter), 10000, replace=F)
inter <- inter[s,]

unm <- alldf[which(alldf$type=="Unmeth"),]
s <- sample(1:nrow(unm), 10000, replace=F)
unm <- unm[s,]

alldf.s <- rbind(cod, prom, inter, unm) #10000 points subsampled for each type

ggplot(alldf.s, aes(x=type, y=ipd, fill=type))+geom_boxplot(outlier.size=.5, alpha=.5)+
  xlim(c("promoter", "coding", "intergenic", "Unmeth"))+
  facet_grid(.~base)+
  ylim(c(0,15))+xlab("Region")+ylab("IPD ratio")
ggsave("coding_inter_prom_2.pdf", useDingbats=F, width=6, height=4)


Aprom <- alldf.s$ipd[which(alldf.s$base=="A" & alldf.s$type=="promoter")]
Acoding <- alldf.s$ipd[which(alldf.s$base=="A" & alldf.s$type=="coding")]
Aint <- alldf.s$ipd[which(alldf.s$base=="A" & alldf.s$type=="intergenic")]
Aunmet <- alldf.s$ipd[which(alldf.s$base=="A" & alldf.s$type=="Unmeth")]

wilcox.test(Aprom, Acoding)$p.value
# 1.56337e-10
wilcox.test(Aprom, Aint)$p.value
# 6.634248e-05
wilcox.test(Acoding, Aint)$p.value
# 4.681991e-26

Cprom <- alldf.s$ipd[which(alldf.s$base=="C" & alldf.s$type=="promoter")]
Ccoding <- alldf.s$ipd[which(alldf.s$base=="C" & alldf.s$type=="coding")]
Cint <- alldf.s$ipd[which(alldf.s$base=="C" & alldf.s$type=="intergenic")]
Cunmet <- alldf.s$ipd[which(alldf.s$base=="C" & alldf.s$type=="Unmeth")]

wilcox.test(Cprom, Ccoding)$p.value
# 0.6014984
wilcox.test(Cprom, Cint)$p.value
# 0.1555036
wilcox.test(Ccoding, Cint)$p.value
# 0.4144237

