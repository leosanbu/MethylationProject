
## Gwide IPD ratios for WHO_F and make plots for particular regions ##

gff <- get_GFF("../GFF_files/WHO_F.gff")
gsize <- 2292467
  
coords <- coords_table(window=5000, step=2500, max=gsize)
if (coords[nrow(coords),2]>gsize){ coords[nrow(coords),2] <- gsize }

modif <- fread("WHO_F/modifications.csv", header=T, sep=",")  # Use 'fread' to load very big tables #
IPD0 <- modif[which(modif$strand==0),]
IPD0 <- IPD0$ipdRatio
IPD1 <- modif[which(modif$strand==1),]
IPD1 <- IPD1$ipdRatio

extrainfo <- read.table("WHO_F_chr_features.bed", header=T)
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
          opa <- grep("opA54|opacity|piiC", tmp_gff[,9], ignore.case = TRUE)
          pil <- grep("pilE|pilS", tmp_gff[,9], ignore.case = TRUE)
          maf <- grep("mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
          if (length(opa)>0){
            co <- "yellowgreen"
          }else if (length(pil)>0){
            co <- "red4"
          }else if (length(maf)>0){
            co <- "dodgerblue"
          }else{
            phages <- grep("phage", tmp_gff[,9])
            if (length(phages)>0){
              co <- "darkmagenta"
            }else{
              rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
              if (length(rib)>0){
                co <- "paleturquoise3"
              }else{
                co <- ""
              }
            }        
          }
        }else{
          wh <- c(max(which(genes_end<=end)), min(which(genes_ini>=ini)))
          tmp_gff <- gff[wh,, drop=F]
          opa <- grep("opA54|opacity|piiC", tmp_gff[,9], ignore.case = TRUE)
          pil <- grep("pilE|pilS", tmp_gff[,9], ignore.case = TRUE)
          maf <- grep("mafA|mafB", tmp_gff[,9], ignore.case = TRUE)
          if (length(opa)>0){
            co <- "yellowgreen"
          }else if (length(pil)>0){
            co <- "red4"
          }else if (length(maf)>0){
            co <- "dodgerblue"
          }else{
            phages <- grep("phage", tmp_gff[,9])
            if (length(phages)>0){
              co <- "darkmagenta"
            }else{
              rib <- grep("30S ribosomal protein|50S ribosomal protein", tmp_gff[,9], ignore.case = TRUE)
              if (length(rib)>0){
                co <- "paleturquoise3"
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
  }else{ 
    # no extra info, check GFF
  }
}

getMean <- apply(coords, 1, function(x){
  zero <- mean(IPD0[x[1]:x[2]])
  one <- mean(IPD1[x[1]:x[2]])
  list(zero,one)
})

IPD0win <- unlist(lapply(getMean, function(x) x[[1]]))
IPD1win <- unlist(lapply(getMean, function(x) x[[2]]))

df <- data.frame(start=coords[,1], end=coords[,2], IPD0=IPD0win, IPD1=IPD1win, col=col)
df$col <- as.vector(df$col)
df$col[which(df$col=="")] <- "grey80"

cols <- c("darkmagenta", "darkorange2", "dodgerblue", "goldenrod1", "paleturquoise3", 
          "grey80", "red4", "yellowgreen")
names(cols) <- c("darkmagenta", "darkorange2", "dodgerblue", "goldenrod1", "paleturquoise3", 
                 "grey80", "red4", "yellowgreen")

dfm <- melt(df[,3:5])
sortmedian <- aggregate(value ~ col, data=dfm, FUN=median)
sortmedian <- sortmedian[order(sortmedian$value),]
  
ggplot(dfm, aes(x=col, y=value, fill=col))+
  geom_boxplot(alpha=.5)+
  scale_fill_manual(values=cols)+
  scale_x_discrete(limits=sortmedian[,1], 
                   labels=c("GGI", "mafA/mafB", "GI", "ribosomal", "pilE/pilS",
                            "phage-related", "opA54/piiC", "rest"))+
  xlab("Genomic region")+
  ylab("IPD ratio")
ggsave("IPD_by_category.pdf", useDingbats=F, width=6, height=5)

restvalues <- dfm$value[which(dfm$col=="grey80")]
co <- cols[which(cols!="grey80")]
pvals <- c()
for (i in 1:length(co)){
  sub <- dfm[which(dfm$col==co[i]),3]
  test <- wilcox.test(sub, restvalues)$p.value
  pvals <- c(pvals, test)
}
names(pvals) <- co

# pvals
# darkmagenta    darkorange2     dodgerblue     goldenrod1 paleturquoise3           red4    yellowgreen 
# 1.014354e-02   4.015098e-26   4.762410e-24   5.088369e-28   8.788221e-10   4.159504e-02   8.183620e-01 
# phage-related  GI             mafA/mafB       GGI           ribosomal      pilE/pilS      opA/piiC

#GGI
pdf("GGI_WHO_F.pdf")
m <- matrix(1:2, 2, 1)
layout(m)
par(mar=c(0,4,12,4))
plot(df$IPD0[291:356], type="h", lwd=6, col=df$col[291:356], xaxt="n", xlab="", ylab="Mean IPD0")
par(mar=c(12,4,0,4))
plot(-df$IPD1[291:356], type="h", lwd=6, col=df$col[291:356], xaxt="n", xlab="", ylab="Mean IPD1")
axis(1, at=seq(1, 66, by=10), labels=df$start[291:356][seq(1, 66, by=10)], las=2)
dev.off()

pdf("GI1_WHO_F.pdf")
m <- matrix(1:2, 2, 1)
layout(m)
par(mar=c(0,4,12,4))
plot(df$IPD0[387:441], type="h", lwd=7, col=df$col[387:441], xaxt="n", xlab="", ylab="Mean IPD0")
par(mar=c(12,4,0,4))
plot(-df$IPD1[387:441], type="h", lwd=7, col=df$col[387:441], xaxt="n", xlab="", ylab="Mean IPD1")
axis(1, at=seq(1,55, by=10), labels=df$start[387:441][seq(1,55, by=10)], las=2)
dev.off()

pdf("GI2_WHO_F.pdf")
m <- matrix(1:2, 2, 1)
layout(m)
par(mar=c(0,4,12,4))
plot(df$IPD0[595:648], type="h", lwd=7, col=df$col[595:648], xaxt="n", xlab="", ylab="Mean IPD0")
par(mar=c(12,4,0,4))
plot(-df$IPD1[595:648], type="h", lwd=7, col=df$col[595:648], xaxt="n", xlab="", ylab="Mean IPD1")
axis(1, at=seq(1,54, by=10), labels=df$start[595:648][seq(1,54, by=10)], las=2)
dev.off()

pdf("maf_phage_WHO_F.pdf")
m <- matrix(1:2, 2, 1)
layout(m)
par(mar=c(0,4,12,4))
plot(df$IPD0[427:494], type="h", lwd=6, col=df$col[427:494], xaxt="n", xlab="", ylab="Mean IPD0")
par(mar=c(12,4,0,4))
plot(-df$IPD1[427:494], type="h", lwd=6, col=df$col[427:494], xaxt="n", xlab="", ylab="Mean IPD1")
axis(1, at=seq(1,68, by=10), labels=df$start[427:494][seq(1,68, by=10)], las=2)
dev.off()

pdf("ribosomal_prots_WHO_F.pdf")
m <- matrix(1:2, 2, 1)
layout(m)
par(mar=c(0,4,12,4))
plot(df$IPD0[748:797], type="h", lwd=7, col=df$col[748:797], xaxt="n", xlab="", ylab="Mean IPD0")
par(mar=c(12,4,0,4))
plot(-df$IPD1[748:797], type="h", lwd=7, col=df$col[748:797], xaxt="n", xlab="", ylab="Mean IPD1")
axis(1, at=seq(1,50,by=10), labels=df$start[748:797][seq(1,50,by=10)], las=2)
dev.off()
