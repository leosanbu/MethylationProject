library(circlize)

## Script to make Circos plots using the circlize package                       ##
##                                                                              ##
## October 2017                                                                 ##
##                                                                              ##
## Leonor Sánchez-Busó (Infection Genomics)                                     ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire   ##

bedfiles <- c("WHO_G_chr_motif11.bed", "WHO_G_chr_motif10.bed", "WHO_G_chr_motif19.bed", "WHO_G_chr_motif14.bed",
              "WHO_G_chr_motif16.bed", "WHO_G_chr_motif24.bed", "WHO_G_chr_motif6.bed", "WHO_G_chr_motif18.bed")

colors <- c("blue4", "brown", "darkgoldenrod4", "darkgreen", 
            "darkorchid4", "orange3", "lightslateblue", "firebrick4")
names(colors) <- c("motif11", "motif10", "motif19", "motif14",
                   "motif16", "motif24", "motif6", "motif18")

## Initialize genome
## Uncomment the strain in use and replace genome size

gsize <- 2167361

# df <- data.frame(name  = c("WHO_F"), start = c(1), end = c(2292467))
df <- data.frame(name  = c("WHO_G"), start = c(1), end = c(gsize)) #chr 2167361 #pCryptic 4207 #pConjugative 42004
# df <- data.frame(name = c("WHO_K_pBlaTEM"), start = c(1), end = c(4153)) #chr 2169846 #pCryptic 4153
# df <- data.frame(name = c("WHO_L_pBlaTEM"), start = c(1), end = c(39054)) #chr 2168633 #pCryptic 4207 #pConjugative 39054
# df <- data.frame(name = c("WHO_M_pBlaTEM"), start = c(1), end = c(5597)) #chr 2178344 #pCryptic 4153 #pConjugative 39015 #pBlaTEM 5597
# df <- data.frame(name = c("WHO_N_pBlaTEM"), start = c(1), end = c(7449)) #chr 2172826 #pCryptic 4207 #pConjugative 42004 #pBlaTEM 7449
# df <- data.frame(name = c("WHO_O_pBlaTEM"), start = c(1), end = c(5598)) #chr 2169062 #pCryptic 4207 #pConjugative 39015 #pBlaTEM 5598
# df <- data.frame(name = c("WHO_P_pBlaTEM"), start = c(1), end = c(4207)) #chr 2173861 #pCryptic 4207
# df <- data.frame(name = c("WHO_U_pBlaTEM"), start = c(1), end = c(4153)) #chr 2234269 #pCryptic 4153
# df <- data.frame(name = c("WHO_V_pBlaTEM"), start = c(1), end = c(7427)) #chr 2221284 #pCryptic 4153 #pBlaTEM 7427
# df <- data.frame(name = c("WHO_W_pBlaTEM"), start = c(1), end = c(39054)) #chr 2222386 #pCryptic 4153 #pConjugative 39054
# df <- data.frame(name = c("WHO_X_pBlaTEM"), start = c(1), end = c(4153)) #chr 2171112 #pCryptic 4153
# df <- data.frame(name = c("WHO_Y_pBlaTEM"), start = c(1), end = c(4153)) #chr 2228980 #pCryptic 4153
# df <- data.frame(name = c("WHO_Z_pBlaTEM"), start = c(1), end = c(4153)) #chr 2153788 #pCryptic 4153
# df <- data.frame(name = c("FA1090"), start = c(1), end = c(2153788))

allbeds <- read.table(bedfiles[1], header=T) 
for (x in 2:length(bedfiles)){
  tmp <- read.table(bedfiles[x], header=T)
  allbeds <- rbind(allbeds, tmp) 
}

features <- read.table("WHO_G_chr_features.bed", header=T)
features <- features[order(features$start),]


pdf("WHO_G_circos.pdf")
circos.genomicInitialize(df)

## Customize genome track
circos.genomicTrack(features, ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                    circos.genomicRect(region,
                                         ytop = 1-0.1,
                                         ybottom = 0+0.1,
                                         border = as.vector(value$clr), col = as.vector(value$clr))
                    #circos.genomicText(region, value, y = 1, labels.column = 2, facing = "clockwise", adj = c(1, 0.5),
                    #                     posTransform = posTransform.text, cex = 0.8, padding = 0.2, niceFacing = TRUE)
                    }, bg.border = "black", track.height = 0.0400005)

cols <- rep("darkblue", nrow(allbeds))
cols[which(allbeds$value!="U")] <- "orangered"
allbeds$clr <- cols
allbeds$clr <- as.factor(allbeds$clr)

n <- length(bedfiles) # number of tracks
circos.genomicTrack(allbeds, ylim = c(0.5, n + 0.5), 
                    panel.fun = function(region, value, ...) {
                      all_tx = unique(value$motif)
                      for(i in seq_along(all_tx)) {
                        l = value$motif == all_tx[i]
                        # for each motif
                        current_tx_start = min(region[l, 1])
                        current_tx_end = max(region[l, 2])
                        circos.lines(c(1,gsize), c(n - i + 0.5, n - i + 0.5), col = "grey")
                        circos.genomicRect(region[l, , drop = FALSE], 
                                           ytop = n - i + 1 + 0.4, 
                                           ybottom = n - i + 1 - 0.4,
                                           border = as.vector(value$clr[l]), col = as.vector(value$clr[l]))
                      }
                    }, bg.border = NA, track.height = 0.25)

dev.off()
 




