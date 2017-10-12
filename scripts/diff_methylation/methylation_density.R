
## Script to calculate methylation density in chromosomes and plasmids         ##
##                                                                             ##
## October 2017                                                                ##
##                                                                             ##
## Leonor Sánchez-Busó (Infection Genomics)                                    ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire  ##

tab <- read.table("methylation_density.txt", header=T)

## As

pdf("barplot_A_C_density.pdf")
m <- matrix(1:8, 2, 4, byrow = T)
layout(m)
tab.A <- tab[which(tab$Base=="A"),]
rownames(tab.A) <- tab.A$Strain
tab.A <- tab.A[order(tab.A$Chromosome),]
bp <- barplot(tab.A$Chromosome, horiz=TRUE, main="Chromosome", col="lightblue")
text(x=0.1, y=bp[,1], labels = rownames(tab.A))
# tab.Ap1 <- tab.A[which(tab.A$pCryptic>0),]
barplot(tab.A$pCryptic, horiz=TRUE, main="pCryptic", col="lightcoral")
# tab.Ap2 <- tab.A[which(tab.A$pConjugative>0),]
barplot(tab.A$pConjugative, horiz=TRUE, main="pConjugative", col="lightgoldenrod1")
# tab.Ap3 <- tab.A[which(tab.A$pBlaTEM>0),]
barplot(tab.A$pBlaTEM, horiz=TRUE, main="pBlaTEM", col="lightgreen")

tab.C <- tab[which(tab$Base=="C"),]
rownames(tab.C) <- tab.C$Strain
tab.C <- tab.C[order(tab.C$Chromosome),]
bp <- barplot(tab.C$Chromosome, horiz=TRUE, main="Chromosome", col="lightblue")
text(x=0.1, y=bp[,1], labels = rownames(tab.C))
# tab.Cp1 <- tab.C[which(tab.C$pCryptic>0),]
barplot(tab.C$pCryptic, horiz=TRUE, main="pCryptic", col="lightcoral")
# tab.Cp2 <- tab.C[which(tab.C$pConjugative>0),]
barplot(tab.C$pConjugative, horiz=TRUE, main="pConjugative", col="lightgoldenrod1")
# tab.Cp3 <- tab.C[which(tab.C$pBlaTEM>0),]
barplot(tab.C$pBlaTEM, horiz=TRUE, main="pBlaTEM", col="lightgreen")
dev.off()
