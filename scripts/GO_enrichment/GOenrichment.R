args <- commandArgs(TRUE)
suppressMessages(library(topGO))

## Script to perform a GO enrichment of a list of selected geneNames           ##
##                                                                             ##
## September 2017                                                              ##
##                                                                             ##
## Leonor Sánchez-Busó (Infection Genomics)                                    ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire  ##

gene2GO <- readMappings(args[1]) #'pan_genome_reference.interproscan.GO.txt'
geneNames <- read.table(args[2], sep="\t", row.names=1) #'geneNames.txt'
targetGenes <- scan(args[3], what="\n")
outputFileName <- args[4]

for (i in 1:length(gene2GO)){
  code <- names(gene2GO)[i]
  name <- as.vector(geneNames[code,1])
  names(gene2GO)[i] <- name
}

geneList <- factor(as.integer(names(gene2GO) %in% targetGenes))
names(geneList) <- names(gene2GO)

# The 'ontology' argument can be 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
finalRes <- matrix(1:8, 1, 8)
ontology <- c("BP", "MF", "CC")
for (on in ontology){
  # Build topGOdata object
  GOdata <- new("topGOdata", ontology = on, allGenes=geneList, annot = annFUN.gene2GO, gene2GO = gene2GO, nodeSize=10)
  
  # Stats for the GO terms
  GOstats <- termStat(GOdata)
  
  # Run Fisher tests
  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
  resultWeight <- runTest(GOdata, algorithm="weight01", statistic="fisher")
  
  allRes <- GenTable(GOdata, classic = resultFisher, weight = resultWeight,
                     orderBy = "weight", ranksOf = "weight", topNodes = 20)
  allRes <- cbind(rep(on, nrow(allRes)), allRes)
  
  finalRes <- rbind(finalRes, as.matrix(allRes))
}
finalRes <- finalRes[-1,]
colnames(finalRes)[1] <- "Ontology"
rownames(finalRes) <- 1:nrow(finalRes)
finalRes <- as.data.frame(finalRes)

write.table(finalRes, file=outputFileName, quote = F, sep="\t", row.names=F)

