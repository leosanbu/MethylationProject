args <- commandArgs(TRUE)
suppressMessages(library(Rcpp))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))

## Script to extract the IPD values for each base of a motif from a list of motifs ##
##                                                                                 ##
## July 2017                                                                       ##
##                                                                                 ##
## Leonor Sánchez-Busó (Infection Genomics)                                        ##
## Wellcome Sanger Institute, Wellcome Genome Campus, Hinxton, Cambridgeshire      ##

## Motifs coordinates are expected to be extracted using 'fuzznuc' (EMBOSS), i.e. 
## fuzznuc -sequence $i.fna -pattern @ngono_motifs.pat -outfile $i.fuzznuc -complement 1     
## http://emboss.sourceforge.net/apps/cvs/emboss/apps/fuzznuc.html                

## Run this script using:
## Rscript methylation_queryMotif.R path/to/fuzznuc/files paths_to_modifications.txt query_motifs.txt 
##
## IMPORTANT: order of arguments matters!
## path/to/fuzznuc/files:      path to where fuzznuc output files are
## paths_to_modifications.txt: table containing name<TAB>path_to_modifications.csv_file
## query_motifs.txt:           list of motifs to work with

# Get input files #
fuzzpath <- args[1]
pathfile <- args[2]
motifslist <- args[3]

### Functions ###

# Get IPDs - Rcpp #
cppFunction('List getIPDs(NumericMatrix x, NumericMatrix y){
            int nrowx = x.nrow();
            int nrowy = y.nrow();
            NumericVector indices = y(_,2);
            List IPD0list(nrowx);
            List IPD1list(nrowx);
            
            for (int z = 0; z < nrowx; z++){
            int start = x(z,0);
            int end = x(z,1);
            int n = end-start+1;
            
            NumericVector IPD0(n);
            NumericVector IPD1(n);
            int c0 = 0;
            int c1 = 0;
            for (int w = 0; w < nrowy; w++){
            if ((start <= y(w,0)) && (y(w,0) <= end)){
            if (y(w,1)==0){ 
            IPD0[c0] = y(w,2); 
            c0 += 1;
            }else{ 
            IPD1[c1] = y(w,2); 
            c1 += 1;
            }
            }
            }
            IPD0list[z] = IPD0;
            IPD1list[z] = IPD1;
            }
            
            List IPD;
            IPD["IPD0"] = IPD0list;
            IPD["IPD1"] = IPD1list;
            return IPD;
            }')


readMotifs <- function(motifslist){
  queryFile <- scan(motifslist, what="\n")
  queryMat <- matrix(1:2, 1, 2)
  name <- ""
  seq <- ""
  for (i in 1:length(queryFile)){
    g <- grep("^>", queryFile[i])
    if (length(g)>0){
      if (seq != ""){
        queryMat <- rbind(queryMat, c(name, seq))
      }
      name = queryFile[i]
      seq = ""
    }else{
      if (name != ""){
        seq = queryFile[i]
      }
    }
  }
  queryMat <- rbind(queryMat, c(name, seq))
  queryMat <- queryMat[-1,]
  queryMat[,1] <- gsub(">", "", queryMat[,1])
  return(queryMat)
}
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
  # g <- grep("contig", fuzz[gseq])
  # if (length(g)==0){
  #   names.gseq <- c()
  #   for (x in 1:length(gseq)){
  #     names.gseq <- c(names.gseq, paste0("contig00000", x))
  #   }
  # }else{
    names.gseq <- gsub("# Sequence: ", "", fuzz[gseq])
    names.gseq <- gsub(".*contig", "contig", names.gseq)
    names.gseq <- gsub(" .*", "", names.gseq)
    names.gseq <- gsub("\\..*", "", names.gseq)
  # }
  names(fuzzlist) <- names.gseq
  return(fuzzlist)
}
mainFunction <- function(motifs, mod.file, motif){
  
  mod <- mod.file
  queryMotif <- motif
  
  # Make input matrices numeric #
  motif.num <- cbind(as.numeric(motifs[[1]]), as.numeric(motifs[[2]]))
  modif.num <- cbind(as.numeric(mod[[2]]),as.numeric(mod[[3]]),as.numeric(mod[[9]]))
  
  # Get IPDs using Rcpp function #
  system.time(IpdScore <- getIPDs(motif.num,modif.num))
  
  IPD0 <- do.call(rbind, IpdScore[[1]])
  IPD1 <- do.call(rbind, IpdScore[[2]])
  colnames(IPD0) <- strsplit(queryMotif, "")[[1]]
  colnames(IPD1) <- strsplit(queryMotif, "")[[1]]
  
  IPD0 <- cbind(rep("IPD0", nrow(motifs)), motifs, IPD0)
  IPD1 <- cbind(rep("IPD1", nrow(motifs)), motifs, IPD1)
  colnames(IPD0)[1:7] <- c("IPDstrand", "MotifStart", "MotifEnd", "MotifStrand", "MotifName", ".", "Motif")
  colnames(IPD1)[1:7] <- c("IPDstrand", "MotifStart", "MotifEnd", "MotifStrand", "MotifName", ".", "Motif")
  return(list(IPD0, IPD1))
}

## MAIN ##
fuzznucfiles <- dir(fuzzpath, pattern=".fuzznuc$")
paths <- read.table(pathfile)
queryMat <- readMotifs(motifslist)

for (i in 1:length(fuzznucfiles)){
  cat(fuzznucfiles[i], "\n")
  
  # Load modifications file #
  strain <- gsub(fuzznucfiles[i], pattern="\\..*", replacement="")
  g <- which(paths[,1]==strain)
  file <- as.character(paths[g,2])
  if (!file.exists(file)){
    zip <- grep(".gz", file)
    if (length(zip)>0){
      system(paste("gunzip", file))
      file <- gsub(paths[g], pattern=".gz", replacement="")      
    }
  }
  
  # Read modifications files from PacBio #    
  modifications <- fread(file, header=T, sep=",")  # Use 'fread' to load very big tables #
  modifications[[1]] <- gsub(" .*", "", modifications[[1]])
  modifications[[1]] <- gsub(".*\\|", "", modifications[[1]])
  
  # Parse fuzznuc file #
  fuzznuc <- parseFuzznuc(fuzznucfiles[i])
  
  # Check how many contigs! #
  contigs <- names(fuzznuc)#unique(modifications[[1]])
  
  # Run main function #
  for (x in 1:nrow(queryMat)){
    motifName <- queryMat[x,1]
    motifSeq <- queryMat[x,2]
      
    # Create a directory if it doesn't exist
    if (!file.exists(motifName)){
      system(paste("mkdir", motifName))
    }
    
    # Run only if file doesn't exist
    check.file <- paste0(motifName, "/", strain, "_", motifName, "_IPDs.txt")
    if (!file.exists(check.file)){
      point <- 0
      for (con in 1:length(contigs)){
        contig <- gsub(".*contig", "contig", contigs[con])
        contig <- gsub(" .*", "", contig)
        mod <- modifications[which(modifications[[1]]==contigs[con]),] # Subset the modifications file to the contig #
        #g <- grep(contig, modifications[[1]])
        #if (length(g)>0){ mod <- modifications[g,] }else{ stop("Contig ", contig, " not found in modifications.csv file\n") }
        fuzztig <- fuzznuc[[contig]] # Subset the fuzznuc file to the contig #
        motifs <- fuzztig[which(fuzztig[[4]]==motifName),] # Select motif #
        if (is.null(dim(motifs))){ motifs <- t(as.matrix(motifs)) }
        
        if (nrow(motifs)>0){
          point <- point+1
          tmpResults <- mainFunction(motifs, mod, motifSeq)
          if (point==1){ 
            printOut <- rbind(tmpResults[[1]], tmpResults[[2]])
            printOut <- cbind(rep(contig, nrow(printOut)), printOut)
          }else{
            tmpprint <- rbind(tmpResults[[1]], tmpResults[[2]])
            tmpprint <- cbind(rep(contig, nrow(tmpprint)), tmpprint)
            printOut <- rbind(printOut, tmpprint)
          }
        }
      }
      
      # Write ipd ratios of this motif/strain to a text file #
      colnames(printOut)[1] <- "contig"
      
      # For motifs in reverse strand, reverse values (move this bit into previous scripts)
      g0 <- which(printOut$MotifStrand=="-" & printOut$IPDstrand=="IPD0")
      g1 <- which(printOut$MotifStrand=="-" & printOut$IPDstrand=="IPD1")
      
      if (length(g0)>0){
        printOut$IPDstrand[g0] <- "IPD1"
        printOut[g0,9:ncol(printOut)] <- rev(printOut[g0,9:ncol(printOut)])
      }
      if (length(g1)>0){ 
        printOut$IPDstrand[g1] <- "IPD0" 
        printOut[g1,9:ncol(printOut)] <- rev(printOut[g1,9:ncol(printOut)])
      }
      
      write.table(printOut, paste0(motifName, "/", strain, "_", motifName, "_IPDs.txt"), sep="\t", quote = F, row.names=F) 
      
    }
  }
}
