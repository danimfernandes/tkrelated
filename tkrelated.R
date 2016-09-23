## Created by: Daniel Fernandes - dani.mag.fernandes@gmail.com
## Presented on:
## "The Identification of a 1916 Irish Rebel: New Approach for Estimating Relatedness From Low Coverage Homozygous Genomes"
## 2016

##################### CybRsex- SIMULATING RELATEDNESS #############################
#---------------------------------------------------------------------------------#
# This script generates random sets of related and unrelated individuals based on #
# SNP allele frequencies. These must be provided in PLINK format.		              #
# Using the third-party software SPAGeDI, relatedness coefficients for these      #
# individuals are calculated and then plotted as histograms.			                #
#											                                                            #
# INPUT:										                                                      #
# plink.frq -allele frequencies in PLINK tableformat	                  			    #

# Set working directory
setwd("foo/bar")

######################## SET OF UNRELATED INDIVIDUALS #######################
# --------------------------------------------------------------------------#
unrelatedfunc = function(file,samples.tag,numUnrelatedPairs,identifier,run.spagedi=TRUE,reduce.SNPs=FALSE) {

  ## Arguments' description:
  ##
  ## file - name/location of the plink frequencies file to simulate individuals from
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## numUnrelatedPairs - the number of unrelated pairs of individuals to generate
  ## identifier - a numeric identifier, for example, 1. Can be used to generate many independent sets who can be joined together at the end
  ## run.spagedi - boolean for whether to automatically run SPAGeDI or not
  ## reduce.SNPs - if your SNP dataset is too large it might crash SPAGeDI, so you can provide a number of SNPs to randomly subset from the original file


  ################# Read allele frequencies
  alFreq = read.csv(file,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  #A1 is minor allele
  alFreq$NCHROBS=NULL
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ##Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
  maf=format(as.numeric(alFreq$MAF),signif=3)
  alFreq$MAF=maf
  af2=format((1-as.numeric(maf)),signif=3)
  alFreq$AFa2=af2
  alFreq$CHR=NULL
  
  if(typeof(reduce.SNPs)=="double") {
    alFrequ=alFreq[sample(1:nrow(alFreq),reduce.SNPs,replace=FALSE),]
  }  else { 
    alFrequ=alFreq 
  }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  numIndividuals = numUnrelatedPairs*2
  seqIndividuals = seq(1:numIndividuals)
  alFrequ=alFrequNumericals
  for(i in seqIndividuals) {
    alFrequ[paste("IND",i,"a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
    alFrequ[paste("IND",i,"b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
  }

  ################# Create new data.frame with forced homozygous
  alFrequFH=alFrequ[1:5]
  for(i in seqIndividuals){
    alFrequFH[paste0("IND",i)] = ifelse(as.character(runif(length(alFrequ$SNP),0.00,1.00)) < 0.5, c(alFrequ[paste("IND",i,"a",sep="")])[[1]], c(alFrequ[paste("IND",i,"b",sep="")])[[1]])
  }
  alFrequFH$SNP=NULL
  alFrequFH$A1=NULL
  alFrequFH$A2=NULL
  alFrequFH$MAF=NULL
  alFrequFH$AFa2=NULL

  ################# Create and export SPAGEDI's data file
  ## Remove previous SPAGeDI output files 
  for(iterit in list.files(getwd())) {
    if(grepl(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,"_spgout.txt"),iterit)== TRUE) {
      system(paste0("rm -f ",iterit))
    }
  }
  ## Create input file
  maxCols=length(alFrequ$SNP)
  sink(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,".txt"))
  row1=c(numIndividuals,0,0,maxCols,3,2)
  cat(row1,sep="\t")
  cat("\n")
  row2="0"
  cat(row2,sep="\t")
  cat("\n")
  row3=c("IND",alFrequ$SNP)
  cat(row3,sep="\t")
  cat("\n")
  indList=paste0("IND",seqIndividuals)
  iter=4  ## First row in SPAGEDI input file for samples
  for(i in indList){
    colNum=grep(paste0("^",i,"$"),c(names(alFrequFH)))
    cat(c(i,paste(alFrequFH[,colNum],alFrequFH[,colNum],sep="")),sep="\t")
    cat("\n")
  }
  cat("END")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  row1=c(paste(alFrequ$SNP,"2",sep="\t"))
  row2=c(paste(alFrequ$A1,alFrequ$MAF,sep="\t"))
  row3=c(paste(alFrequ$A2,alFrequ$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_freq",identifier,".txt"),append=FALSE)
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  sink(paste0("SP_CMDS_Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,".txt"))
  cat(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,".txt\n"))
  cat(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,"_spgout.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0("Unrelated",numUnrelatedPairs,"_",samples.tag,"_freq",identifier,".txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,".txt" ))
    cmd = paste0("spagedi < SP_CMDS_Unrelated",numUnrelatedPairs,"_",samples.tag,"_in",identifier,".txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}
########### Function to read output from SPAGeDI and plot it according to provided samples.tag
unrelatedplot = function(samples.tag) {
  wdfiles=list.files(getwd())
  it=1
  for(filet in wdfiles) {
    if(grepl("_spgout.txt$",filet) == TRUE) {
      if(grepl("^Unrelated",filet) == TRUE) {
        if(grepl(paste0("_",samples.tag,"_"),filet) == TRUE) {
          print(paste0("Reading SPAGeDi output file: ",filet))
          assign(paste0("histoUN","_",samples.tag,"_",it),read.csv(filet,header = FALSE))
          it=it+1
        }
      }
    }
  }
  
  ## Read until: "PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form"
  unrelated = c()
  enviro=ls()
  for(i in enviro) {
    if(grepl(paste0("^histoUN","_",samples.tag,"_"),i)==TRUE){
      it=1
      print(paste0("Reading genetic distances from: ",i))
      lineInd= which(grepl(" individuals$",eval(parse(text=paste0(i,"[,1]"))))==TRUE)
      numInd = as.numeric(strsplit(as.character(eval(parse(text=paste0(i,"[",lineInd,",1]"))))," ")[[1]][1])
      lineStart = as.numeric(which("PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form" == eval(parse(text=paste0(i,"[,1]")))) + 3)
      lineEnd = length(eval(parse(text=paste0(i,"[,1]"))))
      assign(paste0("Dist",i),eval(parse(text=paste0("data.frame(",i,"[",lineStart,":",lineEnd,",1],stringsAsFactors=FALSE)"))))
      write.table(eval(parse(text=paste0("Dist",i))),paste0(i,"Distances"),quote = FALSE,row.names = FALSE,col.names = FALSE)
      assign(paste0("Distances_in_",i),read.table(paste0(i,"Distances"),header = FALSE))
      tempdf = get(paste0("Distances_in_",i))
      tempdf[,7] = paste(tempdf[,3],tempdf[,4],sep="-")
      ## The to_remove loop works on the supposition that the pair names are like IND3-4I and IND3-4II, with the full name of the first pair
      ## included in the second pair (IND3-4I exists in IND3-4II)
      it=1
      to_keep = c("1")
      for(x in seq(1:(numInd/2))) {
        it=it+2
        keepPos = which(paste(it,it+1,sep="-") == tempdf[,7])
        to_keep = c(to_keep,keepPos)
      }
      posis=seq(1:length(tempdf[,1]))
      to_remove = setdiff(posis, as.numeric(to_keep))
      tempdf=tempdf[-to_remove,]
      assign(paste0("Distances_in_",i),tempdf)
      unrelated=c(unrelated,tempdf[,6])
    }
  }
  
  write.table(unrelated,paste0("histoUN_",samples.tag,"Distances_Final"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  xlimitsFS=c(min(unrelated)-0.05,max(unrelated)+0.05)
  hist(unrelated,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsFS,ylab=paste0("N=",length(unrelated)), main=paste0("Unrelated individuals for ",samples.tag))
}


############################ SET OF FULL SIBLINGS ###########################
# ------------------------------------------------------------------------- #
fullsibsfunc = function(file,samples.tag,numFullSiblingsPairs,identifier,run.spagedi=TRUE,reduce.SNPs=FALSE) {
  
  ## Arguments' description:
  ##
  ## file - name/location of the plink frequencies file to simulate individuals from
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## numFullSiblingPairs - the number of full-sibling pairs of individuals to generate
  ## identifier - a numeric identifier, for example, 1. Can be used to generate many independent sets who can be joined together at the end
  ## run.spagedi - boolean for whether to automatically run SPAGeDI or not
  ## reduce.SNPs - if your SNP dataset is too large it might crash SPAGeDI, so you can provide a number of SNPs to randomly subset from the original file

  ################# Read allele frequencies
  alFreq = read.csv(file,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  #A1 is minor allele
  alFreq$NCHROBS=NULL
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ##Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
  maf=format(as.numeric(alFreq$MAF),signif=3)
  alFreq$MAF=maf
  af2=format((1-as.numeric(maf)),signif=3)
  alFreq$AFa2=af2
  alFreq$CHR=NULL
  
  if(typeof(reduce.SNPs)=="double") {
    alFrequ=alFreq[sample(1:nrow(alFreq),reduce.SNPs,replace=FALSE),]
  }
  else {
    alFrequ=alFreq 
  }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  numIndividuals = numFullSiblingsPairs*2 
  seqIndividuals = seq(1:numIndividuals)
  alFrequ=alFrequNumericals
  for(i in seqIndividuals) {
    alFrequ[paste("IND",i,"a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
    alFrequ[paste("IND",i,"b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
  }
  
  ################# Choosing one allele at random from parent
  fullSib=alFrequ[1:5]
  iter=1
  for(i in seqIndividuals){
    nextInd=iter+1
    if(nextInd <= numIndividuals) {
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"b",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"b",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
    }
    iter=iter+2
  }
  fullSib$SNP=NULL
  fullSib$A1=NULL
  fullSib$A2=NULL
  fullSib$MAF=NULL
  fullSib$AFa2=NULL
  
  ################### Create a list of the pairs of siblings
  listFullSibs = list()
  numSiblingPairs = numFullSiblingsPairs
  seqSiblings = seq(1:numSiblingPairs)
  iter=1
  for(i in seqSiblings) {
    nextInd=iter+1
    ind = paste0("FSIB",iter,"-",nextInd,"I") 
    ind2 = paste0("FSIB",iter,"-",nextInd,"II") 
    listFullSibs = c(listFullSibs, ind,ind2)
    iter=iter+2
  }
  
  ################# Force homozygous
  fullSibFH=alFrequ[1:5]
  for(i in listFullSibs){
    try((fullSibFH[paste0(i)] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(fullSib[paste(i,"_a",sep="")])[[1]], c(fullSib[paste(i,"_b",sep="")])[[1]])),silent=TRUE)
  }
  fullSibFH$SNP=NULL
  fullSibFH$A1=NULL
  fullSibFH$A2=NULL
  fullSibFH$MAF=NULL
  fullSibFH$AFa2=NULL
  
  ################# Create and export SPAGEDI's data file
  ## Remove previous SPAGeDI output files 
  for(iterit in list.files(getwd())) {
    if(grepl(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt"),iterit)== TRUE) {
      system(paste0("rm -f ",iterit))
    }
  }
  ## Create input file
  maxCols=length(alFrequ$SNP)
  sink(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  row1=c(length(listFullSibs),0,0,maxCols,3,2)
  cat(row1,sep="\t")
  cat("\n")
  row2="0"
  cat(row2,sep="\t")
  cat("\n")
  row3=c("IND",alFrequ$SNP)
  cat(row3,sep="\t")
  cat("\n")
  iter=4  ## First row in SPAGEDI input file for samples
  for(i in listFullSibs){
    colNum=grep(paste0("^",i,"$"),c(names(fullSibFH)))
    cat(c(i,paste(fullSibFH[,colNum],fullSibFH[,colNum],sep="")),sep="\t")
    cat("\n")
  }
  cat("END")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  row1=c(paste(alFrequ$SNP,"2",sep="\t"))
  row2=c(paste(alFrequ$A1,alFrequ$MAF,sep="\t"))
  row3=c(paste(alFrequ$A2,alFrequ$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"),append=FALSE)
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  sink(paste0("SP_CMDS_FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  cat(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,".txt\n"))
  cat(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0("FullSibs",numFullSiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,".txt" ))
    cmd = paste0("spagedi < SP_CMDS_FullSibs",numFullSiblingsPairs,"_",samples.tag,"_in",identifier,".txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}
########### Function to read output from SPAGeDI and plot it according to provided samples.tag
fullsibsplot = function(samples.tag) {
  wdfiles=list.files(getwd())
  it=1
  for(filet in wdfiles) {
    if(grepl("_spgout.txt$",filet) == TRUE) {
      if(grepl("^FullSibs",filet) == TRUE) {
        if(grepl(paste0("_",samples.tag,"_"),filet) == TRUE) {
          print(paste0("Reading SPAGeDi output file: ",filet))
          assign(paste0("histoFS","_",samples.tag,"_",it),read.csv(filet,header = FALSE))
          it=it+1
        }
      }
    }
  }
  
  ## Read until: "PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form"
  fullsibs = c()
  enviro=ls()
  for(i in enviro) {
    if(grepl(paste0("^histoFS","_",samples.tag,"_"),i)==TRUE){
      it=1
      print(paste0("Reading genetic distances from: ",i))
      lineInd= which(grepl(" individuals$",eval(parse(text=paste0(i,"[,1]"))))==TRUE)
      numInd = as.numeric(strsplit(as.character(eval(parse(text=paste0(i,"[",lineInd,",1]"))))," ")[[1]][1])
      lineStart = as.numeric(which("PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form" == eval(parse(text=paste0(i,"[,1]")))) + 3)
      lineEnd = length(eval(parse(text=paste0(i,"[,1]"))))
      assign(paste0("Dist",i),eval(parse(text=paste0("data.frame(",i,"[",lineStart,":",lineEnd,",1],stringsAsFactors=FALSE)"))))
      write.table(eval(parse(text=paste0("Dist",i))),paste0(i,"Distances"),quote = FALSE,row.names = FALSE,col.names = FALSE)
      assign(paste0("Distances_in_",i),read.table(paste0(i,"Distances"),header = FALSE))
      tempdf = get(paste0("Distances_in_",i))
      tempdf[,7] = paste(tempdf[,3],tempdf[,4],sep="-")
      ## The to_remove loop works on the supposition that the pair names are like IND3-4I and IND3-4II, with the full name of the first pair
      ## included in the second pair (IND3-4I exists in IND3-4II)
      it=1
      to_keep = c("1")
      for(x in seq(1:(numInd/2))) {
        it=it+2
        keepPos = which(paste(it,it+1,sep="-") == tempdf[,7])
        to_keep = c(to_keep,keepPos)
      }
      posis=seq(1:length(tempdf[,1]))
      to_remove = setdiff(posis, as.numeric(to_keep))
      tempdf=tempdf[-to_remove,]
      assign(paste0("Distances_in_",i),tempdf)
      fullsibs=c(fullsibs,tempdf[,6])
    }
  }
  
  write.table(fullsibs,paste0("histoFS_",samples.tag,"Distances_Final"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  xlimitsFS=c(min(fullsibs)-0.05,max(fullsibs)+0.05)
  hist(fullsibs,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsFS,ylab=paste0("N=",length(fullsibs)), main=paste0("Full Siblings for ",samples.tag))
}


######################## SET FOR HALF SIBLINGS #######################
# ------------------------------------------------------------------------- #
halfsibsfunc = function(file,samples.tag,numHalfSiblingsPairs,identifier,run.spagedi=TRUE,reduce.SNPs=FALSE) {

  ## Arguments' description:
  ##
  ## file - name/location of the plink frequencies file to simulate individuals from
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## numHalfSiblingPairs - the number of half-sibling pairs of individuals to generate
  ## identifier - a numeric identifier, for example, 1. Can be used to generate many independent sets who can be joined together at the end
  ## run.spagedi - boolean for whether to automatically run SPAGeDI or not
  ## reduce.SNPs - if your SNP dataset is too large it might crash SPAGeDI, so you can provide a number of SNPs to randomly subset from the original file

  ################# Read allele frequencies
  alFreq = read.csv(file,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  ## A1 is minor allele
  alFreq$NCHROBS=NULL
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ## Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
  maf=format(as.numeric(alFreq$MAF),signif=3)
  alFreq$MAF=maf
  af2=format((1-as.numeric(maf)),signif=3)
  alFreq$AFa2=af2
  alFreq$CHR=NULL
  
  if(typeof(reduce.SNPs)=="double") {
    alFrequ=alFreq[sample(1:nrow(alFreq),reduce.SNPs,replace=FALSE),]
  }
  else {
    alFrequ=alFreq 
  }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  numIndividuals = numHalfSiblingsPairs*3 
  numIndividuals = numIndividuals + 1 ##Correcting for permutations
  seqIndividuals = seq(1:numIndividuals)
  for(i in seqIndividuals) {
    alFrequ[paste("IND",i,"a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
    alFrequ[paste("IND",i,"b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
  }
  
  ################# Choosing one allele at random from parent
  halfSib=alFrequ[1:5]
  iter=1
  for(i in seqIndividuals){
    nextInd=iter+1
    nextNextInd=nextInd+1
    if(nextInd <= numIndividuals) {
      halfSib[paste("HSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"b",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextNextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextNextInd,"a",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextNextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"b",sep="")])[[1]], c(alFrequ[paste("IND",nextNextInd,"b",sep="")])[[1]])
    }
    iter=iter+2
  }
  halfSib$SNP=NULL
  halfSib$A1=NULL
  halfSib$A2=NULL
  halfSib$MAF=NULL
  halfSib$AFa2=NULL
  
  ################### Create a list of the pairs of siblings
  listHalfSibs = list()
  numSiblingPairs = numHalfSiblingsPairs
  seqSiblings = seq(1:numSiblingPairs)
  iter=1
  for(i in seqSiblings) {
    nextInd=iter+1
    nextNextInd=nextInd+1
    ind = paste0("HSIB",iter,"-",nextInd,"I") 
    ind2 = paste0("HSIB",iter,"-",nextNextInd,"II") 
    listHalfSibs = c(listHalfSibs, ind,ind2)
    iter=iter+2
  }
  
  ######### Force homozygous
  halfSibFH=alFrequ[1:5]
  for(i in listHalfSibs){
    try((halfSibFH[paste0(i)] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(halfSib[paste(i,"_a",sep="")])[[1]], c(halfSib[paste(i,"_b",sep="")])[[1]])),silent=TRUE)
  }
  halfSibFH$SNP=NULL
  halfSibFH$A1=NULL
  halfSibFH$A2=NULL
  halfSibFH$MAF=NULL
  halfSibFH$AFa2=NULL
  
  ################# Create and export SPAGEDI's data file
  ## Remove previous SPAGeDI output files 
  for(iterit in list.files(getwd())) {
    if(grepl(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt"),iterit)== TRUE) {
      system(paste0("rm -f ",iterit))
    }
  }
  ## Create input file
  maxCols=length(alFrequ$SNP)
  sink(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  row1=c(length(listHalfSibs),0,0,maxCols,3,2)
  cat(row1,sep="\t")
  cat("\n")
  row2="0"
  cat(row2,sep="\t")
  cat("\n")
  row3=c("IND",alFrequ$SNP)
  cat(row3,sep="\t")
  cat("\n")
  iter=4  ## First row in SPAGEDI input file for samples
  for(i in listHalfSibs){
    colNum=grep(paste0("^",i,"$"),c(names(halfSibFH)))
    cat(c(i,paste(halfSibFH[,colNum],halfSibFH[,colNum],sep="")),sep="\t")
    cat("\n")
  }
  cat("END")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  row1=c(paste(alFrequ$SNP,"2",sep="\t"))
  row2=c(paste(alFrequ$A1,alFrequ$MAF,sep="\t"))
  row3=c(paste(alFrequ$A2,alFrequ$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"),append=FALSE)
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  sink(paste0("SP_CMDS_HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  cat(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,".txt\n"))
  cat(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0("HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,".txt" ))
    cmd = paste0("spagedi < SP_CMDS_HalfSibs",numHalfSiblingsPairs,"_",samples.tag,"_in",identifier,".txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}
########### Function to read output from SPAGeDI and plot it according to provided samples.tag
halfsibsplot = function(samples.tag) {
  wdfiles=list.files(getwd())
  it=1
  for(filet in wdfiles) {
    if(grepl("_spgout.txt$",filet) == TRUE) {
      if(grepl("^HalfSibs",filet) == TRUE) {
        if(grepl(paste0("_",samples.tag,"_"),filet) == TRUE) {
          print(paste0("Reading SPAGeDi output file: ",filet))
          assign(paste0("histoHS","_",samples.tag,"_",it),read.csv(filet,header = FALSE))
          it=it+1
        }
      }
    }
  }
  
  ## Read until: "PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form"
  halfsibs = c()
  enviro=ls()
  for(i in enviro) {
    if(grepl(paste0("^histoHS","_",samples.tag,"_"),i)==TRUE){
      it=1
      print(paste0("Reading genetic distances from: ",i))
      lineInd= which(grepl(" individuals$",eval(parse(text=paste0(i,"[,1]"))))==TRUE)
      numInd = as.numeric(strsplit(as.character(eval(parse(text=paste0(i,"[",lineInd,",1]"))))," ")[[1]][1])
      lineStart = as.numeric(which("PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form" == eval(parse(text=paste0(i,"[,1]")))) + 3)
      lineEnd = length(eval(parse(text=paste0(i,"[,1]"))))
      assign(paste0("Dist",i),eval(parse(text=paste0("data.frame(",i,"[",lineStart,":",lineEnd,",1],stringsAsFactors=FALSE)"))))
      write.table(eval(parse(text=paste0("Dist",i))),paste0(i,"Distances"),quote = FALSE,row.names = FALSE,col.names = FALSE)
      assign(paste0("Distances_in_",i),read.table(paste0(i,"Distances"),header = FALSE))
      tempdf = get(paste0("Distances_in_",i))
      tempdf[,7] = paste(tempdf[,3],tempdf[,4],sep="-")
      ## The to_remove loop works on the supposition that the pair names are like IND3-4I and IND3-4II, with the half name of the first pair
      ## included in the second pair (IND3-4I exists in IND3-4II)
      it=1
      to_keep = c("1")
      for(x in seq(1:(numInd/2))) {
        it=it+2
        keepPos = which(paste(it,it+1,sep="-") == tempdf[,7])
        to_keep = c(to_keep,keepPos)
      }
      posis=seq(1:length(tempdf[,1]))
      to_remove = setdiff(posis, as.numeric(to_keep))
      tempdf=tempdf[-to_remove,]
      assign(paste0("Distances_in_",i),tempdf)
      halfsibs=c(halfsibs,tempdf[,6])
    }
  }
  
  write.table(halfsibs,paste0("histoHS_",samples.tag,"Distances_Final"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  xlimitsHS=c(min(halfsibs)-0.05,max(halfsibs)+0.05)
  hist(halfsibs,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsHS,ylab=paste0("N=",length(halfsibs)), main=paste0("Half Siblings for ",samples.tag))
}


################### RELATEDNESS ANALYSIS - PLINK2SPAGEDI ###########################
#----------------------------------------------------------------------------------#
# Using low-coverage, forced homozygous, SNP data in the PLINK format, this script #
# transforms it into a SPAGEDI-ready format for relatedness estimations of two     #
# samples.                                                                         #
# Allele frequencies should be provided, and can be retrieved by using the command #
# --freq with PLINK on a dataset of interest.                                      #
# 										                                                             #
# INPUT:                                                                           #
# foo.ped   -PLINK genotype file for sample1					                             #
# foo.map   -PLINK map file for sample1						                                 #
# bar.ped   -PLINK genotype file for sample2					                             #
# bar.map   -PLINK map file for sample2						                                 #
# plink.frq -allele frequencies in PLINK tableformat				                       #
#											                                                             #
# Note: SNP names in *.map files and in *.frq file must be exactly the same	       #


############# Relatedness test from forced homozygote SNP data ##############
# --------------------------------------------------------------------------#
relatedHomozSNP = function(sample1,sample2,freqs,run.spagedi=TRUE) {
  
  ################# Create tranposable read.csv function (from http://stackoverflow.com/questions/17288197/reading-a-csv-file-organized-horizontally)
  read.tcsv = function(file, header=TRUE, sep=",", ...) {
    n = max(count.fields(file, sep=sep), na.rm=TRUE)
    x = readLines(file)
    .splitvar = function(x, sep, n) {
      var = unlist(strsplit(x, split=sep))
      length(var) = n
      return(var)
    }
    x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
    x = apply(x, 1, paste, collapse=sep) 
    out = read.csv(text=x, sep=sep, header=header, ...)
    return(out)
  }
  
  ####### Sample1
  print(paste0("Reading sample 1: ",sample1))
  sample1_ped_temp=read.tcsv(paste0(sample1,".ped"),as.is=TRUE,header=FALSE,stringsAsFactors = FALSE,sep=" ")
  
  #######IMPORTANT - This next line is not always needed. If you get an error about end of line while reading tcsv, add a line in a text editor to the PED file and then use this to remove it.
  #sample1_ped_temp = head(sample1_ped_temp,-1)
  
  sample1_map=read.csv(paste0(sample1,".map"),sep="",stringsAsFactors = FALSE,header=FALSE)
  snps_samp1=length(sample1_map$V1)
  sample1_ped=data.frame(matrix(ncol=4,nrow=snps_samp1))
  sample1_ped$X1=sample1_map$V1
  sample1_ped$X2=sample1_map$V4
  i=7
  ter=1
  for(allele in sample1_ped_temp$V1) {
    while(i<=length(sample1_ped_temp$V1)) {
      biall=paste(sample1_ped_temp[i,1],sample1_ped_temp[i+1,1])
      sample1_ped$X3[ter]=sample1_ped_temp[i,1]
      sample1_ped$X4[ter]=sample1_ped_temp[i+1,1]
      ter=ter+1
      i=i+2
    }
  }
  sample1_ped$X5=paste(sample1_ped$X1,"-",sample1_ped$X2,sep="")
  sample1_ped$X6=sample1_map$V2
  ## Remove SNPs with no data (due to QC during PED creation)
  sample1_ped = sample1_ped[sample1_ped$X3 != 0,]
  
  
  ####### Sample 2
  print(paste0("Reading sample 2: ",sample2))
  sample2_ped_temp=read.tcsv(paste0(sample2,".ped"),as.is=TRUE,header=FALSE,stringsAsFactors = FALSE,sep=" ")
  
  #######IMPORTANT - This next line is not always needed. If you get an error about end of line while reading tcsv, add a line in a text editor to the PED file and then use this to remove it.
  #sample2_ped_temp = head(sample2_ped_temp,-1)
  
  sample2_map=read.csv(paste0(sample2,".map"),sep="",stringsAsFactors = FALSE,header=FALSE)
  snps_samp2=length(sample2_map$V1)
  sample2_ped=data.frame(matrix(ncol=4,nrow=snps_samp2))
  sample2_ped$X1=sample2_map$V1
  sample2_ped$X2=sample2_map$V4
  i=7
  ter=1
  for(allele in sample2_ped_temp$V1) {
    while(i<=length(sample2_ped_temp$V1)) {
      biall=paste(sample2_ped_temp[i,1],sample2_ped_temp[i+1,1])
      sample2_ped$X3[ter]=sample2_ped_temp[i,1]
      sample2_ped$X4[ter]=sample2_ped_temp[i+1,1]
      ter=ter+1
      i=i+2
    }
  }
  sample2_ped$X5=paste(sample2_ped$X1,"-",sample2_ped$X2,sep="")
  sample2_ped$X6=sample2_map$V2
  ## Remove SNPs with no data (due to QC during PED creation)
  sample2_ped = sample2_ped[sample2_ped$X3 != 0,]
  
  
  #################### Get common SNPs
  print("Finding common SNPs")
  all_SNPs=c(sample1_ped$X5,sample2_ped$X5)
  common_SNPs_logic=duplicated(all_SNPs)
  counter=1
  common_SNPs=data.frame(snp=as.character(),stringsAsFactors = FALSE)
  for(i in all_SNPs) {
    if(common_SNPs_logic[counter]==TRUE) {
      new_row=data.frame(snp=i)
      common_SNPs = rbind(common_SNPs,new_row)}
    counter=counter+1
  }
  
  ## New data.frame with data from common SNPs (most time-consuming step)
  commonDataFrame=data.frame(chr=as.numeric(),pos=as.numeric(),snp=as.character(),S1x=as.character(),S1y=as.character(),S2x=as.character(),S2y=as.character())
  for(i in common_SNPs$snp) {
    posIn1=which(i==sample1_ped$X5)
    posIn2=which(i==sample2_ped$X5)
    newRow=data.frame(sample1_ped[posIn1,1],sample1_ped[posIn1,6],sample1_ped[posIn1,2],sample1_ped[posIn1,3],sample1_ped[posIn1,4],sample2_ped[posIn2,3],sample2_ped[posIn2,4])
    commonDataFrame=rbind(commonDataFrame,newRow)
  }
  colnames(commonDataFrame) = c("chr","snp","pos","S1x","S1y","S2x","S2y")
  
  ################## Read allele frequencies
  print(paste0("Reading allele frequencies from file: ",freqs))
  alFreq = read.csv(freqs,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  #A1 is minor allele
  alFreq$NCHROBS=NULL
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ## Add frequencies to commonDataFrame
  ## WARNING: SNP names in frequency and map files have to be exactly the same
  counter=1
  for(i in commonDataFrame$snp) {
    posInAlFreq=which(i==alFreq$SNP)
    commonDataFrame$al1freq[counter]=alFreq[posInAlFreq,5]
    commonDataFrame$al2freq[counter]=alFreq[posInAlFreq,6]
    counter=counter+1
  }
  ## Remove SNPs with fixed alleles
  commonDataFrame = commonDataFrame[commonDataFrame$al1freq != 0,]
  ## Create alFreqCommon
  alFreqCommon=data.frame(colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  counter=1
  for(i in commonDataFrame$snp) {
    posInAlFreq=which(i==alFreq$SNP)
    newRow=alFreq[posInAlFreq,]
    alFreqCommon=rbind(alFreqCommon,newRow)
    counter=counter+1
  }
  
  ################# Swap alleles' letter by numericals
  commonDataFrameNumericals = commonDataFrame
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Test if both samples have all homozygous SNPs
  all(commonDataFrameNumericals$S1x == commonDataFrameNumericals$S1y)
  all(commonDataFrameNumericals$S2x == commonDataFrameNumericals$S2y)
  ## Do the same for alFreqCommon
  alFreqCommonNumericals = alFreqCommon
  alFreqCommonNumericals[] = lapply(alFreqCommonNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFreqCommonNumericals[] = lapply(alFreqCommonNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFreqCommonNumericals[] = lapply(alFreqCommonNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFreqCommonNumericals[] = lapply(alFreqCommonNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Work with 3 decimal places for frequencies and make sure MAF+AFa2=1
  maf=format(as.numeric(alFreqCommonNumericals$MAF),signif=3)
  alFreqCommonNumericals$MAF=maf
  af2=format((1-as.numeric(maf)),signif=3)
  alFreqCommonNumericals$AFa2=af2
  
  ################# Remove SNPs with zero frequency alleles
  print("Removing SNPs with zero frequency alleles")
  toRemove=c()
  ## Check in individual 1
  it=1
  for(i in commonDataFrameNumericals$S1x) {
    if((i!=alFreqCommonNumericals$A1[it])==TRUE && (i!=alFreqCommonNumericals$A2[it])==TRUE) {
      toRemove=c(toRemove,alFreqCommonNumericals$SNP[it])
      print(alFreqCommonNumericals[it,])}
    it=it+1}
  ## Check in individual 2
  it=1
  for(i in commonDataFrameNumericals$S2x) {
    if((i!=alFreqCommonNumericals$A1[it])==TRUE && (i!=alFreqCommonNumericals$A2[it])==TRUE) {
      toRemove=c(toRemove,alFreqCommonNumericals$SNP[it])
      print(alFreqCommonNumericals[it,])}
    it=it+1}
  ## Remove those SNPs
  for(i in toRemove) {
    pos=which(i==commonDataFrameNumericals$snp)
    commonDataFrameNumericals=commonDataFrameNumericals[-pos,]
    alFreqCommonNumericals=alFreqCommonNumericals[-pos,]
  }
  
  ################# Create SPAGEDI's data file
  print("Creating SPAGeDI's data file")
  maxCols=length(commonDataFrameNumericals$snp)
  row1=c(2,0,0,maxCols,3,2)
  row2="0"
  row3=c("IND",commonDataFrameNumericals$snp)
  row4=c(sample1,paste(commonDataFrameNumericals$S1x,commonDataFrameNumericals$S1y,sep=""))
  row5=c(sample2,paste(commonDataFrameNumericals$S2x,commonDataFrameNumericals$S2y,sep=""))
  row6="END"
  ## Export SPAGEDI data file
  sink(paste0(sample1,"_",sample2,"_spagedi_input.txt"))
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  cat("\n")
  cat(row4,sep="\t")
  cat("\n")
  cat(row5,sep="\t")
  cat("\n")
  cat(row6,sep="\t")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  print("Creating SPAGeDI's allele frequency  file")
  row1=c(paste(commonDataFrameNumericals$snp,"2",sep="\t"))
  row2=c(paste(alFreqCommonNumericals$A1,alFreqCommonNumericals$MAF,sep="\t"))
  row3=c(paste(alFreqCommonNumericals$A2,alFreqCommonNumericals$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0(sample1,"_",sample2,"_spagedi_frequencies.txt"))
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  print("Creating SPAGeDI's commands file")
  sink(paste0("SP_CMDS_Test",sample1,"_",sample2,"_in.txt"))
  cat(paste0(sample1,"_",sample2,"_spagedi_input.txt\n"))
  cat(paste0(sample1,"_",sample2,"_spagedi_out.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0(sample1,"_",sample2,"_spagedi_frequencies.txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    file.remove(paste0(sample1,"_",sample2,"_spagedi_out.txt\n"))
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_Test",sample1,"_",sample2,"_in.txt"))
    cmd = paste0("spagedi < SP_CMDS_Test",sample1,"_",sample2,"_in.txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
  
}
########### Transforms SPAGeDI Frequencies File into useable PLINK format for relatedness tests
freqCommonSNPs = function(spagedi.freq.file,dataset) {
  ids=strsplit(spagedi.freq.file,"_spagedi_frequencies.txt")
  freqsFile=read.table(paste0(spagedi.freq.file),as.is=TRUE,header=FALSE,stringsAsFactors = FALSE,sep="\t",nrows=1)
  to_rem=seq(from=2,to=length(freqsFile),by=2)
  freqsFile=freqsFile[,-to_rem]
  freqsFileT=t(freqsFile)
  write.table(freqsFileT,paste0("common_alleles_",ids),quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE)
  plinkCMD = paste0("plink --bfile ",dataset," --extract common_alleles_",ids," --make-bed --out common_",ids)
  system(plinkCMD,wait=TRUE)
  plinkCMD2 = paste0("plink --bfile common_",ids," --freq")
  system(plinkCMD2,wait=TRUE)
  file.rename(from="plink.frq", to=paste0("common_freqs_",ids,".frq"))
}













