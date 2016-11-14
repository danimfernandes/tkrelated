## TKrelated & CybRsex
## v.1.6
## 10 November 2016

## Maintained by: Daniel Fernandes - dani.mag.fernandes@gmail.com
## Written by Daniel Fernandes (main body, simulations) and John Finarelli (posterior probabilities)
## github.com/danimag/tkrelated
## Please cite as:
## "The Identification of a 1916 Irish Rebel: New Approach for Estimating Relatedness From Low Coverage Homozygous Genomes"
## 2016
## http://dx.doi.org/10.1101/076992

##CHANGELOG:
##    Added posterior probabilities
##    Multiple speed improvements on "relatedHomozSNP" function (up to 24X faster)
##    Fixed bug for "reduce.SNPs"=FALSE
##    Fixed bug due to presence of fixed alleles

##TODO:
##    Long-Term: Port to python


######################## TKrelated - RELATEDNESS ANALYSIS ###########################
#----------------------------------------------------------------------------------#
# Using low-coverage, forced homozygous, SNP data in the PLINK format, this        #
# function transforms it into a SPAGEDI-ready format for relatedness estimations   #
# of two samples.                                                                  #
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

######################## Relatedness test from forced homozygote SNP data ##############
# --------------------------------------------------------------------------#
relatedHomozSNP = function(sample1,sample2,freqs,run.spagedi=TRUE) {

  ## Arguments' description:
  ##
  ## sample1 - name/location of text plink sample 1
  ## sample2 - name/location of text plink sample 2
  ## freqs - name/location of the plink frequencies file from desired population
  ## run.spagedi - boolean for whether to automatically run SPAGeDI or not
  
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
  ## Get allele information from PED file
  sample1_ped_temp = data.frame(X1=sample1_ped_temp[7:length(sample1_ped_temp$V1),],stringsAsFactors = F)
  alA = seq(1,length(sample1_ped_temp[,1]),by=2)
  alB = seq(2,length(sample1_ped_temp[,1]),by=2)
  sample1_ped$X3 = sample1_ped_temp[alA,]
  sample1_ped$X4 = sample1_ped_temp[alB,]
  sample1_ped$X5=paste(sample1_ped$X1,"-",sample1_ped$X2,sep="")
  sample1_ped$X6=sample1_map$V2
  ## Remove SNPs with no data (due to QC during PED creation)
  sample1_ped = sample1_ped[sample1_ped$X3 != 0,]
  ## Remove SNPs with different names but same coordinates. Randomly selects all but one to remain
  iprev = "NA"
  to_rem_ped1=c()
  posia=c()
  for(i in sample1_ped$X5) {
    if(i == iprev) {
      posi=which(i == sample1_ped$X5)
      posia=sample(posi,length(posi)-1)
      to_rem_ped1=c(to_rem_ped1,posia)
    }
    iprev=i
  }
  try(expr = (sample1_ped = sample1_ped[-to_rem_ped1,]) ,silent = TRUE)
  
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
  ## Get allele information from PED file
  sample2_ped_temp = data.frame(X1=sample2_ped_temp[7:length(sample2_ped_temp$V1),],stringsAsFactors = F)
  alA2 = seq(1,length(sample2_ped_temp[,1]),by=2)
  alB2 = seq(2,length(sample2_ped_temp[,1]),by=2)
  sample2_ped$X3 = sample2_ped_temp[alA2,]
  sample2_ped$X4 = sample2_ped_temp[alB2,]
  sample2_ped$X5=paste(sample2_ped$X1,"-",sample2_ped$X2,sep="")
  sample2_ped$X6=sample2_map$V2
  ## Remove SNPs with no data (due to QC during PED creation)
  sample2_ped = sample2_ped[sample2_ped$X3 != 0,]
  ## Remove SNPs with different names but same coordinates. Randomly selects all but one to remain
  iprev = "NA"
  to_rem_ped2=c()
  posia=c()
  for(i in sample2_ped$X5) {
    if(i == iprev) {
      posi=which(i == sample2_ped$X5)
      posia=sample(posi,length(posi)-1)
      to_rem_ped2=c(to_rem_ped2,posia)
    }
    iprev=i
  }
  try(expr = (sample2_ped = sample2_ped[-to_rem_ped2,]) ,silent = TRUE)

  ##### Get common SNPs
  print("Finding common SNPs")
  all_SNPs = c(sample1_ped$X5,sample2_ped$X5)
  common_SNPs = data.frame(snp = all_SNPs[which(duplicated(all_SNPs) == TRUE)], stringsAsFactors = F)
  
  ## New data.frame with data from common SNPs
  posIn1=which(sample1_ped$X5 %in% common_SNPs$snp)
  posIn2=which(sample2_ped$X5 %in% common_SNPs$snp)
  commonDataFrame = cbind(sample1_ped[posIn1,1],sample1_ped[posIn1,6],sample1_ped[posIn1,2],sample1_ped[posIn1,3],sample1_ped[posIn1,4],sample2_ped[posIn2,3],sample2_ped[posIn2,4])
  commonDataFrame = data.frame(commonDataFrame)
  colnames(commonDataFrame) = c("chr","snp","pos","S1x","S1y","S2x","S2y")
  
  ################## Read allele frequencies
  print(paste0("Reading allele frequencies from file: ",freqs))
  alFreq = read.csv(freqs,stringsAsFactors=FALSE,sep="",colClasses = c(character(),character(),character(),character(),numeric(),numeric()))
  #A1 is minor allele
  alFreq$NCHROBS=NULL
  alFreq$AFa2=format((1-as.numeric(alFreq$MAF)),digits=3)
  alFreq[5]=lapply(alFreq[5],round,3)
  ## Add frequencies to commonDataFrame
  ## WARNING: SNP names in frequency and map files have to be an exact match
  posInAlFreq = which(alFreq$SNP %in% commonDataFrame$snp)
  commonDataFrame$al1freq=alFreq[posInAlFreq,5]
  commonDataFrame$al2freq=alFreq[posInAlFreq,6]
  ## Remove SNPs with fixed alleles
  commonDataFrame = commonDataFrame[commonDataFrame$al1freq != 0,]
  ## Create alFreqCommon
  posInAlFreq = which(alFreq$SNP %in% commonDataFrame$snp)
  alFreqCommon = alFreq[posInAlFreq,]

  ## Export tables with common allele IDs and frequencies
  write.table(commonDataFrame$snp,paste0("common_alleles_",sample1,"_",sample2),quote = FALSE,sep="\t", row.names = FALSE, col.names = FALSE)
  alFreqCommonToExport = alFreqCommon
  names(alFreqCommonToExport)[6] = "NCHROBS"
  alFreqCommonToExport$NCHROBS = NA
  write.table(alFreqCommonToExport,paste0("commons_freqs_",sample1,"_",sample2,".frq"), quote=FALSE,row.names = FALSE,col.names = TRUE,sep="\t")
  
  ################# Swap alleles' letter by numericals
  commonDataFrameNumericals = commonDataFrame
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  commonDataFrameNumericals[] = lapply(commonDataFrameNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Test if both samples have all homozygous SNPs, and force so if FALSE
  if(all(commonDataFrameNumericals$S1x == commonDataFrameNumericals$S1y) == FALSE) {
    quais1 = which(commonDataFrameNumericals$S1x != commonDataFrameNumericals$S1y)
    for(i in quais1) {
      chosenAl=sample(c(commonDataFrameNumericals[i,4],commonDataFrameNumericals[i,5]))[1]
      commonDataFrameNumericals[i,4] = chosenAl
      commonDataFrameNumericals[i,5] = chosenAl
    }    
  }
  if(all(commonDataFrameNumericals$S2x == commonDataFrameNumericals$S2y) == FALSE) {
    quais2 = which(commonDataFrameNumericals$S2x != commonDataFrameNumericals$S2y)
    for(i in quais2) {
      chosenAl=sample(c(commonDataFrameNumericals[i,6],commonDataFrameNumericals[i,7]))[1]
      commonDataFrameNumericals[i,6] = chosenAl
      commonDataFrameNumericals[i,7] = chosenAl
    }    
  }
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
      }
    it=it+1}
  ## Check in individual 2
  it=1
  for(i in commonDataFrameNumericals$S2x) {
    if((i!=alFreqCommonNumericals$A1[it])==TRUE && (i!=alFreqCommonNumericals$A2[it])==TRUE) {
      toRemove=c(toRemove,alFreqCommonNumericals$SNP[it])
      }
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
    try(ile.remove(paste0(sample1,"_",sample2,"_spagedi_out.txt\n")),silent=TRUE)
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_Test",sample1,"_",sample2,"_in.txt"))
    cmd = paste0("spagedi < SP_CMDS_Test",sample1,"_",sample2,"_in.txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}

######################## CybRsex- SIMULATING RELATEDNESS #############################
#---------------------------------------------------------------------------------#
# This script generates random sets of related and unrelated individuals based on #
# SNP allele frequencies. These must be provided in PLINK format.		              #
# Using the third-party software SPAGeDI, relatedness coefficients for these      #
# individuals are calculated and then plotted as histograms.			                #
#											                                                            #
# INPUT:										                                                      #
# plink.frq -allele frequencies in PLINK tableformat	                  			    #

######################## SET OF UNRELATED INDIVIDUALS (TRC 0%) #######################
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
  print("Generating virtual individuals based on allele frequencies")
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

######################## SET OF FULL SIBLINGS (TRC 50%) ###########################
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
  }  else {   alFrequ=alFreq  }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  print("Generating virtual individuals based on allele frequencies")
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
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
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
  hist(fullsibs,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsFS,ylab=paste0("N=",length(fullsibs)), main=paste0("First Order for ",samples.tag))
}

######################## SET FOR HALF SIBLINGS (TRC 25%) #######################
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
  }  else {    alFrequ=alFreq   }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  print("Generating virtual individuals based on allele frequencies")
  numIndividuals = numHalfSiblingsPairs*2 
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
      halfSib[paste("HSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextNextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      halfSib[paste("HSIB",iter,"-",nextNextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextNextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextNextInd,"b",sep="")])[[1]])
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
  hist(halfsibs,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsHS,ylab=paste0("N=",length(halfsibs)), main=paste0("Second Order for ",samples.tag))
}

######################## SET FOR 3/4 SIBLINGS (TRC 37.5%) #######################
# ------------------------------------------------------------------------- #
threequartersibsfunc = function(file,samples.tag,num34SiblingsPairs,identifier,run.spagedi=TRUE,reduce.SNPs=FALSE) {

  ## Arguments' description:
  ##
  ## file - name/location of the plink frequencies file to simulate individuals from
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## num34siblingPairs - the number of 3/4-sibling pairs of individuals to generate
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
  }  else {    alFrequ=alFreq   }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  print("Generating virtual individuals based on allele frequencies")
  ## UNRELATED Individuals
  numIndividuals = num34SiblingsPairs*3  ## Multiplier depends on number of unique individuals needed per final pair
  seqIndividuals = seq(1:numIndividuals)
  alFrequ=alFrequNumericals
  for(i in seqIndividuals) {
    alFrequ[paste("IND",i,"a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
    alFrequ[paste("IND",i,"b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
  }
  
  ################# Choosing one allele at random from parent to create full siblings
  fullSib=alFrequ[1:5]
  iter=1
  for(i in seqIndividuals){
    nextInd=iter+1
    if(nextInd <= numIndividuals) {
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
    }
    iter=iter+3
  }
  fullSib = fullSib[,-(1:5)]
  
  ################## Generating three-quarter siblings from full siblings and unrelated individuals
  threequarterSib=alFrequ[1:5]
  iter=1
  for(ind in seq(1:(num34SiblingsPairs*2))) {
    nextInd = iter+1
    samePar = iter+2
    if(nextInd <= numIndividuals) {
      threequarterSib[paste("TQSIB",iter,"-",nextInd,"I-",samePar,"_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(fullSib[paste("FSIB",iter,"-",nextInd,"I_a",sep="")])[[1]], c(fullSib[paste("FSIB",iter,"-",nextInd,"I_b",sep="")])[[1]])
      threequarterSib[paste("TQSIB",iter,"-",nextInd,"I-",samePar,"_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",samePar,"a",sep="")])[[1]], c(alFrequ[paste("IND",samePar,"b",sep="")])[[1]])
      threequarterSib[paste("TQSIB",iter,"-",nextInd,"II-",samePar,"_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(fullSib[paste("FSIB",iter,"-",nextInd,"II_a",sep="")])[[1]], c(fullSib[paste("FSIB",iter,"-",nextInd,"II_b",sep="")])[[1]])
      threequarterSib[paste("TQSIB",iter,"-",nextInd,"II-",samePar,"_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",samePar,"a",sep="")])[[1]], c(alFrequ[paste("IND",samePar,"b",sep="")])[[1]])
    }
    iter = iter+3
  }
  threequarterSib = threequarterSib[,-(1:5)]
  
  ################### Create a list of the pairs of three-quarter-siblings
  list34Sibs = list()
  numSiblingPairs = num34SiblingsPairs
  seqSiblings = seq(1:numSiblingPairs)
  iter=1
  for(i in seqSiblings) {
    nextInd=iter+1
    samePar=iter+2
    ind = paste0("TQSIB",iter,"-",nextInd,"I-",samePar) 
    ind2 = paste0("TQSIB",iter,"-",nextInd,"II-",samePar) 
    list34Sibs = c(list34Sibs, ind,ind2)
    iter=iter+3
  }
  
  ################# Force homozygous
  threequarterSibFH=alFrequ[1:5]
  for(i in list34Sibs){
    try((threequarterSibFH[paste0(i)] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(threequarterSib[paste(i,"_a",sep="")])[[1]], c(threequarterSib[paste(i,"_b",sep="")])[[1]])),silent=TRUE)
  }
  threequarterSibFH = threequarterSibFH[,-(1:5)]
  
  ################# Create and export SPAGEDI's data file
  ## Remove previous SPAGeDI output files 
  for(iterit in list.files(getwd())) {
    if(grepl(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt"),iterit)== TRUE) {
      system(paste0("rm -f ",iterit))
    }
  }
  ## Create input file
  maxCols=length(alFrequ$SNP)
  sink(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  row1=c(length(list34Sibs),0,0,maxCols,3,2)
  cat(row1,sep="\t")
  cat("\n")
  row2="0"
  cat(row2,sep="\t")
  cat("\n")
  row3=c("IND",alFrequ$SNP)
  cat(row3,sep="\t")
  cat("\n")
  iter=4  ## First row in SPAGEDI input file for samples
  it=1
  for(i in list34Sibs){
    colNum=grep(paste0("^",i,"$"),c(names(threequarterSibFH)))
    cat(c(i,paste(threequarterSibFH[,colNum],threequarterSibFH[,colNum],sep="")),sep="\t")
    cat("\n")
    it=it+1
  }
  cat("END")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  row1=c(paste(alFrequ$SNP,"2",sep="\t"))
  row2=c(paste(alFrequ$A1,alFrequ$MAF,sep="\t"))
  row3=c(paste(alFrequ$A2,alFrequ$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"),append=FALSE)
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  sink(paste0("SP_CMDS_ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,".txt"))
  cat(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,".txt\n"))
  cat(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,"_spgout.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0("ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_freq",identifier,".txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,".txt" ))
    cmd = paste0("spagedi < SP_CMDS_ThreeQuarterSibs",num34SiblingsPairs,"_",samples.tag,"_in",identifier,".txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}
########### Function to read output from SPAGeDI and plot it according to provided samples.tag
threequartersibsplot = function(samples.tag) {
  
  wdfiles=list.files(getwd())
  it=1
  for(filet in wdfiles) {
    if(grepl("_spgout.txt$",filet) == TRUE) {
      if(grepl("^ThreeQuarter",filet) == TRUE) {
        if(grepl(paste0("_",samples.tag,"_"),filet) == TRUE) {
          print(paste0("Reading SPAGeDi output file: ",filet))
          assign(paste0("histoTQ","_",samples.tag,"_",it),read.csv(filet,header = FALSE))
          it=it+1
        }
      }
    }
  }
  
  ## Read until: "PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form"
  threequartersibs = c()
  enviro=ls()
  for(i in enviro) {
    if(grepl(paste0("^histoTQ","_",samples.tag,"_"),i)==TRUE){
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
      ## Remove loop
      it=1
      to_keep = c()
      for(x in seq(1:(numInd/2))) {
        #it=it+2
        nextIt=it+1
        keepPos = which(tempdf[,3]==it & tempdf[,4]==nextIt)
        to_keep = c(to_keep,keepPos)
        it = it+2
      }
      posis=seq(1:length(tempdf[,1]))
      to_remove = setdiff(posis, as.numeric(to_keep))
      tempdf=tempdf[-to_remove,]
      assign(paste0("Distances_in_",i),tempdf)
      threequartersibs=c(threequartersibs,tempdf[,6])
    }
  }
  
  write.table(threequartersibs,paste0("histoTQ_",samples.tag,"Distances_Final"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  xlimitsFS=c(min(threequartersibs)-0.05,max(threequartersibs)+0.05)
  hist(threequartersibs,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsFS,ylab=paste0("N=",length(threequartersibs)), main=paste0("First-&-Half Order for ",samples.tag))
}

######################## SET FOR FIRST COUSINS (TRC 12.5%) #####################
# ------------------------------------------------------------------------- #
firstcousfunc = function(file,samples.tag,numfirstcousPairs,identifier,run.spagedi=TRUE,reduce.SNPs=FALSE) {
  
  ## Arguments' description:
  ##
  ## file - name/location of the plink frequencies file to simulate individuals from
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## numfirstcousPairs - the number of 3/4-sibling pairs of individuals to generate
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
  }  else {    alFrequ=alFreq   }
  
  ################# Swap alleles' letter by numericals
  alFrequNumericals = alFrequ
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "A$", replacement = 100)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "C$", replacement = 200)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "G$", replacement = 300)#, fixed = TRUE)
  alFrequNumericals[] = lapply(alFrequNumericals, gsub, pattern = "T$", replacement = 400)#, fixed = TRUE)
  ## Rename the data.frame back to original
  alFrequ=alFrequNumericals
  
  ################# Generate individual+allele columns PLUS choose allele according to RANDOM value
  print("Generating virtual individuals based on allele frequencies")
  ## UNRELATED Individuals
  numIndividuals = numfirstcousPairs*4  ## Multiplier depends on number of unique individuals needed per final pair
  seqIndividuals = seq(1:numIndividuals)
  alFrequ=alFrequNumericals
  for(i in seqIndividuals) {
    alFrequ[paste("IND",i,"a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
    alFrequ[paste("IND",i,"b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < alFrequ$MAF, alFrequ$A1, alFrequ$A2)
  }
  
  ################# Choosing one allele at random from parent to create full siblings
  fullSib=alFrequ[1:5]
  iter=1
  for(i in seqIndividuals){
    nextInd=iter+1
    if(nextInd <= numIndividuals) {
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"I_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",iter,"a",sep="")])[[1]], c(alFrequ[paste("IND",iter,"b",sep="")])[[1]])
      fullSib[paste("FSIB",iter,"-",nextInd,"II_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",nextInd,"a",sep="")])[[1]], c(alFrequ[paste("IND",nextInd,"b",sep="")])[[1]])
    }
    iter=iter+4
  }
  fullSib = fullSib[,-(1:5)]
  
  ################## Generating first cousins from full siblings and unrelated individuals
  firstCous=alFrequ[1:5]
  iter=1
  for(ind in seq(1:(numfirstcousPairs*2))) {
    nextInd = iter+1
    parOne = iter+2
    parTwo = iter+3
    if(nextInd <= numIndividuals) {
      firstCous[paste("FCOU",iter,"-",nextInd,"I-",parOne,"_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(fullSib[paste("FSIB",iter,"-",nextInd,"I_a",sep="")])[[1]], c(fullSib[paste("FSIB",iter,"-",nextInd,"I_b",sep="")])[[1]])
      firstCous[paste("FCOU",iter,"-",nextInd,"I-",parOne,"_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",parOne,"a",sep="")])[[1]], c(alFrequ[paste("IND",parOne,"b",sep="")])[[1]])
      firstCous[paste("FCOU",iter,"-",nextInd,"II-",parTwo,"_","a",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(fullSib[paste("FSIB",iter,"-",nextInd,"II_a",sep="")])[[1]], c(fullSib[paste("FSIB",iter,"-",nextInd,"II_b",sep="")])[[1]])
      firstCous[paste("FCOU",iter,"-",nextInd,"II-",parTwo,"_","b",sep="")] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(alFrequ[paste("IND",parTwo,"a",sep="")])[[1]], c(alFrequ[paste("IND",parTwo,"b",sep="")])[[1]])
    }
    iter = iter+4
  }
  firstCous = firstCous[,-(1:5)]
  
  ################### Create a list of the pairs of first cousins
  listFirstCous = list()
  numSiblingPairs = numfirstcousPairs
  seqSiblings = seq(1:numSiblingPairs)
  iter=1
  for(i in seqSiblings) {
    nextInd=iter+1
    parOne=iter+2
    parTwo=iter+3
    ind = paste0("FCOU",iter,"-",nextInd,"I-",parOne) 
    ind2 = paste0("FCOU",iter,"-",nextInd,"II-",parTwo) 
    listFirstCous = c(listFirstCous, ind,ind2)
    iter=iter+4
  }
  
  ################# Force homozygous
  firstCousFH=alFrequ[1:5]
  for(i in listFirstCous){
    try((firstCousFH[paste0(i)] = ifelse(runif(length(alFrequ$SNP),0.00,1.00) < 0.5, c(firstCous[paste(i,"_a",sep="")])[[1]], c(firstCous[paste(i,"_b",sep="")])[[1]])),silent=TRUE)
  }
  firstCousFH = firstCousFH[,-(1:5)]
  
  ################# Create and export SPAGEDI's data file
  ## Remove previous SPAGeDI output files 
  for(iterit in list.files(getwd())) {
    if(grepl(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,"_spgout.txt"),iterit)== TRUE) {
      system(paste0("rm -f ",iterit))
    }
  }
  ## Create input file
  maxCols=length(alFrequ$SNP)
  sink(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,".txt"))
  row1=c(length(listFirstCous),0,0,maxCols,3,2)
  cat(row1,sep="\t")
  cat("\n")
  row2="0"
  cat(row2,sep="\t")
  cat("\n")
  row3=c("IND",alFrequ$SNP)
  cat(row3,sep="\t")
  cat("\n")
  iter=4  ## First row in SPAGEDI input file for samples
  it=1
  for(i in listFirstCous){
    colNum=grep(paste0("^",i,"$"),c(names(firstCousFH)))
    cat(c(i,paste(firstCousFH[,colNum],firstCousFH[,colNum],sep="")),sep="\t")
    cat("\n")
    it=it+1
  }
  cat("END")
  sink()
  
  ################# Create SPAGEDI's allele frequency file
  row1=c(paste(alFrequ$SNP,"2",sep="\t"))
  row2=c(paste(alFrequ$A1,alFrequ$MAF,sep="\t"))
  row3=c(paste(alFrequ$A2,alFrequ$AFa2,sep="\t"))
  ## Export SPAGEDI frequencies file
  sink(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_freq",identifier,".txt"),append=FALSE)
  cat(row1,sep="\t")
  cat("\n")
  cat(row2,sep="\t")
  cat("\n")
  cat(row3,sep="\t")
  sink()
  
  ################ Create SPAGEDI's commands file for our type of analysis
  sink(paste0("SP_CMDS_FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,".txt"))
  cat(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,".txt\n"))
  cat(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,"_spgout.txt\n"))
  cat("\n4\n6\n\n")
  cat(paste0("FirstCousins",numfirstcousPairs,"_",samples.tag,"_freq",identifier,".txt"))
  cat("\n4\n3")
  sink()
  
  if(run.spagedi == TRUE) {
    print(paste0("Running SPAGeDI with commands from: ","SP_CMDS_FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,".txt" ))
    cmd = paste0("spagedi < SP_CMDS_FirstCousins",numfirstcousPairs,"_",samples.tag,"_in",identifier,".txt")
    system(cmd, ignore.stdout = TRUE, wait=TRUE)
  }
}
########### Function to read output from SPAGeDI and plot it according to provided samples.tag
firstcousplot = function(samples.tag) {
  
  wdfiles=list.files(getwd())
  it=1
  for(filet in wdfiles) {
    if(grepl("_spgout.txt$",filet) == TRUE) {
      if(grepl("^FirstCousins",filet) == TRUE) {
        if(grepl(paste0("_",samples.tag,"_"),filet) == TRUE) {
          print(paste0("Reading SPAGeDi output file: ",filet))
          assign(paste0("histoFC","_",samples.tag,"_",it),read.csv(filet,header = FALSE))
          it=it+1
        }
      }
    }
  }
  
  ## Read until: "PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form"
  firstcous = c()
  enviro=ls()
  for(i in enviro) {
    if(grepl(paste0("^histoFC","_",samples.tag,"_"),i)==TRUE){
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
      ## Remove loop
      it=1
      to_keep = c()
      for(x in seq(1:(numInd/2))) {
        #it=it+2
        nextIt=it+1
        keepPos = which(tempdf[,3]==it & tempdf[,4]==nextIt)
        to_keep = c(to_keep,keepPos)
        it = it+2
      }
      posis=seq(1:length(tempdf[,1]))
      to_remove = setdiff(posis, as.numeric(to_keep))
      tempdf=tempdf[-to_remove,]
      assign(paste0("Distances_in_",i),tempdf)
      firstcous=c(firstcous,tempdf[,6])
    }
  }
  
  write.table(firstcous,paste0("histoFC_",samples.tag,"Distances_Final"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  xlimitsFS=c(min(firstcous)-0.05,max(firstcous)+0.05)
  hist(firstcous,col = "orangered",breaks=25,xlab="Relatedness Coefficient",xlim=xlimitsFS,ylab=paste0("N=",length(firstcous)), main=paste0("Third Order for ",samples.tag))
}

######################## REAL ALL ESTIMATIONS ###################
########### Function to read all outputs from SPAGeDI and plot them according to provided samples.tag
readAll = function(samples.tag,plot.out=TRUE,stats.out=TRUE,pdf.w=25,pdf.h=15,pdf.cex=30) {

  ## Arguments' description:
  ##
  ## samples.tag - a user-defined tag to identify the simulation, for example, a tag identifying the pair of individuals the frequencies were extracted from (BobJane)
  ## plot.out - boolean for whether to create an output PDF of the plot or not
  ## stats.out - boolean for whether to create an output TXT tab-spaced file with posterior probabilities
  ## pdf.w - width of the PDF file in inches
  ## pdf.h - height of the PDF file in inches
  ## pdf.cex - point size for when exporting to PDF
  
  library(plotrix)
  
  ## Read in the five relatedness simulations
  histoUN = read.table(paste0("histoUN_",samples.tag,"Distances_Final"))  
  histoHS = read.table(paste0("histoHS_",samples.tag,"Distances_Final"))  
  histoFS = read.table(paste0("histoFS_",samples.tag,"Distances_Final"))
  histoTQ = read.table(paste0("histoTQ_",samples.tag,"Distances_Final"))
  histoFC = read.table(paste0("histoFC_",samples.tag,"Distances_Final"))
  
  n_size = length(histoUN$V1)+length(histoHS$V1)+length(histoFS$V1)+length(histoFC$V1)+length(histoTQ$V1)
  x_lims = c(min(c(histoHS$V1,histoFS$V1,histoUN$V1)-0.05),max(c(histoHS$V1,histoFS$V1,histoUN$V1)+0.05))
  
  ## Read in the relatedness coefficient for this pair of samples and number of SNPs used
  fileIn = read.csv(paste0(samples.tag,"_spagedi_out.txt"),stringsAsFactors = F)
  coeff = which("PAIRWISE SPATIAL AND GENETIC DISTANCES written in column form" == fileIn)+3
  coefficientRelatedness = as.numeric(tail(strsplit(as.character(fileIn[coeff,]),"\t")[[1]], n=1))
  snps = grep("loci:",fileIn[,1])
  commsnps = as.numeric(strsplit(fileIn[snps,1],"loci:")[[1]][1])
  
  ## Plot multiple histogram
  if(plot.out==TRUE) {
    pdf(file=paste0(getwd(),"/",samples.tag,"readAll.pdf"),width=pdf.w,height=pdf.h, pointsize=pdf.cex,onefile = FALSE) 
  }
  par(mar=c(4.5,4.5,3,2),mgp=c(3.2,2,1))
  breaksUN = seq(min(histoUN$V1),max(histoUN$V1),0.01)
  hist(histoUN$V1,xlim=c(-0.15,0.40),col=rgb(42,77,110,190, maxColorValue = 255), breaks=length(breaksUN), xlab="Halved Relatedness Coefficient", ylab=paste0("N=",as.numeric(n_size)), main=paste0("Relatedness Coefficients for ",samples.tag),xaxt="n")
  hist(histoUN$V1,xlim=c(-0.15,0.40),col=rgb(42,77,110,190, maxColorValue = 255), breaks=length(breaksUN), xlab="Halved Relatedness Coefficient", ylab=paste0("N=",as.numeric(n_size)), main=paste0("Relatedness Coefficients for ",samples.tag),xaxt="n", ylim=c(0,round(par("usr")[4]+(par("usr")[4]*0.1),-1)))
  breaksFS = seq(min(histoFS$V1),max(histoFS$V1),0.01)
  hist(histoFS$V1,add=TRUE,col=rgb(168,56,59,190, maxColorValue = 255),breaks=length(breaksFS))
  breaksHS = seq(min(histoHS$V1),max(histoHS$V1),0.01)
  hist(histoHS$V1,add=TRUE,col=rgb(154,166,55,190,maxColorValue = 255),breaks=length(breaksHS))
  breaksTQ = seq(min(histoTQ$V1),max(histoTQ$V1),0.01)
  hist(histoTQ$V1,add=TRUE,col=rgb(240,175,19,190,maxColorValue = 255),breaks=length(breaksTQ))
  breaksFC = seq(min(histoFC$V1),max(histoFC$V1),0.01)
  hist(histoFC$V1,add=TRUE,col=rgb(53,158,144,190,maxColorValue = 255),breaks=length(breaksFC))
  axis(side=1, at=c(-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4))
  legend("topright",bty="n",legend=c("Unrelated","Third Order","Second Order","3/4 Siblings","First Order"),pch=21,pt.bg=c(rgb(42,77,110,190, maxColorValue = 255),rgb(53,158,144,190,maxColorValue = 255),rgb(154,166,55,190,maxColorValue = 255),rgb(240,175,19,190,maxColorValue = 255),rgb(168,56,59,190, maxColorValue = 255)),col="black",y.intersp=1,cex=1,pt.cex=2)

  ## Draw line and text for relatedness coefficient
  segments(x0=coefficientRelatedness, y0=-0.5, y1= par("usr")[4]/20, col="black",lwd=12)
  segments(x0=coefficientRelatedness, y0=-0.5, y1= par("usr")[4]/20, col="gold",lwd=10)
  text(coefficientRelatedness, -par("usr")[4]/50, paste0("r=",format(coefficientRelatedness,digits = 3)))
  
  ## Calculate posterior probabilities
  normal.likelihood <- function(params,x){
    mu <- params[1]
    var <- params[2]
    logl <- -0.5*log(var) - (1/(2*var)) * ((x-mu)**2)
    return(logl)
  }
  paramsUN <- c((sum(as.numeric(histoUN$V1))/length(as.numeric(histoUN$V1))),((length(as.numeric(histoUN$V1))-1)/length(as.numeric(histoUN$V1)))*(sd(as.numeric(histoUN$V1))**2))
  paramsFC <- c((sum(as.numeric(histoFC$V1))/length(as.numeric(histoFC$V1))),((length(as.numeric(histoFC$V1))-1)/length(as.numeric(histoFC$V1)))*(sd(as.numeric(histoFC$V1))**2))
  paramsHS <- c((sum(as.numeric(histoHS$V1))/length(as.numeric(histoHS$V1))),((length(as.numeric(histoHS$V1))-1)/length(as.numeric(histoHS$V1)))*(sd(as.numeric(histoHS$V1))**2))
  paramsTQ <- c((sum(as.numeric(histoTQ$V1))/length(as.numeric(histoTQ$V1))),((length(as.numeric(histoTQ$V1))-1)/length(as.numeric(histoTQ$V1)))*(sd(as.numeric(histoTQ$V1))**2))
  paramsFS <- c((sum(as.numeric(histoFS$V1))/length(as.numeric(histoFS$V1))),((length(as.numeric(histoFS$V1))-1)/length(as.numeric(histoFS$V1)))*(sd(as.numeric(histoFS$V1))**2))
  UNnl = normal.likelihood(paramsUN,coefficientRelatedness)
  FCnl = normal.likelihood(paramsFC,coefficientRelatedness)
  HSnl = normal.likelihood(paramsHS,coefficientRelatedness)
  TQnl = normal.likelihood(paramsTQ,coefficientRelatedness)
  FSnl = normal.likelihood(paramsFS,coefficientRelatedness)

  log.like = c(UNnl,FCnl,HSnl,TQnl,FSnl)
  likelihoods = exp(log.like)
  post.probs = sprintf("%.3f",(as.numeric(likelihoods/(sum(likelihoods)))),3)
  
  ## Add posterior probabilities to plot
  text(median(histoUN$V1),(round(par("usr")[4])/4),labels = post.probs[1],col="white")
  text(median(histoFC$V1),(round(par("usr")[4])/4),labels = post.probs[2],col="white")
  text(median(histoHS$V1),(round(par("usr")[4])/4),labels = post.probs[3],col="white")
  text(median(histoTQ$V1),(round(par("usr")[4])/4),labels = post.probs[4],col="white")
  text(median(histoFS$V1),(round(par("usr")[4])/4),labels = post.probs[5],col="white")
  
  ## Export probabilities to text file
  if(stats.out==TRUE) {
    tblUp = data.frame(c(samples.tag,"SNPs","Homozygous relatedness coefficient",""),c("",commsnps,coefficientRelatedness,""))
    postProbTbl = data.frame(c("Relationship","Unrelated","Third Order","Second Order","3/4 Siblings","First Order"),c("Posterior Probability",post.probs))
    colnames(tblUp) = c("A1","B1")
    colnames(postProbTbl) = c("A1","B1")
    merged = rbind(tblUp,postProbTbl)
    write.table(merged,file=paste0("PosteriorProbs_",samples.tag),quote=FALSE,sep="\t",row.names=F,col.names = F)
  }
  
  if(plot.out==TRUE) {
    dev.off()
  }
  
}

######################## EXAMPLE OF USAGE #####################

## Example of usage:
## 1000 simulations for each order 
## (divided by 5 iterations for time efficiency)

# ## Set working directory
# setwd("foo/bar")
# ## Run relatedness test for "ind1" and "ind2"
# relatedHomozSNP("ind1","ind2","my_freqs.frq",run.spagedi = TRUE)
# ## Generate unrelated individuals and simulate relatedness
# unrelatedfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,1,run.spagedi = TRUE)
# unrelatedfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,2,run.spagedi = TRUE)
# unrelatedfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,3,run.spagedi = TRUE)
# unrelatedfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,4,run.spagedi = TRUE)
# unrelatedfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,5,run.spagedi = TRUE)
# ## Use the plot function to compile the 5 individual runs and export it to be used on readAll
# unrelatedplot("ind1_ind2")
# ## ### Generate first order relatives and simulate relatedness
# fullsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,1,run.spagedi = TRUE)
# fullsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,2,run.spagedi = TRUE)
# fullsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,3,run.spagedi = TRUE)
# fullsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,4,run.spagedi = TRUE)
# fullsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,5,run.spagedi = TRUE)
# ## Use the plot function to compile the 5 individual runs and export it to be used on readAll
# fullsibsplot("ind1_ind2")
# ## ### ### Generate second order relatives and simulate relatedness
# halfsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,1,run.spagedi = TRUE)
# halfsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,2,run.spagedi = TRUE)
# halfsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,3,run.spagedi = TRUE)
# halfsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,4,run.spagedi = TRUE)
# halfsibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,5,run.spagedi = TRUE)
# ## Use the plot function to compile the 5 individual runs and export it to be used on readAll
# halfsibsplot("ind1_ind2")
# ## ### ### Generate 1.5 order relatives and simulate relatedness
# threequartersibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,1,run.spagedi = TRUE)
# threequartersibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,2,run.spagedi = TRUE)
# threequartersibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,3,run.spagedi = TRUE)
# threequartersibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,4,run.spagedi = TRUE)
# threequartersibsfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,5,run.spagedi = TRUE)
# ## Use the plot function to compile the 5 individual runs and export it to be used on readAll
# threequartersibsplot("ind1_ind2")
# ## ### ### Generate third order relatives and simulate relatedness
# firstcousfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,1,run.spagedi = TRUE)
# firstcousfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,2,run.spagedi = TRUE)
# firstcousfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,3,run.spagedi = TRUE)
# firstcousfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,4,run.spagedi = TRUE)
# firstcousfunc("commons_freqs_ind1_ind2.frq","ind1_ind2",200,5,run.spagedi = TRUE)
# ## Use the plot function to compile the 5 individual runs and export it to be used on readAll
# firstcousplot("ind1_ind2")
# ## Collect data from the 5 different orders and plot them together
# readAll("ind1_ind2")
