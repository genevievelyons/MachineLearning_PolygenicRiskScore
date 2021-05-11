library(readr)
library(SNPRelate)
library(dplyr)
library(SeqArray)

#### 1
GWAS_load_clumped <- function(phenopath){
  GWAS = as.data.frame(read_table2(phenopath))
  return(as.data.frame(GWAS))
}


##### 2
GWAS_extract_chromosome <- function(GWASdf, chromosome){
  #we rename the clumped file columns as the GWAS and filter by chromosome
  cn <- colnames(GWASdf)
  if("CHR" %in% cn){
    GWASchr = as.data.frame(GWASdf) %>% rename(PValue = P, Position = BP) %>% filter(CHR == chromosome)
  }else if ("Chromosome" %in% cn){
    GWASchr = as.data.frame(GWASdf) %>% filter(Chromosome == chromosome)
  }else{
    print("Unknown fields detected")
    return(c())
  }
  #we select the most significant snps by percentile
  head(GWASchr)
  # p_threshold = 1e-9
  p_threshold = 10  # in this iteration we don't need to filter it here, Genevieve filters before this
  snppos = GWASchr %>% filter(PValue < p_threshold) %>% .[["Position"]]
  return(snppos)
}  

#### 3

gds_load_chr_TOPMed <- function(chromosome){
  dirpath <- ""
  chrpath = paste("freeze.8.chr",
                  as.character(chromosome), ".pass_and_fail.gtonly.minDP0.gds", sep = "")
  gds.fn <- paste0(dirpath,chrpath)
  genofile_TOPMed <- seqOpen(gds.fn)
  return(genofile_TOPMed)
}

#### 4
#this loads the snp positions and the snp.id from the GDS file into a dataframe[,c("position", "id")]
gds_extract_snp_id_position_TOPMed <- function(genofile_TOPMed, chromosome, snppos){
  seqSetFilterChrom(genofile_TOPMed, chromosome)
  snpposall_TOPMed <- seqGetData(genofile_TOPMed, "position")
  snpidall_TOPMed = seqGetData(genofile_TOPMed, "variant.id")
  sallele_TOPMed <- seqGetData(genofile_TOPMed, "allele")
  
  snpssnppos_TOPMed <- bind_cols(snpposall_TOPMed,snpidall_TOPMed,sallele_TOPMed)
  colnames(snpssnppos_TOPMed) <- c("position", "id","allele")
  snpssnppos_TOPMed <- snpssnppos_TOPMed[which(snpssnppos_TOPMed$position %in% snppos),]
  return(snpssnppos_TOPMed)
}


#### 5
#this extracts the SNPs from the gds file into a matrix
gds_to_matrix_snp_TOPMed <- function(genofile_TOPMed, chromosome, snpidlist, sample_ids_TOPMed, reverse_allele=FALSE){
  # sample_ids_TOPMed=train$sample_id
  seqResetFilter(genofile_TOPMed)
  seqSetFilter(genofile_TOPMed, sample.id=as.character(sample_ids_TOPMed))
  seqSetFilter(genofile_TOPMed, variant.id=snpidlist$id)
  if(!reverse_allele){
    data_TOPMed <- seqGetData(genofile_TOPMed, "$dosage")
  }else{
    data_TOPMed <- seqGetData(genofile_TOPMed, "$dosage_alt")
  }
  
  samp.id = seqGetData(genofile_TOPMed, "sample.id")
  iids <- match(samp.id,sample_ids_TOPMed)
  rownames(data_TOPMed)<-sample_ids_TOPMed[iids]

  tpos = seqGetData(genofile_TOPMed, "position")
  iids <- match(tpos,snpidlist$position)
  
  colnames(data_TOPMed)<-paste(chromosome,snpidlist$position[iids], sep = "_")
  return(data_TOPMed)
}

phenotypeFile_TOPMedTrain <- "TOPMed_freeze.8_HCHS_SOL_FHS_ARIC_CHS_MESA_CARDIA_CFS_JHS_Phenotypes_Agefixed_MedicationFixedNoTCTGfix_DedupYesRel_train_v1.RData"
load(phenotypeFile_TOPMedTrain)
phenotypeFile_TOPMedTrainNR <- "TOPMed_freeze.8_HCHS_SOL_FHS_ARIC_CHS_MESA_CARDIA_CFS_JHS_Phenotypes_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_train_v1.RData"
load(phenotypeFile_TOPMedTrainNR)
phenotypeFile_TOPMedTestNR <- "TOPMed_freeze.8_HCHS_SOL_FHS_ARIC_CHS_MESA_CARDIA_CFS_JHS_Phenotypes_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_validate_v1.RData"
load(phenotypeFile_TOPMedTestNR)

outdirname <- ""
workdirname <- ""
filelist <- c( "SleepDur_clumped_SNPs_p1e-4.csv",
               "SBP_clumped_SNPs_p1e-4.csv",
               "TC_clumped_SNPs_p1e-4.csv",
               "TG_clumped_SNPs_p1e-4.csv",
               "Height_clumped_SNPs_p1e-4.csv",
               "SleepDur_unclumped_SNPs_p1e-4_RAW.csv",
               "SBP_unclumped_SNPs_p1e-4_RAW.csv",
               "TC_unclumped_SNPs_p1e-4_RAW.csv",
               "TG_unclumped_SNPs_p1e-4_RAW.csv",
               "Height_unclumped_SNPs_p1e-4_RAW.csv",
               "Height_unclumped_SNPs_p1e-4.csv",
               "SleepDur_unclumped_SNPs_p1e-4.csv",
               "SBP_unclumped_SNPs_p1e-4.csv",
               "TC_unclumped_SNPs_p1e-4.csv",
               "TG_unclumped_SNPs_p1e-4.csv")

chromosomelist = c(as.character(seq(1,22)))
chunks <- c()
for( f in filelist){ # f=filelist[1]
  chunks <- c(chunks,paste(f,chromosomelist,sep=":"))
}

print(paste0("Running ",Sys.getenv("LSB_JOBINDEX")))
args<-c(Sys.getenv("LSB_JOBINDEX"))
if (args==""){
  print("Can't get LSB_JOBINDEX")
  quit("no")
  # args<-c("1")
}

curchunk <- chunks[as.numeric(args[1])]

filename <- unlist(strsplit(curchunk,":",fixed = TRUE))[1]
chromosome <- as.numeric(unlist(strsplit(curchunk,":",fixed = TRUE))[2])

outfilename_TOPMedtrain <- paste0(outdirname,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupYesRel_genotypes_encoded_chr_",chromosome,".csv.gz"),filename))
outfilename_TOPMedtrainNR <- paste0(outdirname,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_genotypes_encoded_chr_",chromosome,".csv.gz"),filename))
outfilename_TOPMedtest <- paste0(outdirname,gsub(".csv",paste0("_TOPMedTest_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_genotypes_encoded_chr_",chromosome,".csv.gz"),filename))
if((file.exists(outfilename_TOPMedtrain))&(file.exists(outfilename_TOPMedtest))){
  print("Already processed")
  quit("no")
}

GWASfilename = paste0(workdirname,filename)

GWASdf = GWAS_load_clumped(GWASfilename)

print(GWASfilename)
snppos = GWAS_extract_chromosome(GWASdf, chromosome)
print(paste0("Chromosome ",chromosome,", extracting ",length(snppos)," SNPs"))

genofile_TOPMed = gds_load_chr_TOPMed(chromosome)

snpidlist = as.data.frame(gds_extract_snp_id_position_TOPMed(genofile_TOPMed, chromosome, snppos))
print(paste0("Was able to get the positions of ",dim(snpidlist)[1]," out of ",length(snppos),"SNPs in GWAS file"))

totidsa <- c(train$sample_id,validate$sample_id)
data_TOPMed <- gds_to_matrix_snp_TOPMed(genofile_TOPMed, chromosome, snpidlist, totidsa)

# make sure for every SNP, minor allele is 1/2 (%1,2 < %0), if not, reverse, throw away all ones with 100% 0/1
nzeros <- colSums(data_TOPMed==0)
propOnes <- 1-nzeros/nrow(data_TOPMed)
pg <- propOnes[(propOnes>0)&(propOnes<1)]
data_TOPMed <- data_TOPMed[,names(pg)]
cn <- names(pg)
for (i in 1:length(cn)){ #i=1
  if(pg[cn[i]]>0.5){
    xt <- data_TOPMed[,cn[i]]
    xt[xt==0] <- 5;xt[xt==2] <- 0;xt[xt==5] <- 2
    data_TOPMed[,cn[i]] <- xt
  }
}


data_TOPMed_train <- data_TOPMed[train$sample_id,]
data_TOPMed_trainNR <- data_TOPMed[trainNR$sample_id,]
data_TOPMed_test <- data_TOPMed[validate$sample_id,]

seqClose(genofile_TOPMed)

print("Training set size")
print(dim(data_TOPMed_train))
write.csv(data_TOPMed_train,gzfile(outfilename_TOPMedtrain), row.names = TRUE)

print("Training NR set size")
print(dim(data_TOPMed_trainNR))
write.csv(data_TOPMed_trainNR,gzfile(outfilename_TOPMedtrainNR), row.names = TRUE)

print("Test set size")
print(dim(data_TOPMed_test))
write.csv(data_TOPMed_test,gzfile(outfilename_TOPMedtest), row.names = TRUE)

