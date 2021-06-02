library(readr)
library(dplyr)
library(glmnet)
library(doMC)
registerDoMC(cores=2)

#load data
PCADJ = TRUE
fstart <- 12

datadirname <- "/yourdir/data/"
load(paste0(datadirname,"TOPMed_freeze.8_HCHS_SOL_FHS_ARIC_CHS_MESA_CARDIA_CFS_JHS_Phenotypes_Agefixed_MedicationFixedNoTCTGfix_DedupYesRel_train_v1.RData"))

pdata <- train

covariates = c("Age","Sex","study","race")
if (PCADJ){
  covariates = c(covariates,paste0("EV",seq(1,5,1)))
}

#adjust for covariates
print("Adjusting raw data for Covariates...")
phenotypesN <- colnames(pdata)
for (i in fstart:length(phenotypesN)){ # i <- fstart
  pheno.1 = phenotypesN[i]
  print(pheno.1)
  model.formula <- as.formula(paste(paste(pheno.1,"~"), paste0(covariates,collapse="+")))
  if(length(which(!is.na(pdata[pheno.1])))<50){
    print("All NA, Skippig")
    pdata[,i]<-NA
    next
  }
  r2 <- pdata[,i]
  result <- tryCatch({
    mod <- lm(model.formula, data = pdata[,c(pheno.1,covariates)],na.action="na.exclude")
    r2 <- resid(mod,na.action="na.exclude")
  }, warning = function(war) {
    print("Adjusting failed")
  }, error = function(err) {
    print("Adjusting failed")
  }) # END tryCatch
  
  pdata[,i]<-r2
}

# #ranknorm
# print("Rank-normalizing...")
# for (i in fstart:length(pdata)) {
#   x <-pdata[,i]
#   pdata[,i] <- scale(qnorm(rank(x, ties.method = "random",na.last="keep")/(length(x)+1)))
# }
rownames(pdata) <- pdata$sample_id
train <- pdata

###################################
datadirname <- "/yourdir/data/"
load(paste0(datadirname,"TOPMed_freeze.8_HCHS_SOL_FHS_ARIC_CHS_MESA_CARDIA_CFS_JHS_Phenotypes_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_train_v1.RData"))

pdata <- trainNR

covariates = c("Age","Sex","study","race")
if (PCADJ){
  covariates = c(covariates,paste0("EV",seq(1,5,1)))
}

#adjust for covariates
print("Adjusting raw data for Covariates...")
phenotypesN <- colnames(pdata)
for (i in fstart:length(phenotypesN)){ # i <- fstart
  pheno.1 = phenotypesN[i]
  print(pheno.1)
  model.formula <- as.formula(paste(paste(pheno.1,"~"), paste0(covariates,collapse="+")))
  if(length(which(!is.na(pdata[pheno.1])))<50){
    print("All NA, Skippig")
    pdata[,i]<-NA
    next
  }
  r2 <- pdata[,i]
  result <- tryCatch({
    mod <- lm(model.formula, data = pdata[,c(pheno.1,covariates)],na.action="na.exclude")
    r2 <- resid(mod,na.action="na.exclude")
  }, warning = function(war) {
    print("Adjusting failed")
  }, error = function(err) {
    print("Adjusting failed")
  }) # END tryCatch
  
  pdata[,i]<-r2
}

#ranknorm
print("Rank-normalizing...")
for (i in fstart:length(pdata)) {
    x <-pdata[,i]
    pdata[,i] <- scale(qnorm(rank(x, ties.method = "random",na.last="keep")/(length(x)+1)))
}
rownames(pdata) <- pdata$sample_id
trainNR <- pdata


########################################## Load the genotypes ###########################
workdirname_data <- "/yourdir/SNP_Lists/extract_distribAgeMedNoRel/"
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
parts <- c()
for(f in filelist){ # f=filelist[6]
  tchunks <- paste(f,chromosomelist,sep=":")
  chunks <- c(chunks,tchunks)
  if(length(grep("unclumped",f,fixed = T))>0){
    numparts <- 5
  }else{
    numparts <- 3
  }
  parts <- c(parts,split(tchunks, sort((1:length(tchunks)) %% numparts)))
}

print(paste0("Running ",Sys.getenv("LSB_JOBINDEX")))
args<-c(Sys.getenv("LSB_JOBINDEX"))
if (args==""){
  stop("Can't get LSB_JOBINDEX")
  # args<-c("8")
}

cparts <- parts[as.numeric(args[1])][[1]]
curchunk <- cparts[1]
filename <- unlist(strsplit(curchunk,":",fixed = TRUE))[1]
if(length(grep("Height_",filename,fixed = TRUE))>0){
  pheno.1 <- "height_baseline_1"
}else if(length(grep("TG_",filename,fixed = TRUE))>0){
  pheno.1 <- "triglycerides_1"
}else if(length(grep("TC_",filename,fixed = TRUE))>0){
  pheno.1 <- "total_cholesterol_1"
}else if(length(grep("SleepDur_",filename,fixed = TRUE))>0){
  pheno.1 <- "sleep_duration_1"
}else if(length(grep("SBP_",filename,fixed = TRUE))>0){
  pheno.1 <- "bp_systolic_1"
}else{
  stop("Unrecognized file, quitting...")
}

outfilename <- paste0(workdirname_data,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupYesRel_genotypes_encoded_LASSO_part_",args[1],"_v2.csv.gz"),filename))
print(paste0("Calculating ",outfilename," for phenotype: ",pheno.1))
if(file.exists(outfilename)){
  stop("Already processed...Exiting...")
}
print("All train set")

############################# RelatedPeople #####################################
#load snps in the region for the training samples
genodata_train <- c()
for (i in 1:length(cparts)){ #i=1
  #   args<-c(i)
  curchunk <- cparts[i]
  
  filename <- unlist(strsplit(curchunk,":",fixed = TRUE))[1]
  chromosome <- as.numeric(unlist(strsplit(curchunk,":",fixed = TRUE))[2])
  print(paste0("Loading ",filename," chromosome ",chromosome)) 
  
  infilename_TOPMedtrain <- paste0(workdirname_data,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupYesRel_genotypes_encoded_chr_",chromosome,".csv.gz"),filename))
  
  #load snps in the region for the training samples
  genodata_tr <- read.table(file=gzfile(infilename_TOPMedtrain),sep=",",header = TRUE, row.names = 1,check.names = FALSE)
  
  if(length(genodata_train)<1){
    genodata_train <- genodata_tr
  }else{
    genodata_train <- cbind(genodata_train,genodata_tr)
  }  
}

#################################################### all races
res <- list()

ids <- as.character(train$sample_id[!is.na(train[,pheno.1])])
ids <- ids[ids %in% rownames(genodata_train)]
x_train <- as.matrix(genodata_train[ids,])
y_train <- train[ids,pheno.1]

# get the SNPs
res$numSNPs <- length(colnames(genodata_train))

# lasso (number of features that survive + MSE)
sst <- Sys.time()
lasso_cv <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10, parallel = TRUE)
print(Sys.time() - sst)

tlids <- which(lasso_cv$nzero>5)
if(length(tlids)<40){
  tlids <- c(tlids,1:(40-length(tlids)))  
}else{
  tlids <- tlids[1:40]
}
# testlambdas <- c(lasso_cv$lambda.min,sort(lasso_cv$lambda[tlids]),10),seq(from = 0.0001, to = 0.1, length.out =  39))
testlambdas <- unique(c(lasso_cv$lambda.min,lasso_cv$lambda[tlids]))

for( ll in testlambdas){ # ll = testlambdas[1]
  best_lambda <- ll
  best_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
  lasso_coef <- coef(best_lasso)
  lasso_res <- data.frame(name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1], coefficient = lasso_coef@x)
  res$lasso_NZ <- dim(lasso_res)[1]
  print(paste0("LASSO number of nonzero parameters:",res$lasso_NZ," out of ",res$numSNPs," SNPs"))
  res$lasso_Names <- paste0(lasso_res$name,collapse = ";")
  res$lasso_Coefs <- paste0(lasso_res$coefficient,collapse = ";")
  
  ot <- gsub(".csv.gz",paste0("_ALL_NZ_",res$lasso_NZ,"_Lambda_",round(ll,3),".csv.gz"),outfilename)
  write.table(res,file = gzfile(ot),sep=",", row.names = FALSE)
}


#################################################### races
allrace <- c("Black","White","HL")
for(frace in allrace){ # frace <- "Black"
  res <- list()
  
  ids <- as.character(train$sample_id[!is.na(train[,pheno.1]) & (train[,"race"] == frace)])
  ids <- ids[ids %in% rownames(genodata_train)]
  x_train <- as.matrix(genodata_train[ids,])
  y_train <- train[ids,pheno.1]
  
  # get the SNPs
  res$numSNPs <- length(colnames(genodata_train))
  
  # lasso (number of features that survive + MSE)
  sst <- Sys.time()
  lasso_cv <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10, parallel = TRUE)
  print(Sys.time() - sst)
  
  tlids <- which(lasso_cv$nzero>5)
  if(length(tlids)<40){
    tlids <- c(tlids,1:(40-length(tlids)))  
  }else{
    tlids <- tlids[1:40]
  }
  # testlambdas <- c(lasso_cv$lambda.min,sort(lasso_cv$lambda[tlids]),10),seq(from = 0.0001, to = 0.1, length.out =  39))
  testlambdas <- unique(c(lasso_cv$lambda.min,lasso_cv$lambda[tlids]))
  
  for( ll in testlambdas){ # ll = testlambdas[1]
    best_lambda <- ll
    best_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
    lasso_coef <- coef(best_lasso)
    lasso_res <- data.frame(name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1], coefficient = lasso_coef@x)
    res$lasso_NZ <- dim(lasso_res)[1]
    print(paste0("LASSO number of nonzero parameters:",res$lasso_NZ," out of ",res$numSNPs," SNPs"))
    res$lasso_Names <- paste0(lasso_res$name,collapse = ";")
    res$lasso_Coefs <- paste0(lasso_res$coefficient,collapse = ";")
    
    ot <- gsub(".csv.gz",paste0("_",frace,"_NZ_",res$lasso_NZ,"_Lambda_",round(ll,3),".csv.gz"),outfilename)
    write.table(res,file = gzfile(ot),sep=",", row.names = FALSE)
  }
}


############################# Non RelatedPeople #####################################
print("Non-related train set")
rm(genodata_train);rm(train)

outfilename <- paste0(workdirname_data,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_genotypes_encoded_LASSO_part_",args[1],"_v2.csv.gz"),filename))
print(paste0("Calculating ",outfilename," for phenotype: ",pheno.1))
if(file.exists(outfilename)){
  stop("Already processed...Exiting...")
}

#load snps in the region for the training samples
genodata_train <- c()
for (i in 1:length(cparts)){ #i=1
  #   args<-c(i)
  curchunk <- cparts[i]
  
  filename <- unlist(strsplit(curchunk,":",fixed = TRUE))[1]
  chromosome <- as.numeric(unlist(strsplit(curchunk,":",fixed = TRUE))[2])
  print(paste0("Loading ",filename," chromosome ",chromosome)) 
  
  infilename_TOPMedtrain <- paste0(workdirname_data,gsub(".csv",paste0("_TOPMedTrain_Agefixed_MedicationFixedNoTCTGfix_DedupNoRel_genotypes_encoded_chr_",chromosome,".csv.gz"),filename))
  
  #load snps in the region for the training samples
  genodata_tr <- read.table(file=gzfile(infilename_TOPMedtrain),sep=",",header = TRUE, row.names = 1,check.names = FALSE)
  
  if(length(genodata_train)<1){
    genodata_train <- genodata_tr
  }else{
    genodata_train <- cbind(genodata_train,genodata_tr)
  }  
}

#################################################### all races
res <- list()

ids <- as.character(trainNR$sample_id[!is.na(trainNR[,pheno.1])])
ids <- ids[ids %in% rownames(genodata_train)]
x_train <- as.matrix(genodata_train[ids,])
y_train <- trainNR[ids,pheno.1]

# get the SNPs
res$numSNPs <- length(colnames(genodata_train))

# lasso (number of features that survive + MSE)
sst <- Sys.time()
lasso_cv <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10, parallel = TRUE)
print(Sys.time() - sst)

tlids <- which(lasso_cv$nzero>5)
if(length(tlids)<40){
  tlids <- c(tlids,1:(40-length(tlids)))  
}else{
  tlids <- tlids[1:40]
}
# testlambdas <- c(lasso_cv$lambda.min,sort(lasso_cv$lambda[tlids]),10),seq(from = 0.0001, to = 0.1, length.out =  39))
testlambdas <- unique(c(lasso_cv$lambda.min,lasso_cv$lambda[tlids]))

for( ll in testlambdas){ # ll = testlambdas[1]
  best_lambda <- ll
  best_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
  lasso_coef <- coef(best_lasso)
  lasso_res <- data.frame(name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1], coefficient = lasso_coef@x)
  res$lasso_NZ <- dim(lasso_res)[1]
  print(paste0("LASSO number of nonzero parameters:",res$lasso_NZ," out of ",res$numSNPs," SNPs"))
  res$lasso_Names <- paste0(lasso_res$name,collapse = ";")
  res$lasso_Coefs <- paste0(lasso_res$coefficient,collapse = ";")
  
  ot <- gsub(".csv.gz",paste0("_ALL_NZ_",res$lasso_NZ,"_Lambda_",round(ll,3),".csv.gz"),outfilename)
  write.table(res,file = gzfile(ot),sep=",", row.names = FALSE)
}


#################################################### Black
allrace <- c("Black","White","HL")
for(frace in allrace){ # frace <- "Black"
  print(frace)
  res <- list()
  
  ids <- as.character(trainNR$sample_id[!is.na(trainNR[,pheno.1]) & (trainNR[,"race"] == frace)])
  ids <- ids[ids %in% rownames(genodata_train)]
  x_train <- as.matrix(genodata_train[ids,])
  y_train <- trainNR[ids,pheno.1]
  
  # get the SNPs
  res$numSNPs <- length(colnames(genodata_train))
  
  # lasso (number of features that survive + MSE)
  sst <- Sys.time()
  lasso_cv <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10, parallel = TRUE)
  print(Sys.time() - sst)
  
  tlids <- which(lasso_cv$nzero>5)
  if(length(tlids)<40){
    tlids <- c(tlids,1:(40-length(tlids)))  
  }else{
    tlids <- tlids[1:40]
  }
  # testlambdas <- c(lasso_cv$lambda.min,sort(lasso_cv$lambda[tlids]),10),seq(from = 0.0001, to = 0.1, length.out =  39))
  testlambdas <- unique(c(lasso_cv$lambda.min,lasso_cv$lambda[tlids]))
  
  for( ll in testlambdas){ # ll = testlambdas[1]
    best_lambda <- ll
    best_lasso <- glmnet(x_train, y_train, alpha = 1, lambda = best_lambda)
    lasso_coef <- coef(best_lasso)
    lasso_res <- data.frame(name = lasso_coef@Dimnames[[1]][lasso_coef@i + 1], coefficient = lasso_coef@x)
    res$lasso_NZ <- dim(lasso_res)[1]
    print(paste0("LASSO number of nonzero parameters:",res$lasso_NZ," out of ",res$numSNPs," SNPs"))
    res$lasso_Names <- paste0(lasso_res$name,collapse = ";")
    res$lasso_Coefs <- paste0(lasso_res$coefficient,collapse = ";")
    
    ot <- gsub(".csv.gz",paste0("_",frace,"_NZ_",res$lasso_NZ,"_Lambda_",round(ll,3),".csv.gz"),outfilename)
    write.table(res,file = gzfile(ot),sep=",", row.names = FALSE)
  }
}


