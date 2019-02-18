#install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
# Note that gwas catalog is on b38 and our genetic data appears to be on b37.
# N.B. bgen data array is snps x samples x genotypes
setwd("/home/jm2294/projects/SLE_GRS/")
lapply(c("dplyr", "ggplot2", "rbgen", "data.table"), require, character.only = T)

sle <- read.csv("gwas_catalog/SLE_risk_SNPs_Bentham.csv", as.is=T)

sample <- fread("/scratch/bs446/BAKEOFF151001/genetic_data/imputed/somalogic/round_all_m0.001_i0.3_PostProtein_IBDQC/SOMAS_round_all_SOMAS_round_all_impute_1_interval_filtered2.sample", data.table=F)
sample <- sample[-1,]

# FUNCTION: Convert bgen AA:AB:BB format to B dosage format
bgenToDose <- function(bgen){
  nsnps <- nrow(bgen$variants)
  nids <- length(bgen$samples)
  out <- matrix(nrow = nids, ncol = nsnps)
  rownames(out) <- bgen$samples
  colnames(out) <- rownames(bgen$variants)
  for(i in 1:nsnps){
    out[,i] <- bgen$data[i,,2] + 2*bgen$data[i,,3]
  }
  out2 <- list("variants" = bgen$variants, "dosage" = out)
  return(out2)
}

firstRow <- TRUE
for(i in 1:22){
  bgenPath <- paste0("/home/jm2294/GENETIC_DATA/INTERVAL/bgen/SOMAS_round_all_SOMAS_round_all_impute_",
                     i,
                     "_interval_filtered2_reworked.bgen")
  snpRow <- which(sle$CHR_ID == i)
  if(length(snpRow) == 0){next}   # Skip to next chromosome if no SNPs on this one
  bgen <- bgen.load(bgenPath, rsids = sle$SNPS[snpRow])
  dose <- bgenToDose(bgen)
  if(firstRow) {
    snps <- dose$variants
    dosage <- dose$dosage
    firstRow <- FALSE
  } else {
    snps <- rbind(snps, dose$variants)
    dosage <- cbind(dosage, dose$dosage)
  }
  rm(bgen)
  rm(dose)
}


 
  # Read in genotypes for hit SNPs from bgen files.
  sample <- fread("/gpfs/igmmfs01/eddie/UK-BioBank-proj19655/genotypes/imputed/ukb19655_imp_chr1_v3_s487395.sample", data.table=F)
  sample <- sample[-1,]
  
  # Loop over all chromosomes. This file needs UK Biobank imputed genotypes
  for(i in 1:22){
    bgen <- bgen.load(paste0("/gpfs/igmmfs01/eddie/UK-BioBank-proj19655/genotypes/imputed/ukbb_proj19655_chr",i,".bgen"), rsids = urateHits$SNP)
    dose <- bgenToDose(bgen)
    if(i == 1) {
      snps <- dose$variants
      dosage <- dose$dosage
    } else {
      snps <- rbind(snps, dose$variants)
      dosage <- cbind(dosage, dose$dosage)
    }
    rm(bgen)
    rm(dose)
  }
  rownames(dosage) <- sample$ID_1
  

  
  
    
    