#install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
# Note that gwas catalog is on b38 and our genetic data appears to be on b37.
# N.B. bgen data array is snps x samples x genotypes
# Parameters
genPath <- "/scratch/bp406/data_sets/interval_subset_olink/genotype_files/unrelated_4994_pihat_0.1875_autosomal_imputed_info_0.4_phwe_1e-4_filtered/per_chr/"
genPrefix <- "interval_subset_olink.wga_imputed.chr_"
genSuffix <- ".unrelated_4994.pihat_0.1875_info_0.4_phwe_1e-4_filtered.bgen"
samplePath <- "/home/jp549/post-doc/genetics/r2-test/sample_files/"
samplePrefix <- "o5000-"
sampleSuffix <- "-outlier_in-r2.sample"

setwd("/home/jm2294/projects/SLE_GRS/OLINK")
lapply(c("dplyr", "ggplot2", "rbgen", "data.table", "cowplot"), require, character.only = T)

sle <- read.csv("../gwas_catalog/SLE_risk_SNPs_Bentham.csv", as.is=T)
sle <- sle %>%
  select(rsid = SNPS,
         chr = CHR_ID,
         pos_38 = CHR_POS,
         region = REGION,
         gene = MAPPED_GENE,
         effAllele = effAll,
         OR = OR.or.BETA,
         beta,
         P = P.VALUE)
# All SNPs should be set to the risk-increasing direction. This has already been done by GWAS catalog in this case.
panels <- c("inf1","cvd2","cvd3")
i <- 1
panel <- panels[i]
sample <- fread(data.table = F, paste0(samplePath, samplePrefix,panel,sampleSuffix))
sample <- sample[-1,]
sample <- apply(sample, MARGIN = 2, FUN = as.numeric)

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

# Loop over all 22 autosomes and save B allele dosages for each individual
firstRow <- TRUE
for(i in 1:22){
  bgenPath <- paste0("/home/jm2294/GENETIC_DATA/INTERVAL/bgen/SOMAS_round_all_SOMAS_round_all_impute_",
                     i,
                     "_interval_filtered2_reworked.bgen")
  snpRow <- which(sle$chr == i)
  if(length(snpRow) == 0){next}   # Skip to next chromosome if no SNPs on this one
  bgen <- bgen.load(bgenPath, rsids = sle$rsid[snpRow])
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
snps <- snps %>%
  mutate(rsid = as.character(rsid),
         allele0 = as.character(allele0),
         allele1 = as.character(allele1))

# Add proper sample ids as these are not embedded in the bgen files
rownames(dosage) <- as.character(sample$ID_1)

# Remove X chr (temporary)
sle <- sle %>%
  filter(chr != "X")

# Reorder GWAS summary statistics to match dosage order 
sle <- sle[match(snps$rsid, sle$rsid),]
identical(sle$rsid, snps$rsid) # check TRUE

swapAllele <- which(toupper(snps$allele1) != toupper(sle$effAllele))
dosage[,swapAllele] <- 2 - dosage[,swapAllele]

# Weight dosage by beta and calculate sum
dosewgt <- data.frame(dosage)
for(i in 1:nrow(sle)){
  dosewgt[,i] <- dosewgt[,i] * sle$beta[i]
}
GRS <-  apply(dosewgt, MARGIN = 1, sum)
grs <- data.frame("id" = rownames(dosage), GRS)

# Generate kernel density plot of GRS
timestamp <- strftime(Sys.time(), format = "%Y-%m-%d_%H:%M:%S" )
kernelDensity <- ggplot(grs, aes(GRS)) +
  geom_density()

save_plot(paste0("SLE_GRS_kdens_",timestamp,".png"),
          kernelDensity)

# Save GRS for later use
fwrite(grs, file = "INTERVAL_soma_SLE_GRS.csv")
