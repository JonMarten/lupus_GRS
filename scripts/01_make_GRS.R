#install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
# Note that gwas catalog is on b38 and our genetic data appears to be on b37.
setwd("/home/jm2294/projects/SLE_GRS/")
lapply(c("dplyr", "ggplot2", "rbgen", "data.table"), require, character.only = T)

sle <- read.csv("gwas_catalog/SLE_risk_SNPs_Bentham.csv", as.is=T)

sample <- fread("/scratch/bs446/BAKEOFF151001/genetic_data/imputed/somalogic/round_all_m0.001_i0.3_PostProtein_IBDQC/SOMAS_round_all_SOMAS_round_all_impute_1_interval_filtered2.sample", data.table=F)
sample <- sample[-1,]

i = 1
for(i in 1:22)
bgenPath <- paste0("/home/jm2294/GENETIC_DATA/INTERVAL/bgen/SOMAS_round_all_SOMAS_round_all_impute_",
                   sle$CHR_ID[i],
                   "_interval_filtered2_reworked.bgen")

data <- bgen.load(bgenPath, rsids = sle$SNPS[i])
data$samples <- sample$ID_1
colnames(data$ploidy) <- sample$ID_1