library(data.table)
library(dplyr)
library(stringr)
sle <- fread("SLE_gwas-association-downloaded_2019-02-13-EFO_0002690-withChildTraits.tsv", data.table=F)
snplist <- fread("sle_benth_2015_SNPlist.txt", data.table=F, h = F) %>% pull(V1)# List of 43 SNPs reported in Bentham et al 2015

# Get 43 SNPs reported in Bentham et al. Where there are duplicates, keep only the one with the lowest p-value 
sle_benth <- sle %>%
  filter(PUBMEDID == "26502338" & SNPS %in% snplist) %>% 
  arrange(`P-VALUE`)
sle_benth <- sle_benth[!duplicated(sle_benth$SNPS),]
sle_benth <- sle_benth %>% 
  select(SNPS, CHR_ID, CHR_POS, REGION, MAPPED_GENE, `STRONGEST SNP-RISK ALLELE`, `OR or BETA`, `P-VALUE`) %>%
  arrange(CHR_ID, CHR_POS) %>%
  mutate(beta = log(`OR or BETA`),
         effAll = str_split_fixed(`STRONGEST SNP-RISK ALLELE`, "-", 2)[,2])

write.csv(sle_benth, quote=T, row.names=F, file = "SLE_risk_SNPs_Bentham.csv")
