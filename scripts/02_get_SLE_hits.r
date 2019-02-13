library(data.table)
library(dplyr)
sle <- fread("data/SLE_gwas-association-downloaded_2019-02-13-EFO_0002690-withChildTraits.tsv", data.table=F)

#sle %>% filter(PUBMEDID=="29848360") # Julia et al 2018 (only 5 new hits in Spanish)
sle_benth <- sle %>% 
  filter(`STUDY ACCESSION`=="GCST003156") %>%  # Bentham et al 2015, 7219 cases and 15991 controls
  select(SNPS, CHR_ID, CHR_POS, REGION, MAPPED_GENE, `STRONGEST SNP-RISK ALLELE`, `OR or BETA`, `P-VALUE`) %>%
  arrange(CHR_ID, CHR_POS)

  



