# Read in somalogic protein data and combine into a single matrix
lapply(c("data.table", "dplyr", "stringr"), require, character.only = T)
setwd("/home/jm2294/projects/SLE_GRS")

protPath <- "/scratch/bs446/BAKEOFF151001/protein_data/somalogic/for_gwas/round_all_SNPTEST/Perprotein_pheno/"
prot <- fread(paste0(protPath,"round_all_IFNA2.3497.13.2.sample"), data.table=F)
names(prot) <- c("id","id2","INFNA2")

protList <- fread("somalogic_proteins_filenames.csv", data.table = F)

for( i in 1:nrow(protList)){
  cat(paste0("\nProcessing ",i," of ",nrow(protList),": ",protList$protein[i]))
  prot <- fread(paste0(protPath,protList$filename[i]), data.table=F)
  if(i == 1){
    allProt <- data.frame(prot$V1, prot$V3)
    names(allProt) <- c("id",protList$protein[i])
  } else {
    allProt[,i+1] <- prot$V3
    names(allProt)[i+1] = protList$protein[i]
  }
  rm(prot)
}

fwrite(allProt, sep = ",", row.names = F, file = "Somalogic_all_proteins_agesexrntrans.csv")
