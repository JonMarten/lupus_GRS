# Concatenate Olink panel files
# Parameters
samplePath <- "/home/jp549/post-doc/genetics/r2-test/sample_files/"
samplePrefix <- "o5000-"
sampleSuffix <- "-outlier_in-r2.sample"

setwd("/home/jm2294/projects/SLE_GRS/OLINK")
lapply(c("dplyr", "ggplot2", "data.table", "cowplot"), require, character.only = T)

panels <- c("inf1","cvd2","cvd3")
protPanel <- data.frame()
pheList <- list()
for(j in 1:length(panels)){
  panel <- panels[j]
  sample <- fread(data.table = F, paste0(samplePath, samplePrefix,panel,sampleSuffix))
  sample <- sample[-1,]
  sample <- apply(sample, MARGIN = 2, FUN = as.numeric) %>% data.frame(stringsAsFactors = F)
  pheList[[j]] <- sample
  names(pheList)[j] <- panel
  
  # Make data frame with protein-panel assignments
  tempPanel <- data.frame("protein" = names(sample),"panel" = panel, stringsAsFactors = F) %>% 
    filter(grepl("___", protein))
  protPanel <- rbind(protPanel, tempPanel)
  
  rm(sample, tempPanel, panel)
}

# Identify proteins included on multiple panels
dupeProts <- protPanel %>% group_by(protein) %>% filter(n() > 1) %>% data.frame(stringsAsFactors = F)
which(names(pheList$inf1) %in% dupeProts$protein)

names(pheList$inf1)[which(names(pheList$inf1) %in% dupeProts$protein)] <- paste0(names(pheList$inf1)[which(names(pheList$inf1) %in% dupeProts$protein)],"_inf1")
names(pheList$cvd2)[which(names(pheList$cvd2) %in% dupeProts$protein)] <- paste0(names(pheList$cvd2)[which(names(pheList$cvd2) %in% dupeProts$protein)],"_cvd2")
names(pheList$cvd3)[which(names(pheList$cvd3) %in% dupeProts$protein)] <- paste0(names(pheList$cvd3)[which(names(pheList$cvd3) %in% dupeProts$protein)],"_cvd3")

phe <- full_join(pheList$inf1, pheList$cvd2)
phe <- full_join(phe,pheList$cvd3)
