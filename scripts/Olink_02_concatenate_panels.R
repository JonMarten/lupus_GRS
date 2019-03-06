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

fwrite(protPanel, file = "olink_protein_panel_manifest.csv")

# Identify proteins included on multiple panels
dupeProts <- protPanel %>% group_by(protein) %>% filter(n() > 1) %>% data.frame(stringsAsFactors = F)
which(names(pheList$inf1) %in% dupeProts$protein)

names(pheList$inf1)[which(names(pheList$inf1) %in% dupeProts$protein)] <- paste0(names(pheList$inf1)[which(names(pheList$inf1) %in% dupeProts$protein)],"_inf1")
names(pheList$cvd2)[which(names(pheList$cvd2) %in% dupeProts$protein)] <- paste0(names(pheList$cvd2)[which(names(pheList$cvd2) %in% dupeProts$protein)],"_cvd2")
names(pheList$cvd3)[which(names(pheList$cvd3) %in% dupeProts$protein)] <- paste0(names(pheList$cvd3)[which(names(pheList$cvd3) %in% dupeProts$protein)],"_cvd3")

inf1 <- pheList$inf1 %>%
  rename(plate_inf1 = plate)
cvd2 <- pheList$cvd2 %>%
  rename(plate_cvd2 = plate)
cvd3 <- pheList$cvd3 %>% 
  rename(plate_cvd3 = plate)

# Check columns match other than NAs
for(i in 1:28) {
  print(names(inf1)[i])
  if(identical(inf1[,i], cvd2[,i]) & identical(cvd2[,i], cvd3[,i])) {
    print("Match") 
  } else {
    print(which(inf1[,i] != cvd2[,i]))
  }
}

# Coalesce to remove NAs
phe <- data.frame(matrix(ncol = 28, nrow = nrow(inf1)))
for(i in 1:28){
  phe[,i] <- coalesce(inf1[,i], cvd2[,i], cvd3[,i])
}
names(phe) <- names(inf1)[1:28]

# Combine into single data frame
phe <- cbind(phe, cvd2$plate_cvd2, cvd3$plate_cvd3)
phe <- phe %>%
  select(id = ID_1,
         age,
         sexPulse,
         season,
         plate_inf1,
         plate_cvd2 = "cvd2$plate_cvd2",
         plate_cvd3 = "cvd3$plate_cvd3",
         bleed_to_process_time,
         PC1:PC20)

phe <- cbind(phe, inf1[,29:120], cvd2[,29:120], cvd3[,29:120])

fwrite(phe, "INTERVAL_Olink_inf1_cvd2_cvd3_merged.csv")
