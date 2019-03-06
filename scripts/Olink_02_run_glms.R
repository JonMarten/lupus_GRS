# Run GLMs on GRS
lapply(c("data.table", "dplyr", "stringr","ggplot2","cowplot", "tidyr"), require, character.only = T)
setwd("/home/jm2294/projects/SLE_GRS/OLINK")

# Read in data
olink <- fread("INTERVAL_Olink_inf1_cvd2_cvd3_merged.csv", data.table = F)
grs <- fread("INTERVAL_olink_SLE_GRS.csv",data.table = F)
manifest <- fread("olink_protein_panel_manifest.csv", data.table = F)

# Assign to GRS categories
bins <- seq(min(grs$GRS)-0.001,max(grs$GRS)+0.001,length.out = 6)
grs <- grs %>%
  mutate(GRScat = cut(GRS, breaks = bins))

# get largest bin to use as reference for ORs
refCat <- names(which.max(table(grs$GRScat)))
grs$GRScat <- as.factor(grs$GRScat)
grs$GRScat <- relevel(grs$GRScat, ref = refCat)

phe <- full_join(grs,olink, by = "id")

glmCoefsAll <- data.frame()
for( i in 1:nrow(manifest)){
  cat(paste0("\nProcessing ",i," of ",nrow(manifest),": ",manifest$protName[i]),"\n")
  panel <- manifest$panel[i]
  a <- paste0(names(phe)[4:6], collapse = " + ")
  b <- paste0(names(phe)[10:30], collapse = " + ")
  c <- paste0("plate_", panel)
  covars <- paste0(a," + ", b," + ", c, collapse = " + ")
   
  form <- formula(paste0(manifest$protName[i], " ~ GRScat + ", covars))
  print(form)
  #dat <- phe[,c(3,i+3)]
  protGLM <- glm(data = phe, formula = form)
  coefs <- data.frame(coef(summary(protGLM)))
  coefs <- coefs %>%
    mutate(GRScat = rownames(coefs),
           protein = manifest$protName[i]) %>%
    select(protein, 
           GRScat, 
           estimate = Estimate,
           SE = Std..Error,
           t = t.value,
           P = Pr...t..)
  coefTall <- gather(coefs[,c("estimate","SE","t","P")])
  coefTall$protein <- manifest$protName[i]
  coefTall$GRScat <- coefs$GRScat
  coefTall <- coefTall %>%
    select(protein, GRScat, key, value)
  glmCoefsAll <- rbind(glmCoefsAll, coefTall)
  #rm(plotname, bxplot,dat,form,protGLM,coefs,coefTall)
  rm(form,protGLM,coefs,coefTall)
}
sig <- filter(glmCoefsAll, key =="P" & value < 0.05/nrow(manifest))

# FDR correction. N.B. THIS HARDCODES GRS cats so will be broken for a different trait or threshold. Should be fixed.
Qvals <- glmCoefsAll %>%
  filter(key == "P") %>%
  spread(key = "GRScat", value = "value") 

Qvals[,-c(1,2)] <- apply(Qvals[,-c(1,2)], MARGIN = 2, FUN = p.adjust)
QvalsT <- gather(Qvals[,-c(1:2)],key = "GRScat")
QvalsT$protein <- Qvals$protein
QvalsT$key <- "Q"
glmCoefsAll <- bind_rows(glmCoefsAll, QvalsT)
fwrite(glmCoefsAll, file = "olink_glm_coefs_all.csv", sep=",")

# Make boxplots for all proteins with a significant association
dir.create("boxplots") # Throws a warning if dir already exists but does not terminate script
phe$GRScat <- factor(phe$GRScat, 
                     levels = sort(levels(phe$GRScat))) # reorder factor to match increasing GRS. 
sigProts <- glmCoefsAll %>% 
  filter(key == "Q", value < 0.05) %>%
  pull(protein) %>%
  unique()

for(i in 1:length(sigProts)){
  protCol <- which(names(phe) == sigProts[i])
  bp <- ggplot(phe, aes(x = GRScat, y = phe[,protCol])) +
    geom_boxplot() +
    geom_jitter(alpha = 0.1, size = 0.2, width = 0.2) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(x = "GRS category", y = sigProts[i])
  ggsave(paste0("boxplots/boxplot_",sigProts[i],".png"), plot = bp)
}

# Continuous GRS models
glmCoefsAllCont <- data.frame()
for( i in 1:nrow(manifest)){
  #for( i in 1:10){
  cat(paste0("\nProcessing ",i," of ",nrow(manifest),": ",manifest$protName[i]))
  # glm
  form <- formula(paste0(manifest$protName[i], " ~ GRS + ", covars))
  protGLM <- glm(data = phe, formula = form)
  coefs <- data.frame(coef(summary(protGLM)))
  coefs <- coefs %>%
    mutate(GRS = rownames(coefs),
           protein = manifest$protName[i]) %>%
    select(protein, 
           GRS, 
           estimate = Estimate,
           SE = Std..Error,
           t = t.value,
           P = Pr...t..)
  coefTall <- gather(coefs[,c("estimate","SE","t","P")])
  coefTall$protein <- manifest$protName[i]
  coefTall$GRS <- coefs$GRS
  coefTall <- coefTall %>%
    select(protein, GRS, key, value)
  glmCoefsAllCont <- rbind(glmCoefsAllCont, coefTall)
  #rm(plotname, bxplot,dat,form,protGLM,coefs,coefTall)
  rm(form,protGLM,coefs,coefTall)
}

# FDR correction
Qvals <- glmCoefsAllCont %>%
  filter(key == "P") %>%
  spread(key = "GRS", value = "value") 

Qvals[,-c(1,2)] <- apply(Qvals[,-c(1,2)], MARGIN = 2, FUN = p.adjust)

QvalsT <- gather(Qvals[,3:ncol(Qvals)],key = "key", value = "value")
QvalsT$protein <- Qvals$protein
QvalsT$key <- "Q"
glmCoefsAllCont <- bind_rows(glmCoefsAllCont, QvalsT)
fwrite(glmCoefsAllCont, file = "glm_coefs_all_continuous.csv", sep=",")