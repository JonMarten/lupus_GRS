# Run GLMs on GRS
lapply(c("data.table", "dplyr", "stringr","ggplot2","cowplot", "tidyr"), require, character.only = T)
setwd("/home/jm2294/projects/SLE_GRS")

allProt <- fread("Somalogic_all_proteins_agesexrntrans.csv", data.table = F)
protList <- fread("somalogic_proteins_filenames.csv", data.table = F)
rownames(allProt) <- allProt$id

# Assign to GRS categories
grs <- fread("INTERVAL_soma_SLE_GRS.csv", data.table=F)
bins <- seq(min(grs$GRS)-0.001,max(grs$GRS)+0.001,length.out = 6)
grs <- grs %>%
  mutate(GRScat = cut(GRS, breaks = bins))

# get largest bin to use as reference for ORs
refCat <- names(which.max(table(grs$GRScat)))
grs$GRScat <- as.factor(grs$GRScat)
grs$GRScat <- relevel(grs$GRScat, ref = refCat)

phe <- full_join(grs,allProt, by = "id")

# Make boxplots
#dir.create("boxplots")

glmCoefsAll <- data.frame()
for( i in 1:nrow(protList)){
#for( i in 1:10){
  cat(paste0("\nProcessing ",i," of ",nrow(protList),": ",protList$protein[i]))
  #plotname <- paste0("boxplots/",protList$protein[i],"_boxplot.png")
  #bxplot <- ggplot(phe, aes(x = GRScat, y = phe[,i+3])) +
  #  geom_jitter(alpha = 0.1, size = 0.2) +
  #  geom_boxplot() +
  #  theme(axis.text.x = element_text(angle = 90)) +
  #  labs(x = "GRS category", y = protList$protein[i])
  #ggsave(bxplot, filename = plotname)
  
  # glm
  dat <- phe[,c(3,i+3)]
  form <- formula(paste0(protList$protein[i], " ~ GRScat"))
  protGLM <- glm(data = dat, formula = form)
  coefs <- data.frame(coef(summary(protGLM)))
  coefs <- coefs %>%
    mutate(GRScat = rownames(coefs),
           protein = protList$protein[i]) %>%
    select(protein, 
           GRScat, 
           estimate = Estimate,
           SE = Std..Error,
           t = t.value,
           P = Pr...t..)
  coefTall <- gather(coefs[,c("estimate","SE","t","P")])
  coefTall$protein <- protList$protein[i]
  coefTall$GRScat <- coefs$GRScat
  coefTall <- coefTall %>%
    select(protein, GRScat, key, value)
  glmCoefsAll <- rbind(glmCoefsAll, coefTall)
  #rm(plotname, bxplot,dat,form,protGLM,coefs,coefTall)
  rm(dat,form,protGLM,coefs,coefTall)
}
sig <- filter(glmCoefsAll, key =="P" & value < 0.05/4034)

# FDR correction. N.B. THIS HARDCODES GRS cats so will be broken for a different trait or threshold. Should be fixed.
Qvals <- glmCoefsAll %>%
  filter(key == "P") %>%
  spread(key = "GRScat", value = "value") %>%
  mutate(`GRScat(3.94,5.3]` = p.adjust(`GRScat(3.94,5.3]`),
         `GRScat(5.3,6.65]` = p.adjust(`GRScat(5.3,6.65]`),
         `GRScat(8.01,9.37]` = p.adjust(`GRScat(8.01,9.37]`),
         `GRScat(9.37,10.7]` = p.adjust(`GRScat(9.37,10.7]`),
         `(Intercept)` = p.adjust(`(Intercept)`)
         )
QvalsT <- gather(Qvals[,3:7],key = "GRScat", value = "value")
QvalsT$protein <- Qvals$protein
QvalsT$key <- "Q"
glmCoefsAll <- bind_rows(glmCoefsAll, QvalsT)
fwrite(glmCoefsAll, file = "glm_coefs_all.csv", sep=",")

# Boxplots
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
  ggsave(paste0("boxplot_",sigProts[i],".png"), plot = bp)
}

# Continuous GRS models?
glmCoefsAllCont <- data.frame()
for( i in 1:nrow(protList)){
  #for( i in 1:10){
  cat(paste0("\nProcessing ",i," of ",nrow(protList),": ",protList$protein[i]))
  # glm
  dat <- phe[,c(2,i+3)]
  form <- formula(paste0(protList$protein[i], " ~ GRS"))
  protGLM <- glm(data = dat, formula = form)
  coefs <- data.frame(coef(summary(protGLM)))
  coefs <- coefs %>%
    mutate(GRS = rownames(coefs),
           protein = protList$protein[i]) %>%
    select(protein, 
           GRS, 
           estimate = Estimate,
           SE = Std..Error,
           t = t.value,
           P = Pr...t..)
  coefTall <- gather(coefs[,c("estimate","SE","t","P")])
  coefTall$protein <- protList$protein[i]
  coefTall$GRS <- coefs$GRS
  coefTall <- coefTall %>%
    select(protein, GRS, key, value)
  glmCoefsAllCont <- rbind(glmCoefsAllCont, coefTall)
  #rm(plotname, bxplot,dat,form,protGLM,coefs,coefTall)
  rm(dat,form,protGLM,coefs,coefTall)
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