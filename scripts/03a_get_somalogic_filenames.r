# Script to get file names of Somalogic proteins in a usable order
# ls -l | awk '{print $9}' > ~/projects/SLE_GRS/soma_list.txt

lapply(c("data.table", "stringr"), require, character.only = T)

protList <- fread("soma_list.txt", data.table = F, h = F)
protNames <- str_split_fixed(protList$V1,"_",3)
protFull <- gsub(".sample","",protNames[,3])
protNames <- str_split_fixed(protNames[,3], "\\.",2)
protList$protein <- protNames[,1]
names(protList)[1] <- "filename"
protList$fullName <- protFull

# Give duplicated names their full identifier
dupes <- protList$protein %in% unique(protList$protein[ duplicated(protList$protein)]) 
protList$protein[dupes] <- protList$fullName[dupes]

write.csv(protList, file = "somalogic_proteins_filenames.csv", quote=F, row.names = F)