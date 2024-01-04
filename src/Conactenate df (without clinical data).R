mrna <- read.csv("c_mrna.csv")
mirna <- read.csv("c_mirna.csv")
meth <- read.csv("c_meth.csv")
proteome <- read.csv("c_prot.csv")

#merging the files
mrna_mirna_merged <- merge(mrna, mirna, by = "X")
mrna_mirna_meth_merged <- merge(mrna_mirna_merged, meth, by = "X")
whole_mreged <- merge(mrna_mirna_meth_merged, proteome, by = "X")
dim(whole_mreged)
#saving them as csv 
#set specific column as row names
rownames(whole_mreged) <- whole_mreged$X
whole_mreged$X <- NULL
View(whole_mreged) #296 X 4661
write.csv(whole_mreged, file = "Concatenated_df.csv")

