mrna <- read.csv("Preprocessed_mRNA.csv",row.names = 1)
mirna <- read.csv("Preprocessed_miRNA.csv",row.names = 1)
meth <- read.csv("Preprocessed_DNA_methylation.csv",row.names = 1)
proteome <- read.csv("Preprocessed_Proteome_Profiling.csv",row.names = 1)
######################################extract names 
index_names1 <- row.names(mrna)
index_names2 <- row.names(mirna)
index_names3 <- row.names(meth)
index_names4 <- row.names(proteome)

##########################################Intersection between mrna and mirna (there was space in mirna)
index_names2 <- trimws(index_names2, "left")
mrna_mirna = intersect(index_names1, index_names2)
print(mrna_mirna)            

#######################################there was - instead of . in meth - change it and save new csv file
index_names3 <- gsub("-", ".", index_names3)
index_names4 <- gsub("-", ".", index_names4)

#################################################
#Extract the sample IDs from each indicies before fouth Dot due to different formatting 
sample_ids_csv1 <- substr(index_names1, 1, 16)
sample_ids_csv2 <- substr(index_names2, 1, 16)
sample_ids_csv3 <- substr(index_names3, 1, 16)
sample_ids_csv4 <- substr(index_names4, 1, 16)

mrna_mirna_intersect = intersect(sample_ids_csv1, sample_ids_csv2) #no duplicated samples ####426 samples
mrna_mirna_meth_intersect <- intersect(mrna_mirna_intersect, sample_ids_csv3) ###416 samples
common_intersect  <- intersect(mrna_mirna_meth_intersect, sample_ids_csv4)  ###296 samples along all of them 

############################################################################
#removed duplicated samples from miRNA (2 samples)
for (r in rownames(mirna)){
  if (startsWith(r, "TCGA.09.0366.01A")){
    print(r)
  }
} #"TCGA.09.0366.01A.01R.1986.13" , "TCGA.09.0366.01A.01R.1564.13"

row_names_df_to_remove = c("TCGA.09.0366.01A.01R.1986.13" , "TCGA.09.0366.01A.01R.1564.13")
mirna_f = mirna[!(row.names(mirna) %in% row_names_df_to_remove),]
dim(mirna_f) #[1] 497 437
df_mirna = as.data.frame(mirna_f)
View(df_mirna)


length(sample_ids_csv2)
temp = sample_ids_csv2
temp_unique = unique(temp)
length(temp_unique)
final_mirna = temp_unique[ !temp_unique == 'TCGA.09.0366.01A']
length(final_mirna)
sum(duplicated((final_mirna)))

##############################################################
#Uniform the format of row names (samples) of all csv files 
rownames(mrna) = sample_ids_csv1
rownames(df_mirna) = final_mirna #‘TCGA.09.0366.01A’ is existed twice 
rownames(meth) = sample_ids_csv3
rownames(proteome) = sample_ids_csv4

##############################################################
#make the indicies of each file as first column
write.csv(mrna, file = "c_mrna.csv")
write.csv(df_mirna, file = "c_mirna.csv")
write.csv(meth, file = "c_meth.csv")
write.csv(proteome, file = "c_prot.csv")


