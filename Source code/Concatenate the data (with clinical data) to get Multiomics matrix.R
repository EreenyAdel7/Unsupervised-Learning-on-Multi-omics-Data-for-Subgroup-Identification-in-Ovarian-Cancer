##this to take only samples shared among the concatenated data (mrna,mirna,proteome,methylation) and clinical 
##this will give me Multiomics data that will be used for further analysis

clinical =  read.delim("ov_tcga_pan_can_atlas_2018_clinical_data.tsv",header = TRUE, sep = "\t")
View(clinical)

clinical_df = clinical[ , !names(clinical) %in% 
      c("Study.ID","Patient.ID")]
View(clinical_df)

clinical_names <- clinical_df$Sample.ID
clinical_names = gsub("-", ".", clinical_names)
clinical_names 
colnames(clinical_df)[1] <- "X"

concatenated_df = read.csv("Concatenated_df.csv")
View(concatenated_df)
c_names <- row.names(concatenated_df)
c_names
c_final <- substr(c_names, 1, 15)
c_final
concatenated_df$X= c_final

cintersect = intersect(c_final, clinical_names)
sum(duplicated(cintersect))
cintersect
length(cintersect) #287


############take only shared samples in clinical data
clinical =  read.delim("ov_tcga_pan_can_atlas_2018_clinical_data.tsv",header = TRUE, sep = "\t")
write.csv(clinical, "CSV_Clinical.csv", row.names=FALSE)
csv_clinical = read.csv("CSV_Clinical.csv",row.names = 3)
clinical_names <- rownames(csv_clinical)
clinical_names = gsub("-", ".", clinical_names)
rownames(csv_clinical) = clinical_names
View(csv_clinical)

preprocessed_clinical <- csv_clinical[cintersect, ]
dim(preprocessed_clinical)
View(preprocessed_clinical)
write.csv(preprocessed_clinical, file = "Preprocessed_clinical.csv")
####################


write.csv(concatenated_df, file = "Concatenated_df_final.csv")
concatenated_df_final <- read.csv("Concatenated_df_final.csv",row.names = 2)
concatenated_df_final = concatenated_df_final[ , !names(concatenated_df_final) %in% 
                          c("X.1")]
View(concatenated_df_final)
selected_rows <- concatenated_df_final[cintersect, ]
View(selected_rows)
dim(selected_rows)
write.csv(selected_rows, file = "Multiomics_Matrix.csv") #287 4661
