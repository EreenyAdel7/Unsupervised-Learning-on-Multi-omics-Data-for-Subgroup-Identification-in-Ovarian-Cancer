############################################################### mRNA Expression  
#reading mrna file
mrna <- read.csv("mrnafpkm.csv",row.names = 1)
View(mrna)
class(mrna)

# Transpose the dataframe to have samples as rows and features as columns
df_transposed <- t(mrna)
class(df_transposed)
df_mrna = as.data.frame(df_transposed)

# Calculate the percentage of zero expression for each feature
zero_percentages <- colMeans(df_mrna == 0)
print(zero_percentages)

# Drop features with zero expression in more than 20% of the samples
features_to_drop <- names(zero_percentages[zero_percentages > 0.2])
print(features_to_drop)
class(features_to_drop)
df_filtered <- df_mrna[ , !(names(df_mrna) %in% features_to_drop)]


# Sort features based on standard deviation
df_std <- apply(df_filtered,2, sd)
df_sorted <- df_filtered[, order(df_std, decreasing = TRUE)]

# Select the top 2000 most variable genes (features)
top_features <- names(df_sorted)[1:2000]
df_top <- df_sorted[, top_features]

# Perform min-max normalization
df_normalized <- apply(df_top, 2, function(x) (x - min(x)) / (max(x) - min(x)))

# Export the normalized data to a CSV file
write.csv(df_normalized, file = "Preprocessed_mRNA.csv")
temp <- df_normalized[, !duplicated(colnames(df_normalized))]
dim(temp)

###################################################### miRNA Expression
#reading mrna file 
mirna <- read.csv("mirna.csv",row.names = 2)
View(mirna)

#Chose only reads_per_million_miRNA_mapped 
df_mirna<- mirna[ , !(names(mirna) %in% c("X"))]
selected_columns <- grep("^reads_per_million_miRNA_mapped_", colnames(df_mirna), value = TRUE)
new_df <- df_mirna[, selected_columns]
colnames(new_df) <- sub("^reads_per_million_miRNA_mapped_", " ", colnames(new_df))
View(new_df)

# Transpose the dataframe to have samples as rows and features as columns
mirna_transposed <- t(new_df)
class(df_transposed)
df_mirna = as.data.frame(mirna_transposed)


# Calculate the percentage of zero expression for each feature
zero_percentages <- colMeans(df_mirna == 0)
print(zero_percentages)

# Drop features with zero expression in more than 20% of the samples
features_to_drop <- names(zero_percentages[zero_percentages > 0.2])
mirna_filtered <- df_mirna[ , !(names(df_mirna) %in% features_to_drop)]
dim(mirna_filtered)

# Perform min-max normalization
mirna_normalized <- apply(mirna_filtered, 2, function(x) (x - min(x)) / (max(x) - min(x)))
View(mirna_normalized[1:5,])
temp <- mirna_normalized[, !duplicated(colnames(mirna_normalized))]
dim(temp)


# Export the normalized data to a CSV file
write.csv(mirna_normalized, file = "Preprocesses_miRNA.csv")









