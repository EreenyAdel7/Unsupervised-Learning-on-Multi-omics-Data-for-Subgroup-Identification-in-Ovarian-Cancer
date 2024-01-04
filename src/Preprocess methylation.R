#################Concatenate Methylation Data as I have each sample in separate file
# Implement function for that
read_text_documentz <- function(folder_path) {
  files <- list.files(path = folder_path, recursive = TRUE, full.names = TRUE)
  columns <- list()
  
  for (i in seq_along(files)) {
    file <- files[i]
    if (endsWith(file, ".txt")) {
      content <- readLines(file)
      column_name <- paste0("Document", i)
      columns[[column_name]] <- sapply(strsplit(content, "\t"), `[`, 2)
    }
  }
  
  concatenated_df <- as.data.frame(columns)
  return(concatenated_df)
}

# Concatenate it using the function 
folder_path <- "F://Graduation project//GDCdata//TCGA-OV//harmonized//DNA_Methylation//Methylation_Beta_Value"
textd = read_text_documentz(folder_path)
View(textd)

#Put cloumns names correctly
library(readxl)
excel_data <- read_excel("methylation.csv", sheet = 1)
View(excel_data)
methyl_col = excel_data$sorted_cases

colnames(textd) = methyl_col
View(textd[1:5,])
class(textd)


#Put rows names correctly
read_rows <- function(folder_path) {
  files <- list.files(path = folder_path, recursive = TRUE, full.names = TRUE)
  columns <- list()
  
  for (i in seq_along(files)) {
    file <- files[i]
    if (endsWith(file, ".txt")) {
      content <- readLines(file)
      column_name <- paste0("Document", i)
      columns[[column_name]] <- sapply(strsplit(content, "\t"), `[`, 1)
      break 
    }
  }
  
  methyl_rows <- as.data.frame(columns)
  return(methyl_rows)
}


folder_path <- "F://Graduation project//GDCdata//TCGA-OV//harmonized//DNA_Methylation//Methylation_Beta_Value"
rows_names = read_rows(folder_path)
View(rows_names)

methyl_samples = rows_names$Document1
rownames(textd) = methyl_samples
View(textd[1:5,])
df_methylation = textd
View(df_methylation[1:20,])


###################################Preprocessing the methylation data 
# Transpose the dataframe to have samples as rows and features as columns
df_transposed <- t(df_methylation)
df_methyl = as.data.frame(df_transposed)

# Calculate the percentage of zero expression for each feature
zero_percentages <- colMeans(df_methyl == "NA", na.rm = TRUE)
print(zero_percentages)

# Drop features with zero expression in more than 10% of the samples
features_to_drop <- names(zero_percentages[zero_percentages > 0.1])
meth_filtered <- df_methyl[ , !(names(df_methyl) %in% features_to_drop)]
View(meth_filtered[1:20, ])

# to see if there are NA values still 
zero_percentages <- colMeans(meth_filtered == "NA", na.rm = TRUE)
percentage_sum <- sum(zero_percentages)
print(percentage_sum) # 105.1093 to be replaced by KNN 5

# Implement Na values using KNN done on Python as it did not work here 

###################################################################3
# Sort features based on standard deviation
knn_meth = read.csv("KNN_DNA_Methylation.csv",row.names = 1)
View(knn_meth[1:5,])
meth_std <- apply(knn_meth,2, sd)
meth_sorted <- knn_meth[, order(meth_std, decreasing = TRUE)]

# Select the top 2000 most variable genes (features)
top_meth_features <- names(meth_sorted)[1:2000]
meth_top <- meth_sorted[, top_meth_features]
class(meth_top)

# Export the preprocessed data to a CSV file
write.csv(meth_top, file = "Preprocessed_methylation.csv")
temp_meth <- meth_top[, !duplicated(colnames(meth_top))]
dim(temp_meth)
