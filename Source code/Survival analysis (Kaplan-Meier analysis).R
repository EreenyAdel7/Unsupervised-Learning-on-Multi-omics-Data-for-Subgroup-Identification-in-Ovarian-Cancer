#Kaplan-Meier analysis for seeing the survival analysis between clusters resulted from clustering for both AE reduced data and PCA reduced data

#=====================================================================================
# IP files
ip_clinical_path <- "Preprocessed_clinical.csv" 
ip_cluster_path <- "Class_Consensus_Clustering_K_10.csv" ##done for both "Class_Consensus_Clustering_K_10.csv" get from Autoencoder and PCA
ip_cluster_path 

# Read clinical data
clinical_data <- fread(ip_clinical_path, header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
clinical_data[1:5, ]
View(clinical_data)
colnames(clinical_data)[1] <- "Samples"
View(clinical_data)
colnames(clinical_data)

# Read cluster labels
cons_res <- fread(ip_cluster_path, header = TRUE, data.table = FALSE, stringsAsFactors = FALSE)
cons_res[1:5, ]
View(cons_res)
colnames(cons_res)[1] <- "Samples"
View(cons_res)

####this is only apply for PCA
cons_res$Samples = clinical_data$Samples
View(cons_res)

# Combine clinical data and cluster results
Clinical_cluster <- merge(cons_res, clinical_data, by = "Samples")
View(Clinical_cluster)
#=====================================================================================
# For overall survival
sum(is.na(Clinical_cluster$Overall.Survival.Status))
sum(is.na(Clinical_cluster$Overall.Survival..Months.))
Clinical_cluster$Overall.Survival..Months.[is.na(Clinical_cluster$Overall.Survival..Months.)] <- 0 #replace the na with zero as it should be numeric 
sum(is.na(Clinical_cluster$Overall.Survival..Months.))
class(Clinical_cluster$Overall.Survival.Status) #character --> should be numeric for creating survival object
os_status = Clinical_cluster$Overall.Survival.Status
os_status <- substr(os_status, 1, 1)
os_status 
os_status_num = as.numeric(os_status)
Clinical_cluster$Overall.Survival.Status = os_status_num

# Creating Survival Object
Surv_Obj <- Surv(Clinical_cluster$Overall.Survival..Months.,
                 Clinical_cluster$Overall.Survival.Status)

surv_diff <- survdiff(formula = Surv_Obj ~ Clinical_cluster$Clusters)
p_val <- (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))

fit <- do.call(survfit, list(formula = Surv_Obj ~ Clinical_cluster$Clusters, data = Clinical_cluster))

title = "Overall survival"
survp_os <- ggsurvplot(fit, legend = c(0.8,0.8), pval = TRUE, title = title,
                       surv.median.line="hv", 
                       surv.scale = "percent", ylab = "Surviving",  
                       legend.labs = levels(Clinical_cluster$Clusters), 
                       xlab = "Survival time (Months)",  
                       legend.title = "", font.legend = 12, palette = c("#ac92eb", "#4fc1e8", "#a0d568", "#ffce54", "#ed5564",
                                                                        "#f8a5c2", "#6dcaff", "#00bfbf", "#ffb878", "#8d6b94"))

png_op_os_path <- "OverallSurvival_PCA(10).png"
png(png_op_os_path)
print(survp_os$plot, newpage = FALSE)
dev.off()

#======================================================================================
# For disease free survival - DFS
if (!requireNamespace("imputeTS", quietly = TRUE)) {
  install.packages("imputeTS")
}
library(imputeTS)
sum(is.na(Clinical_cluster$Disease.Free..Months)) #149
sum(is.na(Clinical_cluster$Disease.Free.Status)) #149

# Create a copy of the variable
dfs_status <- Clinical_cluster$Disease.Free.Status
dfs_status <- substr(dfs_status, 1, 1)
dfs_status_num = as.numeric(dfs_status)
# Identify missing values
missing_values <- is.na(dfs_status_num)
kmeans_data <- kmeans(dfs_status_num[!missing_values], centers = 2)
dfs_status_num[missing_values] <- kmeans_data$centers[kmeans_data$cluster]
# Replace the original variable with the imputed values
Clinical_cluster$Disease.Free.Status <- dfs_status_num




# Create a copy of the variable
dfm_status <- Clinical_cluster$Disease.Free..Months
# Identify missing values
missing_values_2 <- is.na(dfm_status)
kmeans_data2 <- kmeans(dfm_status[!missing_values_2], centers = 5)
dfm_status[missing_values_2] <- kmeans_data2$centers[kmeans_data2$cluster]
# Replace the original variable with the imputed values
Clinical_cluster$Disease.Free..Months <- dfm_status


# Creating Survival Object
Surv_Obj_DFS <- Surv(Clinical_cluster$Disease.Free..Months,
                     Clinical_cluster$Disease.Free.Status)

surv_diff <- survdiff(formula = Surv_Obj_DFS ~ Clinical_cluster$Clusters)
p_val <- (1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1))

fit <- do.call(survfit, list(formula = Surv_Obj_DFS ~ Clinical_cluster$Clusters, data = Clinical_cluster))

title = "Disease free survival"
survp_dfs <- ggsurvplot(fit, legend = c(0.8,0.8), pval = TRUE, title = title,
                        surv.median.line="hv", 
                        surv.scale = "percent", ylab = "Surviving",  
                        legend.labs = levels(Clinical_cluster$Clusters), 
                        xlab = "Survival time (Months)",  
                        legend.title = "", font.legend = 12, palette = c("#ac92eb", "#4fc1e8", "#a0d568", "#ffce54", "#ed5564",
                                                                         "#f8a5c2", "#6dcaff", "#00bfbf", "#ffb878", "#8d6b94"))

png_op_dfs_path <- "Disease Free Survival.png"
png(png_op_dfs_path)
print(survp_dfs$plot, newpage = FALSE)
dev.off()
