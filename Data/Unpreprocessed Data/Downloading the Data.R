## mRNA - miRNA - DNA Methylation - Proteome profiling data is downloaded using TCGABiolinks R package.
##here is the instructions and code for downloaing it 


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("writexl")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("maftools")
BiocManager::install("sesame")
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(pheatmap)
library(SummarizedExperiment)
library(sesameData)
library(sesame)


gdcprojects <- getGDCprojects()
getProjectSummary('TCGA-OV')



#########RNA
query_TCGA_1<- GDCquery(project = 'TCGA-OV',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        access = 'open')
less_query_TCGA <- getResults(query_TCGA_1)
GDCdownload(query_TCGA_1)
mrna_prep = GDCprepare(query_TCGA_1,summarizedExperiment = TRUE)
class(mrna_prep)

mrna_matrix = assay(mrna_prep, "fpkm_unstrand")
View(mrna_matrix)
write.csv(mrna_matrix, file = "F:\\Graduation project\\mrnafpkm.csv")

#####miRNA
query_mirna<- GDCquery(project = 'TCGA-OV',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'miRNA-Seq',
                        access = 'open',
                       data.type = 'miRNA Expression Quantification')
miRNA_query_TCGA <- getResults(query_mirna)
GDCdownload(query_mirna)
mirna_prep = GDCprepare(query_mirna, summarizedExperiment = TRUE)

mirna_matrix <- as.matrix(mirna_prep)
View(mirna_prep)
rownames(mirna_matrix) <- rownames(mirna_prep)
colnames(mirna_matrix) <- colnames(mirna_prep)
 # Add metadata if available
mirna_se <- SummarizedExperiment(assays = list(counts = mirna_matrix))


mirna_matrixs = assay(mirna_se)
write.csv(mirna_matrix, file = "F:\\Graduation project\\mirna2.csv")


####Proteome Profiling
query_TCGA_4<- GDCquery(project = 'TCGA-OV',
                        data.category = 'Proteome Profiling',
                        access = 'open',
                        data.type = "Protein Expression Quantification")
proteome_query_TCGA <- getResults(query_TCGA_4)
GDCdownload(query_TCGA_4)
prot_prep = GDCprepare(query_TCGA_4, summarizedExperiment = TRUE)


prot_matrix <- as.matrix(prot_prep)
rownames(prot_matrix) <- rownames(prot_prep)
colnames(prot_matrix) <- colnames(prot_prep)
# Add metadata if available
prot_se <- SummarizedExperiment(assays = list(counts = prot_matrix))
prot_matrix = assay(prot_se)

write.csv(prot_matrix, file = "F:\\Graduation project\\proteomeprofiling.csv")


###DNA Methylation
query_methly <- GDCquery(project = 'TCGA-OV',
                         data.category = 'DNA Methylation',
                         platform = 'Illumina Human Methylation 27',
                         access = 'open',
                         data.type = 'Methylation Beta Value')

output_query_methyl <- getResults(query_methly)

GDCdownload(query_methly)
View(output_query_methyl["cases"])

library(dplyr)

# Assuming your dataframe is named 'df' with columns 'ids' and 'cases'
# Sorting the dataframe based on the 'ids' column
df_sorted <- output_query_methyl %>% arrange(id)

# Fetching the 'cases' column corresponding to the sorted 'ids'
sorted_cases <- df_sorted$cases
sorted_cases <- c("", sorted_cases)
library(writexl)
write_xlsx(data.frame(sorted_cases), "F:\\Graduation project\\methylation.csv")

# Printing the sorted cases
View(sorted_cases)


dna.meth <- GDCprepare(query_methly, summarizedExperiment = TRUE)
assay(dna.meth)  

write.csv(meth_matrix, file = "F:\\Graduation project\\methylation.csv")

save.image(file ="GraduationProject")


