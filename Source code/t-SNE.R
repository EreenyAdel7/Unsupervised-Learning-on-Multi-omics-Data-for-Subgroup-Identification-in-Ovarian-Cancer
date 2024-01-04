library(Rtsne)
library(ggplot2)
library(RColorBrewer)


# Perform t-SNE on the original data
original_df = read.csv("Multiomics_matrix.csv",row.names = 1)
dim(original_df)
tsne_original <- Rtsne(original_df, dims = 2)

# Create a dataframe for t-SNE results in the original dimension
tsne_original_df <- data.frame(tsne_original$Y)

# Add the labels to the t-SNE results
labels_df = read.csv("Class_Consensus_Clustering_K_10.csv",row.names = 1)
tsne_original_df$label <- labels_df$Clusters

# Define a color palette with a sufficient number of colors for the clusters
num_clusters <- length(unique(tsne_original_df$label))
colors <- brewer.pal(num_clusters, "Paired")

# Plot t-SNE for the original dimension
ggplot(tsne_original_df, aes(x = X1, y = X2, color = label)) +
  geom_point() +
  labs(title = "t-SNE Plot - Original Dimension")

####################################################

# Perform t-SNE on autoencoder reduced data
reduced_df = read.csv("AE_reduced_dataset_1.csv",row.names = 1)

tsne_reduced<- Rtsne(reduced_df, dims = 2)

# Create a dataframe for t-SNE results in the original dimension
tsne_reduced_df <- data.frame(tsne_reduced$Y)

# Add the labels to the t-SNE results
tsne_reduced_df$label <- labels_df$Clusters

# Define a color palette with a sufficient number of colors for the clusters
num_clusters <- length(unique(tsne_reduced_df$label))
colors <- brewer.pal(num_clusters, "Paired")

# Plot t-SNE for the original dimension with different colors for each cluster
ggplot(tsne_reduced_df, aes(x = X1, y = X2, color = factor(label))) +
  geom_point() +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot - AutoEncoder Reduced Dimension")

###################################################
# Perform t-SNE on the PCA reduced data
pca_df = read.csv("PCA_reduced_df.csv",row.names = 1)
View(pca_df)

tsne_pca<- Rtsne(pca_df, dims = 2)

# Create a dataframe for t-SNE results in the original dimension
tsne_pca_df <- data.frame(tsne_pca$Y)

# Add the labels to the t-SNE results
labels_df = read.csv("Class_Consensus_Clustering_K_10.csv",row.names = 1)
tsne_pca_df$label <- labels_df$Clusters

# Define a color palette with a sufficient number of colors for the clusters
num_clusters <- length(unique(tsne_pca_df$label))
colors <- brewer.pal(num_clusters, "Paired")

# Plot t-SNE for the original dimension with different colors for each cluster
ggplot(tsne_pca_df, aes(x = X1, y = X2, color = factor(label))) +
  geom_point() +
  scale_color_manual(values = colors) +
  labs(title = "t-SNE Plot - PCA Reduced Dimension")

