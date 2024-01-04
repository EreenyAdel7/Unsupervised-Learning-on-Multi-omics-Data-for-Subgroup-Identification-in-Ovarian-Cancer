## Unsupervised-Learning-on-Multi-omics-Data-for-Subgroup-Identification-in-Ovarian-Cancer
# Project Overview
  In this project, an autoencoder-based approach is followed to non-linearly project high-dimensional multi-omics (mRNA, miRNA, methylation, and protein expression) ovarian cancer (OC) data to a lower-dimensional space. Also, principle component analysis is used for producing a lower-dimensional space to validate the accuracy of the autoencoder by comparing results obtained from analysis on both. The compressed data is then subjected to consensus k-means clustering to identify the clusters. Survival analysis of the resulting clusters revealed an insignificant difference in overall survival and disease-free survival.

# Installation Instructions
- Check this repository.
- For R code parts, install the required dependencies using the following command: pip install -r requirements.txt.

# Usage
- In R, set the working directory of the folder that has the downloaded datasets.
- Dimensionality reduction.
- Consensus k-means clustering.
- Survival analysis on overall survival and disease-free survival 

# Additional information 
- For Python code parts, use scikit-learn and TensorFlow in Colab directly.
- For more instructions, see the code and the documentation provided in this repository.
