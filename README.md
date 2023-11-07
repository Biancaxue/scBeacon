# scBeacon Project Overview
![alt text](https://github.com/Biancaxue/scBeacon/blob/main/project_pipeline.jpg?raw=true)

# Introduction
The cellular components of tumors and their microenvironment play pivotal roles in tumor progression, patient survival, and the response to cancer treatments. Our contribution, scBeacon, is a novel tool that derives cell type signatures by integrating and clustering multiple scRNA-seq datasets to extract signatures for deconvolving unrelated tumor datasets on bulk samples. By employing scBeacon on the TCGA cohort, we find previously unrecognized cellular and molecular attributes within specific tumor categories, many with patient outcome relevance. We developed a tumor cell-type map (TCT) to visually depict the relationships among TCGA samples based on the cell-type inferences.

# Code description
EBI_datasets_processing.R: The code used to process the 62 EMBL-EBI scRNA-seq datasets. We generated clusters using Seurat with default settings. If the dataset is too big, we used Geometric Sketching to subset the datasets.

scBeacon.R: The script that integrates the clusters generated from EBI_datasets_processing.R, then computes ranked centroids for the signature matrix used in deconvolution.

multimodel_analysis.py: In this script, we fit multi-model on TCGA CIBERSORT results. Multi-model fitting will be applied per signature and per tumor type. We only keep the ones that have bi-modality and save the parameters.

survival_analysis.R: run multi-variant survival analysis based on CIBERSORT results
