#this script is for preprocessing datasets downloaded directly from EBI websites, the output files are input for scBeacon 
library(Seurat)
library(Matrix)
library(biomaRt)
library(igraph)
library(reticulate)
library(rsvd)
geosketch <- import('geosketch')


setwd(datasets_dir) 
path <- getwd()
dirs <- list.dirs(path=".",full.names=TRUE, recursive=FALSE)

#down sampling big datasets using geometric sketch, loaded from python code (https://github.com/brianhie/geosketch)
#create louvain clusters using Seurat package
for (dir in dirs){
  setwd(dir)
  message("now processing:  ", dir)
  f.umi <- list.files(path = ".", pattern = "filtered_counts.*mtx$",full.names = TRUE)
  f.gene <- list.files(path = ".", pattern = "filtered_counts.mtx_rows$",full.names = TRUE)
  f.barcode <- list.files(path = ".", pattern = "filtered_counts.mtx_cols$",full.names = TRUE)
  umi <- readMM(file = f.umi)   #read in gene expression matrix
  gene.names <- read.delim(f.gene, header = FALSE, stringsAsFactors = FALSE)
  barcode.names <- read.delim(f.barcode, header = FALSE, stringsAsFactors = FALSE)
  colnames(umi) <- barcode.names$V1
  rownames(umi) <- gene.names$V1
  umi <- as.matrix(umi)
  message("number of cells in original dataset: ", ncol(umi))
  # remove duplicated genes
  all_genes <- rownames(umi)
  duplicate_genes <- all_genes[duplicated(all_genes)]
  if (length(duplicate_genes) > 0){
    umi <- umi[-which(all_genes %in% duplicate_genes), ]
  }
  if (ncol(umi) >= 100000){
    message("too many cells, now use geosketch")
    umi <- t(umi)
    # Get top PCs from randomized SVD.
    s <- rsvd(umi, k=10)
    X.pcs <- s$u %*% diag(s$d)
    # Sketch 10000 cells
    sketch.size <- as.integer(10000) 
    sketch.indices <- geosketch$gs(X.pcs, sketch.size)
    umi <- umi[unlist(sketch.indices)+1, ]
    umi <- t(umi)
  }
  message("now nubmer of cells after geosketch: ", ncol(umi))
  # metadata
  Lab <- substr(unlist(strsplit(f.umi, split=".", fixed=TRUE))[2], 2, 100)
  Platform <- "Platform"
  Species <- "hsapiens"
  DevelopmentalStage <- "DevelopmentalStage"
  Batch <- 1
  Tissue <- "Tissue"
  # Louvain clustering using Seurat
  tpm <- CreateSeuratObject(counts = umi) 
  if (sum(GetAssayData(object = tpm)) %%1 == 0){
    tpm <- NormalizeData(tpm, normalization.method = "LogNormalize", scale.factor = 1000000, verbose=FALSE)
  }
  #deteremine how many PCs to keep 
  genes <- rownames(tpm)
  tpm <- ScaleData(tpm, features = genes, verbose=FALSE, do.scale=FALSE, do.center=FALSE)
  tpm <- FindVariableFeatures(tpm, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
  tpm <- RunPCA(tpm, verbose=FALSE)
  tpm <- FindNeighbors(tpm, dims = 1:10, verbose=FALSE)
  tpm <- FindClusters(tpm, resolution = 0.5, verbose=FALSE)
  cell_anno <- Idents(tpm) #output
  n_cluster <- length(levels(cell_anno))
  message("number of clusters in ", Lab, " : ", n_cluster)
  for (i in levels(cell_anno)){
    cell_id <- cell_anno[cell_anno == i]
    CellCluster <- i
    cluster_mat <- umi[, names(cell_id)]
    N <- length(cell_id)
    cluster_dgtmat <- Matrix(cluster_mat, sparse = TRUE)
    filename <- paste(paste(c(Lab, Platform, Species, DevelopmentalStage, Batch, Tissue, CellCluster, N), collapse = "_"), "_matrix.mtx", sep="")
    writeMM(cluster_dgtmat,file=filename)
    write.table(colnames(cluster_dgtmat), file=gsub("matrix.mtx", "cells.tsv", filename), row.names=FALSE, col.names=FALSE)
    write.table(rownames(cluster_dgtmat), file=gsub("matrix.mtx", "genes.tsv", filename), row.names=FALSE, col.names=FALSE)
    message("processed: ", filename)
  }
}
