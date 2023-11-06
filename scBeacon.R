#Code review for scBeacon pipeline

library(Matrix)
library(biomaRt)
library(viper)
library(igraph)
library(survminer)
library(survival)
library(ggplot2)

# set working directory to where datasets are stored
setwd(data_dir)

# generate full transcriptome ensambl_gene_id, entrezgene ud, gene name mapping table
# this table will be used as a gene id ~ name converter and force all datasets to have the same gene profile
gene_id_converter <- function(species="hsapiens"){
  mart <- useDataset(paste(species, "_gene_ensembl", sep = ""), useMart("ensembl")) #choose hsapiens_gene_ensembl from bioMart database
  table <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name"), mart=mart) #retrieves information from the BioMart database
  table <- table[apply(table, 1, function(x) sum(is.na(x)) == 0), ] #remove NAs in the table
  table <- table[!duplicated(table$external_gene_name), ] #remove duplicates
  return(table)}
table <- gene_id_converter(species="hsapiens") 
table <- table[-which(table$external_gene_name == ""),] # remove genes that don't have gene names

# generate ranked centroid and matrix(optional) from clusters                       
generate_centroids <- function(files, species = "hsapiens"){
  for (f in files){
    message(f)
    index <- unlist(strsplit(f, split="/"))[2] #index is the clusters folder, either kmean or louvain
    exp <- as.matrix(readMM(f)) #read in gene expression matrix
    genes <- read.delim(gsub("matrix.mtx", "genes.tsv", f), header = F, stringsAsFactors = F)
    cells <- read.delim(gsub("matrix.mtx", "cells.tsv", f), header = F, stringsAsFactors = F)
    rownames(exp) <- genes$V1
    colnames(exp) <- cells$V1
    id <- which.max(apply(table, 2, function(x, y) sum(x %in% y), y = rownames(exp))) #id represents the index of gene annotation in table ("ensembl_gene_id", "entrezgene_id", "external_gene_name")
    exp <- exp[rownames(exp) %in% table[, id], ]
    rownames(exp) <- table$external_gene_name[match(rownames(exp), table[, id])] #convert other gene annotations to external_gene_name
    centroid <- structure(matrix(0, nrow(table), ncol(exp)), dimnames = list(table$external_gene_name, colnames(exp))) 
    centroid[rownames(exp), ] <- exp
    centroid <- apply(centroid, 2, rank)/nrow(centroid) #for each sample, rank gene expression profile
    metadata <- sapply(strsplit(f, split = "_"), function(x) paste(x[1:(length(x) - 1)], collapse = "_")) #remove "_matrix.mtx" in file names
    metadata <- sapply(strsplit(metadata, split = "/"), function(x) x[length(x)]) #remove directory and "/" in file names
    #generate ranked matrix
    setwd(paste(data_dir, dataset, paste(index, "ranked_matrix", sep="_"), sep="/"))
    write.table(centroid, file = paste(metadata, "_matrix.csv", sep = ""), sep=",", row.names=TRUE, col.names=TRUE)
    #generate centroids
    setwd(paste(data_dir, dataset, paste(index, "centroids", sep="_"), sep="/"))
    centroid <- matrix(rowMeans(centroid), nrow(centroid), 1, dimnames = list(rownames(centroid), metadata))
    write.table(centroid, file=paste(metadata, "_centroid.csv", sep = ""), sep=",", row.names=TRUE, col.names=TRUE)
    setwd(paste(data_dir, dataset, sep="/"))
  }}

# list the directories of datasets, and run generate_centroids function for each dataset
dirs <- list.dirs(path=".",full.names=TRUE, recursive=FALSE) 
for (dir in dirs){
  message(dir)
  setwd(paste(data_dir, dir, sep="/"))
  unlink("./louvain_centroids", recursive=TRUE) #delete the folder if it already exists
  unlink("./louvain_ranked_matrix", recursive=TRUE)
  dataset <- tail(unlist(strsplit(getwd(), split="/")), 1)
  dir.create(paste(data_dir, dataset, "louvain_ranked_matrix", sep="/"))
  dir.create(paste(data_dir, dataset, "louvain_centroids", sep="/"))
  louvain_clusters <- list.files(path="./louvain", pattern="matrix.mtx", full.names = TRUE) #list matrix.mtx files in louvain clustering folder
  louvain_cell_number <- as.numeric(unlist(lapply(louvain_clusters, function(x) unlist(strsplit(x, split="_"))[8]))) 
  valid_louvain_clusters <- louvain_clusters[which(louvain_cell_number > 50)] #filter down clusters that have cell number lower than 50
  generate_centroids(valid_louvain_clusters, species = "hsapiens")
}

# load centroids to a matrix
load_centroids <- function(files){
  cell.numbers <- as.numeric(unlist(lapply(files, function(f) tail(unlist(strsplit(f, split="_")), 2)[1])))
  files <- files[which(cell.numbers >= 50)]
  centroids <- lapply(files, function(f) read.csv(f))
  centroids <- do.call(cbind, centroids)
  colnames(centroids) <- gsub("_centroid.csv", "", files)
  return(as.matrix(centroids))
}

files <- list.files(path=data_dir, pattern="centroid.csv", recursive=TRUE, full.names = T)
centroids <- load_centroids(files) #centroids is a matrix of gene by cluster centroids

# compute viper pairwise similarity score between centroids
get_similarity <- function(centroids, nn = floor(nrow(centroids) * 0.1))
  return(viperSimilarity(centroids, nn = nn, method = "greater")) 

# compute an empirical cumulative distribution function, then find a threshold of similarity score using quantile of a given probability p
get_threshold <- function(similarity, p = 0.01, method = c("greater", "less")){
  switch(method, greater = {
    threshold <- quantile(ecdf(similarity[upper.tri(similarity)]), 1 - p)}, less = {
      threshold <- quantile(ecdf(similarity[upper.tri(similarity)]), p)})
  return(threshold)}

# compute adjacency matrix, only connect centroids that have similarity score above the threshold, then use louvain clustering to find clusters
cluster_centroids <- function(centroids, p=0.01, nn=floor(nrow(centroids)*0.1), verbose=F){
  similarity <- get_similarity(centroids, nn=nn)
  threshold <- get_threshold(similarity, p=p, method="greater")
  g <- graph_from_adjacency_matrix(as.matrix(similarity > threshold), mode="undirected", diag=F)
  cluster <- cluster_louvain(g)
  return(cluster)}

# meta-clusters information is included in the variable clusters
clusters <- cluster_centroids(centroids, p=0.01, nn=floor(nrow(centroids)*0.1), verbose=F)

# function for computing signature matrix from meta-clusters
generate_signature_matrix <- function(clusters){
  #compute meta-cluster gene expression profile by taking the average of centroids
  clusters.mean <- list()
  names(clusters$membership) <- clusters$names
  for (cluster in unique(clusters$membership)){
    # if there are more than one centroid in a meta-cluster, compute the average expression
    if (length(names(which(clusters$membership==cluster))) > 1){
      exp.mean <- apply(centroids[, names(which(clusters$membership==cluster))], 1, mean)
      clusters.mean[[cluster]] <- exp.mean
    } else { 
      # if there is only one centroid in a meta-cluster, use that centroid
      exp.mean <- centroids[, names(which(clusters$membership==cluster))]
      clusters.mean[[cluster]] <- exp.mean
    }
  }
  clusters.mean <- do.call(cbind, clusters.mean)
  colnames(clusters.mean) <- unique(clusters$membership)
  #differential expression genes
  #top 20% most differential expressed genes from top1 and top2 cluster difference
  diff.exp.list <- list()
  for (gene in rownames(clusters.mean)){
    highest.cluster <- names(sort(clusters.mean[gene, ], decreasing=TRUE))[1]
    second.highest.cluster <- names(sort(clusters.mean[gene, ], decreasing=TRUE))[2]
    diff.exp <- clusters.mean[gene, highest.cluster] - clusters.mean[gene, second.highest.cluster]
    diff.exp.list[[gene]] <- diff.exp
  }
  diff.exp.list <- unlist(diff.exp.list)
  n <- 20
  diff.exp.genes <- names(diff.exp.list[diff.exp.list > quantile(diff.exp.list,prob=1-n/100)])
  sig <- clusters.mean[diff.exp.genes,]
  sig
}
sig <- generate_signature_matrix(clusters)
write.table(sig, file=paste0("signature_matrix_", ncol(sig),".tsv"), sep="\t")

# plot the signature matrix
library(pheatmap)
plot_df <-sig[rownames(sig)[order(apply(sig, 1, which.max))],]
colnames(plot_df) <- unlist(lapply(colnames(plot_df), function(x) paste("X", x, sep="")))
pheatmap(t(plot_df), scale="column", show_rownames=T, show_colnames=F, cluster_rows=F, cluster_cols=F, main="signature matrix")
dev.off()





