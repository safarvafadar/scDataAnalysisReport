
library(dplyr)
library(tidyverse)
#-----------------------

pbmcSeuratObjectNormalized <- RunPCA(pbmcSeuratObjectNormalized, verbose = FALSE)
ElbowPlot(pbmcSeuratObjectNormalized, ndims = 50) 
pbmcSeuratObjectNormalized <- RunUMAP(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- RunTSNE(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- FindNeighbors(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- FindClusters(pbmcSeuratObjectNormalized, resolution = 0.5)
UMAPPlot(pbmcSeuratObjectNormalized, label =  TRUE) 

DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", split.by = "Condition_study")
#test.use = "bimod" : default :"MAST"
pbmcSeuratObjectNormalized.markers <- FindAllMarkers(pbmcSeuratObjectNormalized, assay = "RNA", test.use = "bimod", only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

View(pbmcSeuratObjectNormalized.markers) 
# Write results to file
write.csv(pbmcSeuratObjectNormalized.markers, "H:/tsk1File/data/results/all_markers.csv", quote = F)
#======================================================
# Merge gene annotations to marker results
pbmcSeuratObjectNormalized.markers <- left_join(pbmcSeuratObjectNormalized.markers, 
                         annotations[, c(1:2, 3, 5)], 
                         by = c("gene" = "gene_id"))

View(pbmcSeuratObjectNormalized.markers)    

#.........................................................
# Return top 10 markers for cluster specified 'x'
gen_marker_table <- function(x){
  pbmcSeuratObjectNormalized.markers[pbmcSeuratObjectNormalized.markers$cluster == x, ] %>%
    head(n=10)}

# Create a data frame of results for clusters 0-6
top10_markers <- map_dfr(0:6, gen_marker_table)

View(top10_markers)
#-----------------------------------
#Option2 
# Identify the 10 most highly variable genes
#OK=>بجای بالایی top10 <- head(VariableFeatures(pbmcSeuratObjectNormalized), 10)
#--------------------------------------------

# Write results to file
write.csv(top10_markers, "H:/tsk1File/data/results/top10_markers.csv", quote = F)
write.csv(pbmcSeuratObjectNormalized.markers, "H:/tsk1File/data/results/pbmcSeuratObjectNormalizedMarkers.csv", quote = F)

#................
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmcSeuratObjectNormalized)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#-------------------
#Assigning cell type identity to clusters

DimPlot(  pbmcSeuratObjectNormalized, reduction =  "tsne", label =  TRUE )



FeaturePlot(object = pbmcSeuratObjectNormalized, 
            features = c("CCL21" , "MMRN1"), 
            reduction = "tsne")
#---------------------------------------------------------------------
#To get a better idea of cell type identity we can explore the expression of different identified markers by cluster using the FeaturePlot() function.
#For example, we can look at the cluster 3 markers by cluster:

FeaturePlot(object = pbmcSeuratObjectNormalized, 
            features = c(top10_markers[top10_markers$cluster == 3, "gene"]), 
            reduction = "tsne")

#---------------------------------------------------------------------
#We can also explore the range in expression of specific markers by using violin plots:

# Vln plot - cluster 3
VlnPlot(object = pbmcSeuratObjectNormalized,         features = c("FBLN5", "GJA5"))
#---------------------------------------------------------------------

#These results and plots can help us determine the identity of these clusters or verify what we hypothesize the identity to be after exploring the canonical markers of expected cell types previously.

#Sometimes the list of markers returned don’t sufficiently separate some of the clusters.
#For instance, we had previously identified clusters 0 and 1 , if we would like to determine the genes that are differentially expressed between these specific clusters, we can use the FindMarkers() function.

# Determine differentiating markers for CD4 T cell clusters 0 versus 1
markers_0vs1 <- FindMarkers(object = pbmcSeuratObjectNormalized, ident.1 = 0, ident.2 = 1)

View(markers_0vs1)

#---------------------------------------------------------------------
# Add gene symbols to the DE table
# markers_0vs1$gene <- rownames(markers_0vs1)
# markers_0vs1 <- left_join(markers_0vs1, 
#                           annotations[, c(1:2, 3, 5)], 
#                           by = c("gene" = "gene_id"))
# 
# View(markers_0vs1)
#---------------------------------------------------------------------

saveRDS(pbmc, file = "H:/tsk1File/data/results/pbmcSeuratObjectNormalized_final.rds")
#---------------------------------------------------------------------

# Visualization
p1 <- DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", group.by = "Condition_study")
p2 <- DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", label = TRUE)
(p1+ p2)
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", split.by = "Condition_study")






for(i in c(0:9)){
  canonicalCellTypeMarkerGenes <- FindConservedMarkers(pbmcSeuratObjectNormalized, ident.1 = i, grouping.var = "Condition_study", verbose = TRUE)
  View(head(canonicalCellTypeMarkerGenes))
  head(canonicalCellTypeMarkerGenes)
  write.csv(head(canonicalCellTypeMarkerGenes),paste0("H:/tsk1File/data/results/canonicalCellTypeMarkerGenesInCluster_", i, ".csv") , quote = F)
}