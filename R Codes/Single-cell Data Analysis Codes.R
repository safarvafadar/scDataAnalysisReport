library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)

library(stringr)

library(dplyr)
library(tidyverse)
library(cowplot)

library('BiocManager')
library('multtest')
library('metap')
#---------------------------------------------------
metadata=readRDS("H:/tsk1File/data/esophagus_ID16538_metadataMatrix.rds")
pbmc <- readRDS(file = "H:/tsk1File/data/esophagus_ID16021_dataMatrixNormalized.rds")
#---------------------------------------------------
pbmc_Copy <- data.frame(pbmc)            # Create copy of data
rownames(pbmc_Copy)=pbmc_Copy[,1]
remove(pbmc)
#The data frame should be exactly same as column names, and genes as row names. // replace ".." with ": "
colnames(pbmc_Copy) = gsub("\\.\\.", ": ", colnames(pbmc_Copy), fixed=FALSE) 


metadata_Copy <- data.frame(metadata)            # Create copy of metadata
rownames(metadata_Copy)=metadata_Copy[,1]

pbmc_Copy =pbmc_Copy[,-1]
#---------------------------------------------------


pbmcSeuratObject  <- CreateSeuratObject(counts = pbmc_Copy , meta.data =  metadata_Copy)# , min.cells = 0  
#---------------------------------------------------

plotBefore=FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(pbmcSeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#---------------------------------------------------  


plot1 <- FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#--------------------------------------------------------------

metadata %>% ggplot(aes(x=Condition_study, fill=Condition_study)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
#---------------------------------------------------  
#number of counts.png Figure
metadata %>% 
  ggplot(aes(color=Condition_study, x=nbr_count, fill= Condition_study)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Number of Counts (log 10)") +
  geom_vline(xintercept = 500)

#--------------------------------------------------- 
# Number of Genes.png
metadata %>% 
  ggplot(aes(color=Condition_study, x=nbr_features, fill= Condition_study)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Number of Genes (log 10)") +
  geom_vline(xintercept = 500)
#-----------------------------------------------------
# Correlation between genes detected and number of UMIs.png
metadata %>% 
  ggplot(aes(x=log10(nbr_count), y=log10(nbr_features), color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "blue", high = "green") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = log10(500), color="red", linetype=2) +
  geom_vline(xintercept = log10(10000), color="red", linetype=2) +
  geom_hline(yintercept = log10(250), color="red", linetype=2) +
  geom_hline(yintercept = log10(10000), color="red", linetype=2) +
  facet_wrap(~Condition_study)


#--------------------------------------------------------------------
#Normalizing the data
#--------------------------------------------------------------------  
pbmcSeuratObjectNormalized <- NormalizeData(pbmcSeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

plotAfter=FeatureScatter(pbmcSeuratObjectNormalized, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotBefore+plotAfter
VlnPlot(pbmcSeuratObjectNormalized, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#---------------------------------------------------  
#Identification of highly variable features (feature selection)
#---------------------------------------------------  
pbmcSeuratObjectNormalized <- FindVariableFeatures(pbmcSeuratObjectNormalized, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmcSeuratObjectNormalized), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmcSeuratObjectNormalized)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#---------------------------------------------------  
#Scaling the data
#---------------------------------------------------  
allGenes <- rownames(pbmcSeuratObjectNormalized)
pbmcSeuratObjectNormalized <- ScaleData(pbmcSeuratObjectNormalized, features = allGenes)

#---------------------------------------------------  
#linear dimensional reduction
#---------------------------------------------------
pbmcSeuratObjectNormalized <- RunPCA(pbmcSeuratObjectNormalized, features = VariableFeatures(object = pbmcSeuratObjectNormalized))

VizDimLoadings(pbmcSeuratObjectNormalized, dims = 1:2, reduction = "pca")

#---------------------------------------------------  
#Q2
#---------------------------------------------------  

pbmcSeuratObjectNormalized <- RunPCA(pbmcSeuratObjectNormalized, verbose = FALSE)
ElbowPlot(pbmcSeuratObjectNormalized, ndims = 50) 
pbmcSeuratObjectNormalized <- RunUMAP(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- RunTSNE(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- FindNeighbors(pbmcSeuratObjectNormalized, reduction = "pca", dims = 1:16)
pbmcSeuratObjectNormalized <- FindClusters(pbmcSeuratObjectNormalized, resolution = 0.5)
UMAPPlot(pbmcSeuratObjectNormalized, label =  TRUE) 

DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", split.by = "Condition_study", label =  TRUE)


pbmcSeuratObjectNormalized.allMarkers <- FindAllMarkers(pbmcSeuratObjectNormalized, assay = "RNA", test.use = "bimod", only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.25)
pbmcSeuratObjectNormalized.onlyPositiveMarkers <- FindAllMarkers(pbmcSeuratObjectNormalized, assay = "RNA", test.use = "bimod", only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

dim(pbmcSeuratObjectNormalized.allMarkers) 
dim(pbmcSeuratObjectNormalized.onlyPositiveMarkers) 
# Write results to file
write.csv(pbmcSeuratObjectNormalized.allMarkers, "H:/tsk1File/data/results/all_markers.csv", quote = F)
write.csv(pbmcSeuratObjectNormalized.onlyPositiveMarkers, "H:/tsk1File/data/results/onlyPositiveMarkers.csv", quote = F)
#======================================================

#.........................................................
# Return top 10 markers for cluster specified 'x'
gen_marker_table <- function(x){
  pbmcSeuratObjectNormalized.markers[pbmcSeuratObjectNormalized.markers$cluster == x, ] %>%
    head(n=10)}

# Create a data frame of results for clusters 0-9
top10_markers <- map_dfr(0:9, gen_marker_table)

#--------------------------------------------

# Write results to file
write.csv(top10_markers, "H:/tsk1File/data/results/top10_markers.csv", quote = F)
write.csv(pbmcSeuratObjectNormalized.markers, "H:/tsk1File/data/results/pbmcSeuratObjectNormalizedMarkers.csv", quote = F)

#-------------------

DimPlot(  pbmcSeuratObjectNormalized, reduction =  "tsne", label =  TRUE,split.by = "Condition_study" )



FeaturePlot(object = pbmcSeuratObjectNormalized, 
            features = c("CCL21" , "MMRN1"), split.by = "Condition_study", label =  TRUE,
            reduction = "tsne")

#---------------------------------------------------------------------

# Visualization
p1 <- DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", group.by = "Condition_study", label = TRUE)
p2 <- DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", label = TRUE)
p1+ p2
#To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(pbmcSeuratObjectNormalized, reduction = "umap", split.by = "Condition_study")




for(i in c(0:9)){
  canonicalCellTypeMarkerGenes <- FindConservedMarkers(pbmcSeuratObjectNormalized, ident.1 = i, grouping.var = "Condition_study", verbose = TRUE)
  head(canonicalCellTypeMarkerGenes)
  write.csv(head(canonicalCellTypeMarkerGenes),paste0("H:/tsk1File/data/results/canonicalCellTypeMarkerGenesInCluster_", i, ".csv") , quote = F)
}


#---------------------------------------------------
# test for doc (Cluster 0)
# for cluster 0  (ConservedMarkers)
FeaturePlot(pbmcSeuratObjectNormalized, features = c("MMP2","COL4A2","RND1","IGFBP4","SELE","CD74"), min.cutoff = "q9") 

FeaturePlot(pbmcSeuratObjectNormalized, features = c("TFF3","CCL21","LYVE1","EFEMP1","MMRN1","SCNN1B"), min.cutoff = "q9") 

#---------------------------------------------------
# test for doc (Cluster 5, good avg_log2FC and p_val_adj)
canonicalCe5 <- FindConservedMarkers(pbmcSeuratObjectNormalized, ident.1 = 5, grouping.var = "Condition_study", verbose = FALSE)
View(head(canonicalCe5))

#---------------------------------------------------
#Plot for "CCL21" and "MMRN1" in ten clusters for Q4
FeaturePlot(pbmcSeuratObjectNormalized, features = c("CCL21","MMRN1"), min.cutoff = "q9" , split.by = "Condition_study",keep.scale= "feature",label = TRUE)

#-----------------------------------------
#Q5
#top10MarkersAcrossConditions.png
DotPlot(pbmcSeuratObjectNormalized, features = rownames(top10_markers),
        cols = c("green", "orange"), dot.scale = 6, 
        split.by = "Condition_study") + RotatedAxis()

#-------------------------------
#conserved cell type markers across conditions.png
markersToPlot <- c("MMP2","COL4A2","RND1","IGFBP4","SELE","CD74","CLU","PRCP","APLNR","CLDN5","TFPI","RAMP3","ACKR1","ADGRF5","ADIRF","ANGPT2","ANGPTL2","APCDD1","PPP1R14A","SEMA3G","FBLN5","GJA5","PCSK5","ARL15","HLA-DRA","HLA-DPA1","HLA-DRB1","TFF3","CCL21","LYVE1","EFEMP1","MMRN1","SCNN1B","FN1","AQP1","AKR1C1","CEBPD","CPE","CYP1B1","ANLN","APOBEC3B","ARHGAP11A","ASF1B","ASPM","AURKB","CLEC14A","GSN","ABCG2","PLVAP","TXNIP","SPARC")

DotPlot(pbmcSeuratObjectNormalized, features = markersToPlot,
        cols = c("green", "orange"), dot.scale = 6, 
        split.by = "Condition_study") + RotatedAxis()

