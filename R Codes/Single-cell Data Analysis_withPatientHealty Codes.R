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
remove(pbmc)
remove(metadata)

mainSeuratObject  <- CreateSeuratObject(counts = pbmc_Copy , meta.data =  metadata_Copy, min.cells = 3, min.genes = 500 )

#----------------
FeatureScatter(mainSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(mainSeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#As the number of tumor cells (N = 12,213) greatly exceeded that of normal cells (N = 827),
#for significantly changed cell types,
#the down-sampling analysis to remove the influence of this unbalanced sample sizes, is essential.
#--------------------------------------------------------------------------------
#Create Healthy subset
#--------------------------------------------------------------------------------
subSetHealthyCells=subset(mainSeuratObject, subset = Condition_study == "healthy")# , downsample = 1000
subSetHealthyCells$orig.ident <- "healthy"
dim(subSetHealthyCells)
subSetHealthyCells <- NormalizeData(subSetHealthyCells, verbose = FALSE)
subSetHealthyCells <- FindVariableFeatures(subSetHealthyCells, selection.method = "vst", nfeatures = 2000)

#--------------------------------------------------------------------------------
#Create Patient subset
#--------------------------------------------------------------------------------
subSetPatientCells=subset(mainSeuratObject, subset = Condition_study == "tumor")#, downsample = 90
subSetPatientCells$orig.ident <- "tumor"
dim(subSetPatientCells)
subSetPatientCells <- NormalizeData(subSetPatientCells, verbose = FALSE)
subSetPatientCells <- FindVariableFeatures(subSetPatientCells, selection.method = "vst", nfeatures = 2000)
#---------------------------------------------------
#length(VariableFeatures(subSetHealthyCells))
#length(VariableFeatures(subSetPatientCells))
head(VariableFeatures(subSetHealthyCells), 10)
head(VariableFeatures(subSetPatientCells), 10)
#--------------------------------------------------------------------------------

#To select dim in FindIntegrationAnchors we need run jackStraw
#  nUMI = nCount_RNA
subSetHealthyCells <- ScaleData(object = subSetHealthyCells, vars.to.regress = c("nCount_RNA", "percent.mt"))
subSetHealthyCells <- RunPCA(object = subSetHealthyCells, pc.genes = VariableFeatures(subSetHealthyCells), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# pdf("PCA_Healthy.pdf")
# subSetHealthyCells <- ProjectDim (object = subSetHealthyCells)
# dev.off()
#.................................................
subSetHealthyCells <- JackStraw(object = subSetHealthyCells, num.replicate = 100)
# default num.pc is 50, but in this data is 20
ScoreJackStraw(
  subSetHealthyCells,
  reduction = "pca",
  dims = 1:20,
  score.thresh = 1e-05,
  do.plot = TRUE
)
#JackStrawPlotHealthy.png
JackStrawPlot(object = subSetHealthyCells, dims = 1:20)

#--------------------
#--------------------
subSetPatientCells <- ScaleData(object = subSetPatientCells, vars.to.regress = c("nCount_RNA", "percent.mt"))
subSetPatientCells <- RunPCA(object = subSetPatientCells, pc.genes = VariableFeatures(subSetPatientCells), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# pdf("PCA_Patient.pdf")
# subSetPatientCells <- ProjectDim (object = subSetPatientCells)
# dev.off()
#.................................................
subSetPatientCells <- JackStraw(object = subSetPatientCells, num.replicate = 100)
ScoreJackStraw(
  subSetPatientCells,
  reduction = "pca",
  dims = 1:20,
  score.thresh = 1e-05,
  do.plot = TRUE
)
#JackStrawPlotPatient.png
JackStrawPlot(object = subSetPatientCells, dims = 1:20)



#Perform integration
#IntegrateData
#--------------------------------------------------------------------------------

integrateData.anchors <- FindIntegrationAnchors(object.list = list(subSetHealthyCells, subSetPatientCells), dims = 1:20)
integrateData.combined <- IntegrateData(anchorset = integrateData.anchors, dims = 1:20)
#--------------------------------------------------------------------------------
#Perform an integrated analysis
#Now  run a single integrated analysis on all cells!
#--------------------------------------------------------------------------------
DefaultAssay(integrateData.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
integrateData.combined <- ScaleData(integrateData.combined, verbose = FALSE)
integrateData.combined <- RunPCA(integrateData.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrateData.combined <- RunUMAP(integrateData.combined, reduction = "pca", dims = 1:20)
integrateData.combined <- RunTSNE(integrateData.combined, reduction = "pca", dims = 1:20)
integrateData.combined <- FindNeighbors(integrateData.combined, reduction = "pca", dims = 1:20)
integrateData.combined <- FindClusters(integrateData.combined, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)#default resolution 0.5
# In the refPaper: pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#result:Number of communities: 9 !
#--------------------------------------------------------------------------------
# Visualization
#UMAP_ClustersSplitBy Condition_study.png
p1 <- DimPlot(integrateData.combined, reduction = "umap", group.by = "Condition_study")
p2 <- DimPlot(integrateData.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
#-------------------
#UMAP_ClustersSplitBy Condition_study_splitedByCondition.png
DimPlot(integrateData.combined, reduction = "umap", split.by = "Condition_study",label = TRUE)
DimPlot(integrateData.combined, reduction = "umap", split.by = "Condition_study", group.by = "Final_annotation_general",label = TRUE)
#DimPlot(integrateData.combined, reduction = "umap", split.by = "Final_annotation_general",label = TRUE)


#--------------------------------------------------------------------------------
#Identify conserved cell type markers
#--------------------------------------------------------------------------------

for(i in c(0:10)){
  canonicalCellTypeMarkerGenes <- FindConservedMarkers(integrateData.combined, ident.1 = i, grouping.var = "Condition_study", verbose = TRUE)
  head(canonicalCellTypeMarkerGenes)
  write.csv(head(canonicalCellTypeMarkerGenes),paste0("H:/tsk1File/data/resultsConsiderPatients/CanonicalCellTypeMarkerGenes/canonicalCellTypeMarkerGenesInCluster_", i, ".csv") , quote = F)
}

#...................
#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(integrateData.combined, features = c("SELL", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9", split.by = "Condition_study",label = TRUE)






FeatureScatter(integrateData.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
FeatureScatter(integrateData.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Condition_study")
 #GenesCount_CellsCountInClusters.png
VlnPlot(integrateData.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#---------------------------------------------------  
#Identification of highly variable features (feature selection)
#---------------------------------------------------  
integrateData.combined <- FindVariableFeatures(integrateData.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(integrateData.combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(integrateData.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#---------------------------------------------------  
#Scaling the data
#---------------------------------------------------  
allGenes <- rownames(integrateData.combined)
#integrateData.combined <- ScaleData(integrateData.combined, features = allGenes)

#VizDimLoadings(integrateData.combined, dims = 1:20, reduction = "pca")

#---------------------------------------------------  


integrateData.combined.allMarkers <- FindAllMarkers(integrateData.combined, assay = "RNA", test.use = "bimod", only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.25)
integrateData.combined.onlyPositiveMarkers <- FindAllMarkers(integrateData.combined, assay = "RNA", test.use = "bimod", only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)

dim(integrateData.combined.allMarkers) 
dim(integrateData.combined.onlyPositiveMarkers) 
# Write results to file
write.csv(integrateData.combined.allMarkers, "H:/tsk1File/data/resultsConsiderPatients/all_markers.csv", quote = F)
write.csv(integrateData.combined.onlyPositiveMarkers, "H:/tsk1File/data/resultsConsiderPatients/onlyPositiveMarkers.csv", quote = F)
#======================================================
# 
# #.........................................................
# # Return top 10 markers for cluster specified 'x'
 gen_marker_table <- function(x){
   integrateData.combined.allMarkers[integrateData.combined.allMarkers$cluster == x, ] %>%
     head(n=10)}
# 
# # create a data frame of results for clusters 0-10
 top10_markers <- map_dfr(0:10, gen_marker_table)
# 
#--------------------------------------------
# 
#  Write results to file
 write.csv(top10_markers, "H:/tsk1File/data/resultsConsiderPatients/top10_markers.csv", quote = F)
 write.csv(integrateData.combined.allMarkers, "H:/tsk1File/data/resultsConsiderPatients/integrateData.combinedMarkers.csv", quote = F)

#-------------------


 #CCL21SplitedWithConditions_UMAP.png
FeaturePlot(object = integrateData.combined, 
            features = c("CCL21" , "MMRN1"), split.by = "Condition_study", label =  TRUE,
            reduction = "tsne")

#---------------------------------------------------------------------



#---------------------------------------------------
# test for doc (Cluster 0)
# for cluster 0  (ConservedMarkers)
FeaturePlot(integrateData.combined, features = c("MMP2","COL4A2","RND1","IGFBP4","SELE","CD74"), min.cutoff = "q9") 

FeaturePlot(integrateData.combined, features = c("TFF3","CCL21","LYVE1","EFEMP1","MMRN1","SCNN1B"), min.cutoff = "q9") 

#---------------------------------------------------
#Plot for "CCL21" and "MMRN1" in ten clusters for Q4
#CCL21SplitedWithConditions_PCA.png

FeaturePlot(integrateData.combined, features = c("CCL21","MMRN1"), min.cutoff = "q9" , split.by = "Condition_study",keep.scale= "feature",label = TRUE)
FeaturePlot(integrateData.combined, features = c("CCL21","MMRN1"), min.cutoff = "q9" , split.by = "Final_annotation_general",keep.scale= "feature",label = TRUE)

FeaturePlot(integrateData.combined, features = c("CD74","CCL21","MMRN1","SCNN1B"), min.cutoff = "q9" , split.by = "Condition_study",keep.scale= "feature",label = TRUE)

#======================================================================
#DefaultAssay(integrateData.combined) <- "integrated"
#Identify differentially expressed genes between healthy and cancer patients
#Idents(object = integrateData.combined) <- integrateData.combined@meta.data$Condition_study

# 
# degMarkers <- FindMarkers(integrateData.combined, ident.1 = "healthy", ident.2 = "tumor")
# #  Write results to file
# write.csv(degMarkers, "H:/tsk1File/data/resultsConsiderPatients/findAllMarkersResluts.csv", quote = F)



#-----------------------------------------
#Q5
#topMarkersAcrossConditions.png
DotPlot(integrateData.combined, features = rownames(top10_markers),
        cols = c("green", "orange"), dot.scale = 8, 
        split.by = "Condition_study") + RotatedAxis()

#-------------------------------
#conserved cell type markers across conditions.png
markersToPlot <- c("EPCAM", "SFN", "cytokeratins","CCL14", "TK1","MMP2","COL4A2","RND1","IGFBP4","SELE","CD74","CLU","PRCP","APLNR","CLDN5","CCL21","MMRN1","SCNN1B")

DotPlot(integrateData.combined, features = markersToPlot,
        cols = c("green", "orange"), dot.scale = 8,
        split.by = "Condition_study") + RotatedAxis()

#-------------------------------------------------------
Idents(integrateData.combined)
integrateData.combined <- RenameIdents(integrateData.combined, '0'='capillary','1'='vein','2'='capillary','3'='vein','4'='arteriole','5'='vein', '6'='lymphatic', '7'='capillary','8'='vein','9'='other','10'='vein')
DimPlot(integrateData.combined, reduction = "umap", split.by = "Condition_study",label = TRUE)