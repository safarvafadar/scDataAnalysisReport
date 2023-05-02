library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)

library(stringr)

#---------------------------------------------------
metadata=readRDS("H:/tsk1File/data/esophagus_ID16538_metadataMatrix.rds")
pbmc <- readRDS(file = "H:/tsk1File/data/esophagus_ID16021_dataMatrixNormalized.rds")
#---------------------------------------------------
pbmc_Copy <- data.frame(pbmc)            # Create copy of data
rownames(pbmc_Copy)=pbmc_Copy[,1]

#The data frame should be exactly same as column names, and genes as row names. // replace ".." with ": "
colnames(pbmc_Copy) = gsub("\\.\\.", ": ", colnames(pbmc_Copy), fixed=FALSE) #str_replace_all(colnames(pbmc_Copy), "[[:punct:]]+", ": ")


metadata_Copy <- data.frame(metadata)            # Create copy of metadata
rownames(metadata_Copy)=metadata_Copy[,1]


#---------------------------------------------------
pbmc_Copy =pbmc_Copy[,-1]

pbmcSeuratObject  <- CreateSeuratObject(counts = pbmc_Copy , meta.data =  metadata_Copy)# , min.cells = 0  
#---------------------------------------------------

plotBefore=FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(pbmcSeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(pbmcSeuratObject, features = c('CCL21', 'MMRN1', 'A1BG','A1BG-AS1'))
#---------------------------------------------------  

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmcSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#--------------------------------------------------------------
#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

metadata %>% ggplot(aes(x=Condition_study, fill=Condition_study)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
#---------------------------------------------------  
# Visualize the number UMIs/transcripts per cell
#???
#----------------
metadata %>% 
  ggplot(aes(color=Condition_study, x=nbr_count, fill= Condition_study)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Number of Counts") +
  geom_vline(xintercept = 500)

#--------------------------------------------------- !!!!!!!!!!!!!!!! نمایش درسته؟ محور y چی میگه 
metadata %>% 
  ggplot(aes(color=Condition_study, x=nbr_features, fill= Condition_study)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Number of Genes") +
  geom_vline(xintercept = 500)
#-----------------------------------------------------
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
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
dim(pbmcSeuratObjectNormalized[["RNA"]]@scale.data) # before running is <0 x 0 matrix>

allGenes <- rownames(pbmcSeuratObjectNormalized)
pbmcSeuratObjectNormalized <- ScaleData(pbmcSeuratObjectNormalized, features = allGenes)

# because is scaled before!!! 
#Error in txtProgressBar(min = 0, max = max.block, style = 3, file = stderr()) : must have 'max' > 'min'

dim(pbmcSeuratObjectNormalized[["RNA"]]@scale.data) # after running is <14305 x 13040 matrix>
#---------------------------------------------------  
#linear dimensional reduction
#---------------------------------------------------
# if we dont scale the data, we have error : Error in PrepDR(object = object, features = features, verbose = verbose) : Data has not been scaled. Please run ScaleData and retry
pbmcSeuratObjectNormalized <- RunPCA(pbmcSeuratObjectNormalized, features = VariableFeatures(object = pbmcSeuratObjectNormalized))

# Examine and visualize PCA results a few different ways
print(pbmcSeuratObjectNormalized[["pca"]], dims = 1:5, nfeatures = 5)


VizDimLoadings(pbmcSeuratObjectNormalized, dims = 1:2, reduction = "pca")

DimPlot(pbmcSeuratObjectNormalized, reduction = "pca")

DimHeatmap(pbmcSeuratObjectNormalized, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmcSeuratObjectNormalized, dims = 1:15, cells = 500, balanced = TRUE)

#VizPCA(object = pbmcSeuratObjectNormalized, pcs.use = 1:2)
#---------------------------------------------------  
#Determine the ‘dimensionality’ of the dataset
#---------------------------------------------------  
# pbmcSeuratObjectNormalized <- JackStraw(pbmcSeuratObjectNormalized, num.replicate = 100)
# pbmcSeuratObjectNormalized <- ScoreJackStraw(pbmcSeuratObjectNormalized, dims = 1:20)
# 
# JackStrawPlot(pbmcSeuratObjectNormalized, dims = 1:15)
#---------------------------------------------------  

#---------------------------------------------------  

#---------------------------------------------------  

#---------------------------------------------------  

#---------------------------------------------------  

#---------------------------------------------------  

#---------------------------------------------------  