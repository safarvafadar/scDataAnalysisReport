library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)

#---------------------------------------------------
metadata=readRDS("H:/tsk1File/data/esophagus_ID16538_metadataMatrix.rds")
pbmc <- readRDS(file = "H:/tsk1File/data/esophagus_ID16021_dataMatrixNormalized.rds")

#View(metadata[0:5,0:5])
#View(pbmc[0:5,0:5])
#---------------------------------------------------
#rownames(pbmc)=pbmc[,1]
#--- if has exception
  pbmc_Copy <- data.frame(pbmc)            # Create copy of data
  #View(pbmc_Copy[0:5,0:5])
  rownames(pbmc_Copy)=pbmc_Copy[,1]
#---------------------------------------------------
  pbmc_Copy =pbmc_Copy[,-1]
  #View(pbmc_Copy[0:5,0:5])
  pbmcSeuratObject  <- CreateSeuratObject(counts = pbmc_Copy , meta.data =  metadata)  #count=?1111? pbmc
#---------------------------------------------------

VlnPlot(object = pbmcSeuratObject, features = c('CCL21', 'MMRN1', 'A1BG','A1BG-AS1'))
#---------------------------------------------------  
mtPercentage=PercentageFeatureSet(pbmcSeuratObject, pattern = "^MT.")
  View(mtPercentage)
  # we have NA 
  grep ("^MT-", rownames(pbmcSeuratObject[["RNA"]]),value = T) # check if our dataset has mitochondrial genes
  grep ("^MT.", rownames(pbmcSeuratObject[["RNA"]]),value = T) 
