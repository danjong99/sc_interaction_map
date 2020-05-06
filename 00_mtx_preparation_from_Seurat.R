## Interaction Map
rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clipr")
library(Seurat)
library(dplyr)
library(clipr)

## Data loading
cca_pp_mhcii <- readRDS("/Volumes/Extreme SSD/scRNAseq_analysisfile/10XResults/scRNAseq_PP/scRNA-PP/cca_pp_mhcii.rds")
epiMerged <- readRDS("/Volumes/Extreme SSD/scRNAseq_analysisfile/10XResults/SIgA_production/scRNAseq_epithelialCells/epiMerged.integrated.rds")


DimPlot(cca_pp_mhcii, label = T, reduction = 'tsne')
pp.cm.markers <- FindAllMarkers(object = cca_pp_mhcii, only.pos = T)
top10 <- pp.cm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
FeaturePlot(object = cca_pp_mhcii, features = c("Rora"),
            reduction = 'tsne')

DimPlot(epiMerged, label = T, reduction = 'umap')

## LR_Pair_database loading
lr_pair <- read.csv("../LR_Pair_txt.txt", sep=" ", header=T, stringsAsFactors=F)
colnames(lr_pair)[1:3] <- c("Ligand", "Receptor", "known")

## Conversion Function
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

hmligands_table = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                 values = unique(lr_pair$Ligand),
                 mart = human, attributesL = c("mgi_symbol"), 
                 martL = mouse, uniqueRows=T)
hmreceptors_table = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                 values = unique(lr_pair$Receptor),
                 mart = human, attributesL = c("mgi_symbol"), 
                 martL = mouse, uniqueRows=T)

temp <- hmligands_table[!duplicated(hmligands_table$HGNC.symbol),]
sum(duplicated(hmligands_table$HGNC.symbol))
temp <- temp[!duplicated(temp$MGI.symbol),]
rownames(temp) <- NULL
hmligands_table = temp

temp <- hmreceptors_table[!duplicated(hmreceptors_table$HGNC.symbol),]
sum(duplicated(hmreceptors_table$HGNC.symbol))
temp <- temp[!duplicated(temp$MGI.symbol),]
rownames(temp) <- NULL
hmreceptors_table = temp

rm(temp)

lr_pair = lr_pair[lr_pair$Ligand %in% hmligands_table$HGNC.symbol,]
lr_pair = lr_pair[lr_pair$Receptor %in% hmreceptors_table$HGNC.symbol,]

Lds <- apply(lr_pair,1,FUN = function(x){
  temp = hmligands_table[x[1] == hmligands_table$HGNC.symbol,][2]
  temp <- temp %>% unlist() %>% unname()
  return(temp)
  })

Receps <- apply(lr_pair,1,FUN = function(x){
  temp = hmreceptors_table[x[2] == hmreceptors_table$HGNC.symbol,][2]
  temp <- temp %>% unlist() %>% unname()
  return(temp)
})

lr_pair$Ligand <- Lds
lr_pair$Receptor <- Receps
rownames(lr_pair) <- paste0(lr_pair$Ligand,'_',lr_pair$Receptor)

dim(lr_pair)

#####################################################################################
## Interaction Matrix Preparation from Epithelial cells

fae.cell.ids <- Cells(epiMerged)[grep("F_",Cells(epiMerged))]
fae_epithelium <- subset(x = epiMerged, cells = fae.cell.ids)
DimPlot(fae_epithelium, label = T)

dome_epi_mtx <- GetAssayData(fae_epithelium, assay = "RNA")
temp <- t(as.matrix(dome_epi_mtx))
temp <- cbind(data.frame(ident = Idents(fae_epithelium), stringsAsFactors = F), temp)
dome_epi_mtx <- t(temp)

col_name <- gsub("_GFI1B.","",colnames(dome_epi_mtx))
col_name <- gsub("F_B1_","",col_name)
col_name <- gsub("F_B2_","",col_name)
col_name <- gsub("F_B3_.","",col_name)
colnames(dome_epi_mtx) <- col_name
dome_epi_mtx <- dome_epi_mtx[,!duplicated(colnames(dome_epi_mtx))]
dome_epi_mtx <- t(dome_epi_mtx)
temp_idents <- dome_epi_mtx[,1:2]
dim(dome_epi_mtx)
epithelial_ligand_mtx <- dome_epi_mtx[,(colnames(dome_epi_mtx) %in% lr_pair$Ligand)]
epithelial_receptor_mtx <- dome_epi_mtx[,(colnames(dome_epi_mtx) %in% lr_pair$Receptor)]

epithelial_ligand_mtx <- cbind(as.data.frame(temp_idents), epithelial_ligand_mtx)
epithelial_receptor_mtx <- cbind(as.data.frame(temp_idents), epithelial_receptor_mtx)

write.table(lr_pair, file = "LR_Pair_mus_musculus.txt", sep = '\t')
write.table(epithelial_ligand_mtx, file = "epithelial_ligand_mtx.txt", sep = '\t')
write.table(epithelial_receptor_mtx, file = "epithelial_receptor_mtx.txt", sep = '\t')

#####################################################################################
## Interaction Matrix Preparation from PP Phagocytes

pp_apc_mtx <- GetAssayData(cca_pp_mhcii, assay = "RNA")
temp <- t(as.matrix(pp_apc_mtx))
temp <- cbind(data.frame(ident = Idents(cca_pp_mhcii), stringsAsFactors = F), temp)
pp_apc_mtx <- t(temp)

col_name <- gsub("PP1_","",colnames(pp_apc_mtx))
col_name <- gsub("PP2_","",col_name)
colnames(pp_apc_mtx) <- col_name
pp_apc_mtx <- pp_apc_mtx[,!duplicated(colnames(pp_apc_mtx))]
pp_apc_mtx <- t(pp_apc_mtx)
temp_idents <- pp_apc_mtx[,1:2]

pp_apc_ligand_mtx <- pp_apc_mtx[,colnames(pp_apc_mtx) %in% lr_pair$Ligand]
pp_apc_receptor_mtx <- pp_apc_mtx[,colnames(pp_apc_mtx) %in% lr_pair$Receptor]

pp_apc_ligand_mtx <- cbind(as.data.frame(temp_idents), pp_apc_ligand_mtx)
pp_apc_receptor_mtx <- cbind(as.data.frame(temp_idents), pp_apc_receptor_mtx)

write.table(pp_apc_ligand_mtx, file = "pp_apc_ligand_mtx.txt", sep = '\t')
write.table(pp_apc_receptor_mtx, file = "pp_apc_receptor_mtx.txt", sep = '\t')

