rm(list = ls())

library(Seurat)
library(dplyr)
library(patchwork)

## Function definition
get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])

info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

fae_umis = read.delim("./data_resource/GSE92332_FAE_UMIcounts.txt.gz")
org_rankl_umis = read.delim("./data_resource/GSE92332_Org_RANKL_UMIcounts.txt.gz")

fae_umis = fae_umis[!duplicated(unlist(lapply(rownames(fae_umis), get_field, 1,"_"))),]
rownames(fae_umis) = unlist(lapply(rownames(fae_umis), get_field, 1,"_"))
org_rankl_umis = org_rankl[!duplicated(unlist(lapply(rownames(org_rankl_umis), get_field, 1,"_"))),]
rownames(org_rankl_umis) = unlist(lapply(rownames(org_rankl_umis), get_field, 1,"_"))

info(sprintf("Data dimensions: %s" , paste(dim(fae_umis), collapse = "x")))
info(sprintf("Data dimensions: %s" , paste(dim(org_rankl_umis), collapse = "x")))

fae = CreateSeuratObject(counts = fae_umis, project = "fae_dome", min.cells = 3, 
                         min.features = 200)
org_rankl = CreateSeuratObject(counts = org_rankl_umis, project = "epi_org", min.cells = 3,
                               min.features = 200)

VlnPlot(fae, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(org_rankl, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

fae$orig.ident <- paste("F",fae$orig.ident,sep = "_")
fae$norm_group <- "fae"
org_rankl$orig.ident <- paste("O",org_rankl$orig.ident,sep = "_")
org_rankl$norm_group <- org_rankl$orig.ident

epiMerged = merge(fae, org_rankl, add.cell.ids = c("F","O"))
VlnPlot(epiMerged, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,
        group.by = 'orig.ident')
epiMerged.list = SplitObject(epiMerged, split.by = "norm_group")

for (i in 1:length(epiMerged.list)) {
  epiMerged.list[[i]] <- NormalizeData(epiMerged.list[[i]], verbose = FALSE)
  epiMerged.list[[i]] <- FindVariableFeatures(epiMerged.list[[i]], selection.method = "vst", 
                                              nfeatures = 2000, verbose = FALSE)
}

names(epiMerged.list)
epiMerged.list <- epiMerged.list[c("fae","O_B2","O_B3","O_B4","O_B5","O_B6","O_B7")]
epiMerged.anchors <- FindIntegrationAnchors(object.list = epiMerged.list, dims = 1:30,
                                            reference = c(2,3,4,5,6,7))
epiMerged.integrated <- IntegrateData(anchorset = epiMerged.anchors, dims = 1:30)
DefaultAssay(epiMerged.integrated) <- "integrated"
epiMerged.integrated <- ScaleData(epiMerged.integrated, verbose = FALSE)
epiMerged.integrated <- RunPCA(epiMerged.integrated, npcs = 30, verbose = FALSE)

epiMerged.integrated <- FindNeighbors(object = epiMerged.integrated, dims = 1:30)
epiMerged.integrated <- FindClusters(object = epiMerged.integrated, resolution = 0.5)
head(x = Idents(object = epiMerged.integrated), 5)
epiMerged.integrated <- RunUMAP(epiMerged.integrated, reduction = "pca", dims = 1:30)
epiMerged.integrated <- RunTSNE(epiMerged.integrated, reduction = "pca", dims = 1:30)
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = TRUE,
        label.size = 8, pt.size = 0.7)
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = TRUE,
        label.size = 8, pt.size = 0.7, split.by = 'norm_group')

orig_ident = epiMerged.integrated$orig.ident
org_rankl = subset(epiMerged.integrated, 
                   cells = names(orig_ident[orig_ident != "F_B1" & orig_ident != "F_B2" & orig_ident != "F_B3"]))


exp = sapply(rownames(epiMerged.integrated@meta.data),
       FUN = function(x) strsplit(x,"_")[[1]][4]) %>% unlist() %>% unname()
epiMerged.integrated$exp = exp

ident_table = data.frame(cluster = Idents(epimerged.integrated))
ident_table$sample <- "1"
ident_table$cell.id <- "2"
ident_table$exp <- "3"
ident_table$celltype <- "4"

for (i in c(1:length(rownames(ident_table)))) {
  ad <- sapply(rownames(ident_table)[i], FUN = function(x){
    strsplit(x,"_")[[1]]
  }) %>% unlist() %>% unname()
  
  ident_table$sample[i] <- ad[2]
  ident_table$cell.id[i] <- ad[3]
  ident_table$exp[i] <- ad[4]
  ident_table$celltype[i] <- ad[5]
}

d_table <- ident_table %>% group_by(cluster, celltype) %>% count() %>% summarize(count = n) 
d_table2 <- d_table %>% group_by(cluster) %>% 
  summarize(maxCell = max(count),
            maxCellname = celltype[which(count == max(count))],
            totalCellNo = sum(count),
            maxCellper = round(max(count)/sum(count)*100,1))

##### Align with Query Data
names(epiMerged.list)

fae.query <- epiMerged.list[["fae"]]
fae.anchors <- FindTransferAnchors(reference = org_rankl, query = fae.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = fae.anchors, refdata = org_rankl$celltype, 
                            dims = 1:30)
fae.query <- AddMetaData(fae.query, metadata = predictions)
fae_ident_table = data.frame(ident = fae.query$predicted.id)
fae_ident_table2 = fae_ident_table %>% group_by(ident) %>% count() %>% summarize(count = n)
fae_ident_table <- droplevels(fae_ident_table)
fae_ident_table$ident
celltype <- sapply(names(epiMerged.integrated$celltype[is.na(epiMerged.integrated$celltype)]),
       FUN = function(x){ fae_ident_table[x,][1] }) %>% unlist() %>% unname()
celltype <- as.character(celltype)
ident_table$celltype[1:length(celltype)] <- celltype
epiMerged.integrated$celltype <- ident_table$celltype

ident_table$exp <- gsub("GFI1B.","FAE",ident_table$exp)
ident_table$exp <- factor(ident_table$exp, levels = c("FAE","Control","RANKL.Day3","RANKL.Day6"))
epiMerged.integrated$exp <- ident_table$exp

aaa = as.character(Idents(epiMerged.integrated))
new_cell_type = sapply(aaa, FUN = function(x){ d_table2[d_table2$cluster == x,]["maxCellname"] }) %>% unlist() %>% unname()
epiMerged.integrated$new_celltype <- new_cell_type
epiMerged.integrated$sample <- 'FAE'
epiMerged.integrated$sample[epiMerged.integrated$exp != 'FAE'] <- 'ORG'

DimPlot(object = epiMerged.integrated, reduction = "tsne", label = F,
        label.size = 8, pt.size = 0.7)
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = T,
        label.size = 8, pt.size = 0.7, group.by = "new_celltype") + NoLegend()
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = F,
        label.size = 8, pt.size = 0.7, 
        split.by = 'exp', group.by = 'new_celltype')
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = T,
        label.size = 5, pt.size = 0.7, 
        split.by = 'sample', group.by = 'new_celltype')

DefaultAssay(epiMerged.integrated) <- "RNA"
epiMerged.integrated <- NormalizeData(epiMerged.integrated, verbose = TRUE)
FeaturePlot(object = epiMerged.integrated, features = c("Sox8"), reduction = 'tsne',
            min.cutoff = 'q05', pt.size = 0.6, max.cutoff = 'q95')
saveRDS(epimerged.integrated, file = 'epimerged.integrated.rds')

epimerged.all.markers <- FindAllMarkers(object = epimerged.integrated, only.pos = T)
saveRDS(epimerged.all.markers, file = 'epimerged.all.markers.rds')

top10 <- epimerged.all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
epimerged.small = subset(x = epimerged.integrated, downsample = 100)
DefaultAssay(epimerged.small) <- 'integrated'
DoHeatmap(object = epimerged.small, features = top10$gene) + NoLegend()
FeaturePlot(object = epimerged.integrated, features = c("Ccl20"),
            min.cutoff = 'q05', max.cutoff = 'q95')

