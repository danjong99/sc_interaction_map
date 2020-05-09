library(Seurat)
library(dplyr)

## Function definition
get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])

info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

org_rankl = read.delim("./data_resource/GSE92332_Org_RANKL_UMIcounts.txt.gz")
org_rankl = org_rankl[!duplicated(unlist(lapply(rownames(org_rankl), get_field, 1,"_"))),]
rownames(org_rankl) = unlist(lapply(rownames(org_rankl), get_field, 1,"_"))
info(sprintf("Data dimensions: %s" , paste(dim(org_rankl), collapse = "x")))
org_rankl = CreateSeuratObject(counts = org_rankl, project = "epi_org", min.cells = 3, min.features = 200)
VlnPlot(org_rankl, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


#### with Individual Normalization
org_rankl.list <- SplitObject(org_rankl, split.by = "orig.ident")
for (i in 1:length(org_rankl.list)) {
  org_rankl.list[[i]] <- NormalizeData(org_rankl.list[[i]], verbose = FALSE)
  org_rankl.list[[i]] <- FindVariableFeatures(org_rankl.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}

org_rankl.anchors <- FindIntegrationAnchors(object.list = org_rankl.list, dims = 1:30)
org_rankl.integrated <- IntegrateData(anchorset = org_rankl.anchors, dims = 1:30)
DefaultAssay(org_rankl.integrated) <- "integrated"
org_rankl.integrated <- ScaleData(org_rankl.integrated, verbose = FALSE)
org_rankl.integrated <- RunPCA(org_rankl.integrated, npcs = 30, verbose = FALSE)

org_rankl.integrated <- FindNeighbors(object = org_rankl.integrated, dims = 1:30)
org_rankl.integrated <- FindClusters(object = org_rankl.integrated, resolution = 0.5)
head(x = Idents(object = org_rankl.integrated), 5)
org_rankl.integrated <- RunUMAP(org_rankl.integrated, reduction = "pca", dims = 1:30)
org_rankl.integrated <- RunTSNE(org_rankl.integrated, reduction = "pca", dims = 1:30)
DimPlot(object = org_rankl.integrated, reduction = "tsne", label = TRUE,
        label.size = 8, pt.size = 1)
org_rankl.all.cm.markers <- FindAllMarkers(object = org_rankl.integrated, only.pos = TRUE)
top20 <- org_rankl.all.cm.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

ident_table = data.frame(cluster = Idents(org_rankl.integrated))
ident_table$sample <- "1"
ident_table$cell.id <- "2"
ident_table$exp <- "3"
ident_table$celltype <- "4"

for (i in c(1:length(rownames(ident_table)))) {
  ad <- sapply(rownames(ident_table)[i], FUN = function(x){
    strsplit(x,"_")[[1]]
  }) %>% unlist() %>% unname()
  
  ident_table$sample[i] <- ad[1]
  ident_table$cell.id[i] <- ad[2]
  ident_table$exp[i] <- ad[3]
  ident_table$celltype[i] <- ad[4]
}

org_rankl.integrated$celltype <- ident_table$celltype
DimPlot(object = org_rankl.integrated, reduction = 'umap', label = TRUE,
        label.size = 4, group.by = 'celltype')
DimPlot(object = org_rankl.integrated, reduction = 'umap', label = TRUE,
        label.size = 4)

d_table <- ident_table %>% group_by(cluster, celltype) %>% count() %>% summarize(count = n) 
d_table2 <- d_table %>% group_by(cluster) %>% 
  summarize(maxCell = max(count),
            maxCellname = celltype[which(count == max(count))],
            totalCellNo = sum(count),
            maxCellper = round(max(count)/sum(count)*100,1))

#### Without individual normalization
VlnPlot(org_rankl, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
org_rankl <- subset(x = org_rankl, subset = nFeature_RNA > 200 & nFeature_RNA < 5000)
org_rankl <- NormalizeData(object = org_rankl, normalization.method = "LogNormalize", scale.factor = 10000)
org_rankl <- FindVariableFeatures(object = org_rankl, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = org_rankl), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = org_rankl)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(x = org_rankl)
org_rankl <- ScaleData(object = org_rankl, features = all.genes)
org_rankl <- RunPCA(object = org_rankl, features = VariableFeatures(object = org_rankl))
DimPlot(object = org_rankl, reduction = "pca", group.by = "orig.ident")
DimHeatmap(object = org_rankl, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(object = org_rankl)

org_rankl <- FindNeighbors(object = org_rankl, dims = 1:15)
org_rankl <- FindClusters(object = org_rankl, resolution = 0.5)
head(x = Idents(object = org_rankl), 5)

org_rankl <- RunUMAP(object = org_rankl, dims = 1:15)
org_rankl <- RunTSNE(object = org_rankl, dims = 1:15)

DimPlot(object = org_rankl, reduction = "umap", label = TRUE,
        label.size = 4, pt.size = 0.5)

ident_table = data.frame(cluster = Idents(org_rankl))
ident_table$sample <- "1"
ident_table$cell.id <- "2"
ident_table$exp <- "3"
ident_table$celltype <- "4"

for (i in c(1:length(rownames(ident_table)))) {
  ad <- sapply(rownames(ident_table)[i], FUN = function(x){
    strsplit(x,"_")[[1]]
  }) %>% unlist() %>% unname()
  
  ident_table$sample[i] <- ad[1]
  ident_table$cell.id[i] <- ad[2]
  ident_table$exp[i] <- ad[3]
  ident_table$celltype[i] <- ad[4]
}

org_rankl$celltype <- ident_table$celltype
org_rankl$exp <- ident_table$exp

library(ggplot2)

DimPlot(object = org_rankl, reduction = 'tsne', label = F,
        label.size = 4, group.by = 'celltype')
DimPlot(object = org_rankl, reduction = 'tsne', label = F,
        label.size = 4, group.by = 'celltype', split.by = 'orig.ident')

d_table <- ident_table %>% group_by(cluster, celltype) %>% count() %>% summarize(count = n) 
d_table2 <- d_table %>% group_by(cluster) %>% 
  summarize(maxCell = max(count),
            maxCellname = celltype[which(count == max(count))],
            totalCellNo = sum(count),
            maxCellper = round(max(count)/sum(count)*100,1))
