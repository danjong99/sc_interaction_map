library(Seurat)
library(dplyr)
library(cowplot)

## Function definition
get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])

info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

## Downloading UMI count data
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz", destfile="GSE92332_atlas_UMIcounts.txt.gz")
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_FAE_UMIcounts.txt.gz", destfile="GSE92332_FAE_UMIcounts.txt.gz")
## Reading UMI count data from file
atlas_umis = read.delim("GSE92332_atlas_UMIcounts.txt.gz")
fae_umis = read.delim("GSE92332_FAE_UMIcounts.txt.gz")
org_rankl = read.delim("./data_resource/GSE92332_Org_RANKL_UMIcounts.txt.gz")

atlas_umis = atlas_umis[!duplicated(unlist(lapply(rownames(atlas_umis), get_field, 1,"_"))),]
rownames(atlas_umis) = unlist(lapply(rownames(atlas_umis), get_field, 1,"_"))
fae_umis = fae_umis[!duplicated(unlist(lapply(rownames(fae_umis), get_field, 1,"_"))),]
rownames(fae_umis) = unlist(lapply(rownames(fae_umis), get_field, 1,"_"))
org_rankl = org_rankl[!duplicated(unlist(lapply(rownames(org_rankl), get_field, 1,"_"))),]
rownames(org_rankl) = unlist(lapply(rownames(org_rankl), get_field, 1,"_"))

info(sprintf("Data dimensions: %s" , paste(dim(atlas_umis), collapse = "x")))
info(sprintf("Data dimensions: %s" , paste(dim(fae_umis), collapse = "x")))
info(sprintf("Data dimensions: %s" , paste(dim(org_rankl), collapse = "x")))

atlas = CreateSeuratObject(counts = atlas_umis, project = "epi_atlas", min.cells = 3, min.features = 200)
atlas
fae = CreateSeuratObject(counts = fae_umis, project = "fae_dome", min.cells = 3, min.features = 200)
fae
org_rankl = CreateSeuratObject(counts = org_rankl, project = "epi_org", min.cells = 3, min.features = 200)
org_rankl

atlas[["percent.mt"]] <- PercentageFeatureSet(atlas, pattern = "^mt-")
fae[["percent.mt"]] <- PercentageFeatureSet(fae, pattern = "^mt-")
org_rankl[["percent.mt"]] <- PercentageFeatureSet(org_rankl, pattern = "^mt-")

VlnPlot(atlas.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(fae.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(org_rankl, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

atlas.list <- SplitObject(atlas, split.by = "orig.ident")
for (i in 1:length(atlas.list)) {
  atlas.list[[i]] <- NormalizeData(atlas.list[[i]], verbose = FALSE)
  atlas.list[[i]] <- FindVariableFeatures(atlas.list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}
reference.list <- atlas.list[c("B1","B2","B3","B5","B7","B8","B9","B10")]
atlas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
atlas.integrated <- IntegrateData(anchorset = atlas.anchors, dims = 1:30)
DefaultAssay(atlas.integrated) <- "integrated"
atlas.integrated <- ScaleData(atlas.integrated, verbose = FALSE)
atlas.integrated <- RunPCA(atlas.integrated, npcs = 30, verbose = FALSE)

atlas.integrated <- FindNeighbors(object = atlas.integrated, dims = 1:30)
atlas.integrated <- FindClusters(object = atlas.integrated, resolution = 0.8)
head(x = Idents(object = atlas.integrated), 5)
atlas.integrated <- RunUMAP(atlas.integrated, reduction = "pca", dims = 1:30)
DimPlot(object = atlas.integrated, reduction = "umap", label = TRUE,
        label.size = 8, pt.size = 1)
atlas.all.cm.markers <- FindAllMarkers(object = atlas.integrated, only.pos = TRUE)
top50 <- atlas.all.cm.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
top1 <- atlas.all.cm.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
top10 <- atlas.all.cm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = atlas.integrated, features = top10$gene) + NoLegend()

metadata <- data.frame( cellid = rownames(atlas.integrated@meta.data),
            celltype = unlist(lapply(rownames(atlas.integrated@meta.data), get_field, 3,"_")),
            cluster = atlas.integrated@meta.data$seurat_clusters,
            batch = unlist(lapply(rownames(atlas.integrated@meta.data), get_field, 1,"_")))

overlap = table(metadata[,2:3])            
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]
aheatmap(overlap, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
new.cluster.ids <- c("Enterocyte.Progenitor/Progenitor.Late","TA.Early","TA.G2","Enterocyte.Immature.Proximal",
                     "Stem","Stem","Goblet","Enterocyte.Mature.Proximal","Enterocyte.Immature.Distal",
                     "Enterocyte.Progenitor.Early","TA.G1","Endocrine","Enterocyte.Mature.Distal",
                     "Tuft","Paneth","Endocrine")
names(new.cluster.ids) <- levels(atlas.integrated)
atlas.integrated <- RenameIdents(atlas.integrated, new.cluster.ids)
DimPlot(object = atlas.integrated, reduction = "umap", label = TRUE,
        label.size = 4, pt.size = 0.7)
FeaturePlot(object = atlas.integrated, features = c("Nos2"),
            min.cutoff = 'q01')

CCR <- rownames(GetAssayData(atlas.integrated))[grep("Ccr", rownames(GetAssayData(atlas.integrated)))]
CXCR <- rownames(GetAssayData(atlas.integrated))[grep("Cxcr", rownames(GetAssayData(atlas.integrated)))]
CCL <- rownames(GetAssayData(atlas.integrated))[grep("Ccl", rownames(GetAssayData(atlas.integrated)))]
CXCL <- rownames(GetAssayData(atlas.integrated))[grep("Cxcl", rownames(GetAssayData(atlas.integrated)))]
genelist <- c(CCR, CXCR, CCL, CXCL)

DoHeatmap(object = atlas.integrated, 
          features = genelist) + NoLegend()

saveRDS(atlas.integrated, file = "atlas.integrated.rds")
saveRDS(atlas.all.cm.markers, file = "atlas.all.cm.markers.rds")


fae.list <- SplitObject(fae, split.by = "orig.ident")
for (i in 1:length(fae.list)) {
  fae.list[[i]] <- NormalizeData(fae.list[[i]], verbose = FALSE)
  fae.list[[i]] <- FindVariableFeatures(fae.list[[i]], selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
}
reference.list <- fae.list[c("B1","B2")]
fae.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
fae.integrated <- IntegrateData(anchorset = fae.anchors, dims = 1:30)
DefaultAssay(fae.integrated) <- "integrated"
fae.integrated <- ScaleData(fae.integrated, verbose = FALSE)
fae.integrated <- RunPCA(fae.integrated, npcs = 30, verbose = FALSE)

fae.integrated <- FindNeighbors(object = fae.integrated, dims = 1:30)
fae.integrated <- FindClusters(object = fae.integrated, resolution = 0.4)
head(x = Idents(object = fae.integrated), 5)
fae.integrated <- RunUMAP(fae.integrated, reduction = "pca", dims = 1:30)
fae.integrated <- RunTSNE(fae.integrated, reduction = "pca", dims = 1:30)
DimPlot(object = fae.integrated, reduction = "umap", label = TRUE,
        label.size = 8, pt.size = 1)
fae.all.cm.markers <- FindAllMarkers(object = fae.integrated, only.pos = TRUE)
top50 <- fae.all.cm.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
top1 <- fae.all.cm.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
top10 <- fae.all.cm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = fae.integrated, features = top10$gene) + NoLegend()

fae_metadata <- data.frame( cellid = rownames(fae.integrated@meta.data),
                        celltype = unlist(lapply(rownames(fae.integrated@meta.data), get_field, 3,"_")),
                        cluster = fae.integrated@meta.data$seurat_clusters,
                        batch = unlist(lapply(rownames(fae.integrated@meta.data), get_field, 1,"_")))
         
overlap = table(fae_metadata[,2:3])            
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]
aheatmap(overlap, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)

CCR <- rownames(GetAssayData(fae.integrated))[grep("Ccr", rownames(GetAssayData(fae.integrated)))]
CXCR <- rownames(GetAssayData(fae.integrated))[grep("Cxcr", rownames(GetAssayData(fae.integrated)))]
CCL <- rownames(GetAssayData(fae.integrated))[grep("Ccl", rownames(GetAssayData(fae.integrated)))]
CXCL <- rownames(GetAssayData(fae.integrated))[grep("Cxcl", rownames(GetAssayData(fae.integrated)))]
genelist <- c(CCR, CXCR, CCL, CXCL)

FeaturePlot(object = fae.integrated, features = c("Nos2"),
            min.cutoff = 'q10', cols = c("lightgrey","red"), pt.size = 0.65)
DotPlot(object = fae.integrated, features = c("Ccl6","Ccl9","Ccl20"),
        cols = c("lightgrey","red"), col.min = -0.5)
DotPlot(object = fae.integrated, features = c("Nos2","Il13ra1"),
        cols = c("lightgrey","red"), col.min = -0.5)

DotPlot(object = atlas.integrated, features = c("Ccl6","Ccl9","Ccl20"),
        cols = c("lightgrey","red"), col.min = -0.5)
DotPlot(object = atlas.integrated, features = c("Nos2","Il13ra1"),
        cols = c("lightgrey","red"), col.min = -0.5
        )

DoHeatmap(object = fae.integrated, 
          features = genelist)

saveRDS(fae.integrated, file = "fae.integrated.rds")
saveRDS(fae.all.cm.markers, file = "fae.all.cm.markers.rds")

epiMerged = merge(atlas, fae, add.cell.ids = c("A","F"))
VlnPlot(epiMerged, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
VlnPlot(atlas, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
epiMerged.list = SplitObject(epiMerged, split.by = "orig.ident")

for (i in 1:length(epiMerged.list)) {
  epiMerged.list[[i]] <- NormalizeData(epiMerged.list[[i]], verbose = FALSE)
  epiMerged.list[[i]] <- FindVariableFeatures(epiMerged.list[[i]], selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
}

reference.list <- epiMerged.list[c("B1","B2","B3","B5","B7","B8","B9","B10")]
epiMerged.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
epiMerged.integrated <- IntegrateData(anchorset = epiMerged.anchors, dims = 1:30)
DefaultAssay(epiMerged.integrated) <- "integrated"
epiMerged.integrated <- ScaleData(epiMerged.integrated, verbose = FALSE)
epiMerged.integrated <- RunPCA(epiMerged.integrated, npcs = 30, verbose = FALSE)

epiMerged.integrated <- FindNeighbors(object = epiMerged.integrated, dims = 1:30)
epiMerged.integrated <- FindClusters(object = epiMerged.integrated, resolution = 0.8)
head(x = Idents(object = epiMerged.integrated), 5)
epiMerged.integrated <- RunUMAP(epiMerged.integrated, reduction = "pca", dims = 1:30)
DimPlot(object = epiMerged.integrated, reduction = "umap", label = TRUE,
        label.size = 8, pt.size = 0.7)
epiMerged.all.cm.markers <- FindAllMarkers(object = epiMerged.integrated, only.pos = TRUE)
top50 <- epiMerged.all.cm.markers %>% group_by(cluster) %>% top_n(50, avg_logFC)
top1 <- epiMerged.all.cm.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
top10 <- epiMerged.all.cm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = epiMerged.integrated, features = top10$gene,size = 1,
          angle = 90) + NoLegend()

indexf = unlist(lapply(rownames(epiMerged.integrated@meta.data), get_field, 1,"_")) == "F"
indexa = unlist(lapply(rownames(epiMerged.integrated@meta.data), get_field, 1,"_")) == "A"
cellidf = rownames(epiMerged.integrated@meta.data)[indexf]
cellida = rownames(epiMerged.integrated@meta.data)[indexa]
tempf = subset(x = epiMerged.integrated, cells = cellidf)
tempa = subset(x = epiMerged.integrated, cells = cellida)

plot1 <- DimPlot(object = epiMerged.integrated, reduction = "umap",
        label = T, label.size = 4, cells = cellidf)
plot2 <- DimPlot(object = epiMerged.integrated, reduction = "umap",
        label = T, label.size = 4, cells = cellida)
plot_grid(plot1, plot2)

metadata <- data.frame( cellid = rownames(tempa@meta.data),
                        celltype = unlist(lapply(rownames(tempa@meta.data), get_field, 4,"_")),
                        cluster = tempa@meta.data$seurat_clusters,
                        batch = unlist(lapply(rownames(tempa@meta.data), get_field, 2,"_")))

overlap = table(metadata[,2:3])            
overlap = round(sweep(overlap,2,colSums(overlap),`/`),2)
overlap = overlap[,apply(overlap, 1, FUN=which.max)]
aheatmap(overlap, border_color = list("cell"="white"), txt=overlap, Colv = NA, Rowv = NA)
new.cluster.ids <- c(
  "Enterocyte.Progenitor.Early","TA.Early/G1/G2","Stem","Enterocyte.Progenitor.Late","Enterocyte.Immature.Distal",
  "Stem","TA.Early","Enterocyte.Immature.Proximal","Enterocyte.Mature.Proximal","TA.G1",
  "Goblet","Goblet","Tuft","Enterocyte.Mature.Distal","Endocrine","Paneth","Endocrine"
)
names(new.cluster.ids) <- levels(tempa)
epiMerged.integrated <- RenameIdents(epiMerged.integrated, new.cluster.ids)
disp_order <- levels(epiMerged.integrated)
DimPlot(object = epiMerged.integrated, reduction = "umap", label = TRUE,
              label.size = 4, pt.size = 0.7, order = disp_order[order(disp_order, decreasing = T)])
pa <- DimPlot(object = epiMerged.integrated, reduction = "umap", label = TRUE,
        label.size = 4, pt.size = 0.7, order = disp_order[order(disp_order, decreasing = T)],
        cells = cellida)
pf <- DimPlot(object = epiMerged.integrated, reduction = "umap", label = TRUE,
              label.size = 4, pt.size = 0.7, order = disp_order[order(disp_order, decreasing = T)],
              cells = cellidf)
devtools::install_github('sjessa/ggmin')
pa + ggmin::theme_powerpoint()

saveRDS(epiMerged.integrated, file = "epiMerged.integrated.rds")
saveRDS(epiMerged.all.cm.markers, file = "epiMerged.all.cm.markers.rds")

FeaturePlot(object = epiMerged.integrated, features = c("Nos2"),
            min.cutoff = 'q10', cols = c("lightgrey","red"), pt.size = 0.65)
epiMerged.integrated@meta.data$location <- unlist(lapply(rownames(epiMerged.integrated@meta.data), get_field, 1,"_"))

FeaturePlot(object = epiMerged.integrated, features = c(""),
        split.by = "location", min.cutoff = 'q10', cols = c("lightgrey","red"))
DotPlot(object = epiMerged.integrated, features = c("Nos2"), split.by = "location")
