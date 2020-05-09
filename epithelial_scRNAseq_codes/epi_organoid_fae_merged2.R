library(Seurat)
library(dplyr)

epimerged.integrated <- readRDS("~/Desktop/Data/All/scRNAseq_epithelialCells/epimerged.integrated.rds")
epimerged.integrated@meta.data[,c("seurat_clusters","new_celltype")]
epimerged.integrated[["old.ident"]] <- Idents(object = epimerged.integrated)
#epimerged.integrated <- StashIdent(object = epimerged.integrated, save.name = "old.ident")

uniq.celltype <- unique(epimerged.integrated$old.ident) %>% sort()
meta_table <- epimerged.integrated@meta.data

b = as.character()
for (i in c(1:length(uniq.celltype))){
  a = which(meta_table$old.ident == uniq.celltype[i])[1]
  b = append(b, meta_table[a,]$new_celltype)
}
new_cell_table = data.frame(clusters = uniq.celltype, celltype = b)

new.cluster.ids <- new_cell_table$celltype
new.cluster.ids <- as.character(new.cluster.ids)
names(new.cluster.ids) <- new_cell_table$clusters
epimerged.integrated <- RenameIdents(epimerged.integrated, new.cluster.ids)
epimerged.integrated$sample <- 'FAE'
epimerged.integrated$sample[epimerged.integrated$exp != 'FAE'] <- 'ORG'
DimPlot(epimerged.integrated, reduction = "umap", 
        label = TRUE, pt.size = 0.5, split.by = 'sample') + NoLegend()
DimPlot(epimerged.integrated, reduction = "umap", 
        label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(object = epimerged.integrated, features = c("Ccl20"), reduction = 'umap',
            min.cutoff = 'q05', pt.size = 0.6, max.cutoff = 'q95', split.by = 'sample')

epimerged.all.markers <- FindAllMarkers(object = epimerged.integrated, only.pos = T)
top5 <- epimerged.all.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DefaultAssay(epimerged.integrated) <- 'integrated'
DoHeatmap(object = epimerged.integrated, features = top5$gene)
