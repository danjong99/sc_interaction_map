rm(list = ls())

library(Seurat)
library(dplyr)
library(patchwork)
library(sctransform)

## Function definition
get_field = function(string,field=1,delim="_", fixed=T) return(strsplit(string,delim, fixed=fixed)[[1]][field])

info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

fae_umis = read.delim("./data_resource/GSE92332_FAE_UMIcounts.txt.gz")
org_rankl = read.delim("./data_resource/GSE92332_Org_RANKL_UMIcounts.txt.gz")

fae_umis = fae_umis[!duplicated(unlist(lapply(rownames(fae_umis), get_field, 1,"_"))),]
rownames(fae_umis) = unlist(lapply(rownames(fae_umis), get_field, 1,"_"))
org_rankl = org_rankl[!duplicated(unlist(lapply(rownames(org_rankl), get_field, 1,"_"))),]
rownames(org_rankl) = unlist(lapply(rownames(org_rankl), get_field, 1,"_"))

info(sprintf("Data dimensions: %s" , paste(dim(fae_umis), collapse = "x")))
info(sprintf("Data dimensions: %s" , paste(dim(org_rankl), collapse = "x")))

fae = CreateSeuratObject(counts = fae_umis, project = "fae_dome", min.cells = 3, 
                         min.features = 200)
org_rankl = CreateSeuratObject(counts = org_rankl, project = "epi_org", min.cells = 3,
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

for (i in 2:length(epiMerged.list)) {
  epiMerged.list[[1]] <- SCTransform(epiMerged.list[[1]], verbose = TRUE)
}

epiMerged.features <- SelectIntegrationFeatures(object.list = epiMerged.list, nfeatures = 3000)
epiMreged.list <- PrepSCTIntegration(object.list = epiMerged.list, anchor.features = epiMerged.features)
reference_dataset <- c(2,3,4,5,6,7)

epiMerged.anchors <- FindIntegrationAnchors(object.list = epiMerged.list, normalization.method = "SCT", 
                                       anchor.features = epiMerged.features, reference = reference_dataset)
epiMerged.integrated <- IntegrateData(anchorset = epiMerged.anchors, normalization.method = "SCT")

epiMerged.integrated <- RunPCA(object = epiMerged.integrated, verbose = FALSE)
epiMerged.integrated <- RunUMAP(object = epiMerged.integrated, dims = 1:30)

epiMerged.list[[8]]

ident_table = data.frame(cluster = Idents(epiMerged.integrated))
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

epiMerged.integrated$exp <- ident_table$exp
DimPlot(object = epiMerged.integrated, reduction = "tsne", label = F,
        label.size = 8, pt.size = 0.7, split.by = 'exp')

d_table <- ident_table %>% group_by(cluster, celltype) %>% count() %>% summarize(count = n) 
d_table2 <- d_table %>% group_by(cluster) %>% 
  summarize(maxCell = max(count),
            maxCellname = celltype[which(count == max(count))],
            totalCellNo = sum(count),
            maxCellper = round(max(count)/sum(count)*100,1))

epiMerged.integrated$celltype <- ident_table$celltype


##### Align with Query Data
names(epiMerged.list)

fae.query <- epiMerged.list[["fae"]]
fae.anchors <- FindTransferAnchors(reference = epiMerged.integrated, query = fae.query, 
                                   dims = 1:30)
predictions <- TransferData(anchorset = fae.anchors, refdata = epiMerged.integrated$celltype, 
                            dims = 1:30)
fae.query <- AddMetaData(fae.query, metadata = predictions)

fae_ident_table = data.frame(ident = fae.query$predicted.id)
fae_ident_table2 = fae_ident_table %>% group_by(ident) %>% count() %>% summarize(count = n)

library(sctransform)

