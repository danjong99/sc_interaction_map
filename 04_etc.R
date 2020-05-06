pop_pop_int <- read.delim("/Volumes/Extreme SSD/scRNAseq_analysisfile/10XResults/Receptor_Ligand_Interaction/sc_interaction_map/results/03_interaction_cyto_merged_consol.txt")
epi_ident <- read.delim("/Volumes/Extreme SSD/scRNAseq_analysisfile/10XResults/Receptor_Ligand_Interaction/sc_interaction_map/results/epi_ident.txt")
pp_apc_ident <- read.delim("/Volumes/Extreme SSD/scRNAseq_analysisfile/10XResults/Receptor_Ligand_Interaction/sc_interaction_map/results/pp_apc_ident.txt")

rownames(pp_apc_ident) <- pp_apc_ident$cluster
pp_apc_ident$oldcluster <- c("0","1","2","3","4")
colnames(pp_apc_ident) <- c("index","ident","cluster")
pp_apc_ident <- pp_apc_ident[c("cluster","index","ident")]

epi_ident$inds <- paste0("epi_R_",epi_ident$new_cluster)
new_epi_ident <- epi_ident[,c("new_cluster","inds","new_ident")]
new_epi_ident <- new_epi_ident[!duplicated(new_epi_ident$new_cluster),]
rownames(new_epi_ident) <- new_epi_ident$inds
colnames(new_epi_ident) <- c("cluster","index","ident")
combined_index <- rbind(pp_apc_ident, new_epi_ident)
write.table(rbind(pp_apc_ident, new_epi_ident),
            file = "combined_ident_table.txt", sep = '\t')

temp <- combined_index[ord_cells,][c("ident","cluster")]
rownames(temp) <- NULL
pop_pop_int <- pop_pop_int[,1:4]
pop_pop_int <- cbind(pop_pop_int,temp)
write.table(pop_pop_int, file = "03_interaction_cyto_merged_consol.txt", sep = '\t', row.names = F)
