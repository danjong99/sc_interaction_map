pop_pop_int <- read.delim("~/Desktop/Data/All/scRNAseq_epithelialCells/sc_interaction_map/results_pp_dEpi_orgEpi/03_interaction_cyto_merged_consol.txt")
epi_ident <- read.delim("~/Desktop/Data/All/scRNAseq_epithelialCells/sc_interaction_map/results_pp_dEpi_orgEpi/epi_ident.txt")
pp_apc_ident <- read.delim("~/Desktop/Data/All/scRNAseq_epithelialCells/sc_interaction_map/results_pp_dEpi_orgEpi/pp_apc_ident.txt")

rownames(pp_apc_ident) <- pp_apc_ident$cluster
pp_apc_ident$oldcluster <- c("0","1","2","3","4")
colnames(pp_apc_ident) <- c("index","ident","cluster")
pp_apc_ident <- pp_apc_ident[c("cluster","index","ident")]
write.table(pp_apc_ident,file = './results_pp_dEpi_orgEpi/pp_apc_ident.txt', sep = '\t', row.names = F)

epi_ident$inds <- paste0("epi_R_",epi_ident$cluster)
new_epi_ident <- epi_ident[,c("cluster","inds","ident")]
new_epi_ident <- new_epi_ident[!duplicated(new_epi_ident$new_cluster),]
rownames(new_epi_ident) <- new_epi_ident$inds
colnames(new_epi_ident) <- c("cluster","index","ident")
combined_index <- rbind(pp_apc_ident, new_epi_ident)
write.table(combined_index,
            file = "./results_pp_dEpi_orgEpi/combined_ident_table.txt", sep = '\t')

l_ord_cells <- as.character(pop_pop_int$L_ident)
r_ord_cells <- as.character(pop_pop_int$R_ident)

temp1 <- combined_index[l_ord_cells,][c("ident","cluster")]
temp2 <- combined_index[r_ord_cells,][c("ident","cluster")]
rownames(temp1) <- NULL; rownames(temp2) <- NULL
colnames(temp1) <- c("L_celltype","L_cluster")
colnames(temp2) <- c("R_celltype","R_cluster")
pop_pop_int <- pop_pop_int[,1:4]
pop_pop_int <- cbind(pop_pop_int,temp1,temp2)
colnames(pop_pop_int)
pop_pop_int <- pop_pop_int[c("L_ident","L_celltype","L_cluster","R_ident","R_celltype","R_cluster","soap","soap_log2")]
write.table(pop_pop_int, file = "./results_pp_dEpi_orgEpi/03_interaction_cyto_merged_consol.txt", sep = '\t', row.names = F)

table1 <- read.delim("~/Desktop/Data/All/scRNAseq_epithelialCells/sc_interaction_map/results_pp_dEpi_orgEpi/01_interaction_cyto_merged.txt")
table1[,c("L_celltype","L_cluster","R_celltype","R_cluster")] <- NULL
combined_index <- combined_ident_table
combined_index <- rbind(combined_index, combined_index)
a <- as.character(combined_index$index)
a[1:(length(a)/2)] <- gsub("R","L",a[1:(length(a)/2)])
rownames(combined_index) <- a
combined_index$index <- a

l_ord_cells <- as.character(table1$L_ident)
r_ord_cells <- as.character(table1$R_ident)

temp1 <- combined_index[l_ord_cells,][c("ident","cluster")]
temp2 <- combined_index[r_ord_cells,][c("ident","cluster")]
rownames(temp1) <- NULL; rownames(temp2) <- NULL
colnames(temp1) <- c("L_celltype","L_cluster")
colnames(temp2) <- c("R_celltype","R_cluster")
table1 <- cbind(table1,temp1,temp2)
colnames(table1)
table1 <- table1[c("L_Gene","L_ident","L_celltype","L_cluster","L_MeanUMI",
                   "R_Gene","R_ident","R_celltype","R_cluster","R_MeanUMI","Str")]
write.table(table1, file = "./results_pp_dEpi_orgEpi/01_interaction_cyto_merged.txt", sep = '\t', row.names = F)


table2 <- read.delim("~/Desktop/Data/All/scRNAseq_epithelialCells/sc_interaction_map/results_pp_dEpi_orgEpi/02_interaction_cyto_merged_new_Lmerged_0.20.txt")
table2[,c("L_celltype","L_cluster","R_celltype","R_cluster")] <- NULL
l_ord_cells <- as.character(table2$L_ident)
r_ord_cells <- as.character(table2$R_ident)

temp1 <- combined_index[l_ord_cells,][c("ident","cluster")]
temp2 <- combined_index[r_ord_cells,][c("ident","cluster")]
rownames(temp1) <- NULL; rownames(temp2) <- NULL
colnames(temp1) <- c("L_celltype","L_cluster")
colnames(temp2) <- c("R_celltype","R_cluster")
table2 <- cbind(table2,temp1,temp2)
colnames(table2)
table2 <- table2[c("L_Gene","L_ident","L_celltype","L_cluster","L_MeanUMI",
                   "R_Gene","R_ident","R_celltype","R_cluster","R_MeanUMI","Str")]
write.table(table1, file = "./results_pp_dEpi_orgEpi/02_interaction_cyto_merged_new_Lmerged_0.20.txt", sep = '\t', row.names = F)

