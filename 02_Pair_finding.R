source("./Fxns.R")

### Candidates Genes Selection
pp_candidates_L <- Gene_selector(pp_apc_ligand_mtx, pp_apc_ligand_mean_mtx)
pp_candidates_R <- Gene_selector(pp_apc_receptor_mtx, pp_apc_receptor_mean_mtx)
epi_candidates_L <- Gene_selector(epithelial_ligand_mtx, epithelial_ligand_mean_mtx)
epi_candidates_R <- Gene_selector(epithelial_receptor_mtx, epithelial_receptor_mean_mtx)

## Labeling change
pp_candidates_L$ident <- paste0('pp_L_', pp_candidates_L$ident)
pp_candidates_R$ident <- paste0('pp_R_', pp_candidates_R$ident)
epi_candidates_L$ident <- paste0('epi_L_', epi_candidates_L$ident)
epi_candidates_R$ident <- paste0('epi_R_', epi_candidates_R$ident)

## Significant Pair Finding
pp_L_epi_R <- Pair_finder(pp_candidates_L, epi_candidates_R, lr_pair)
pp_L_pp_R <- Pair_finder(pp_candidates_L, pp_candidates_R, lr_pair)
epi_L_pp_R <- Pair_finder(epi_candidates_L, pp_candidates_R, lr_pair) 
epi_L_epi_R <- Pair_finder(epi_candidates_L, epi_candidates_R, lr_pair)

merged <- rbind(pp_L_epi_R, pp_L_pp_R, epi_L_pp_R, epi_L_epi_R)
merged$Str <- (as.numeric(merged$L_MeanUMI)*as.numeric(merged$L_MeanUMI))/max((as.numeric(merged$L_MeanUMI)*as.numeric(merged$L_MeanUMI)))
merged$L_Gene <- paste(merged$L_ident, merged$L_Gene, sep='_')
merged$R_Gene <- paste(merged$R_ident, merged$R_Gene, sep='_')
write.table(merged, './results/01_interaction_cyto_merged.txt', sep='\t', quote=F, row.names=F)


