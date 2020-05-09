#Ligand merger
testt <- merged

testt$L_MeanUMI <- as.numeric(testt$L_MeanUMI)
testt$R_MeanUMI <- as.numeric(testt$R_MeanUMI)
testt$Str <- as.numeric(testt$Str)

testt$L_Gene <- NULL
testt$R_ident <- NULL
#testt$L_ident_R_Gene <- paste(testt$R_Gene, testt$L_ident, sep='_')
#testt$R_Gene <- NULL
#testt$L_ident <- NULL

papapa <- testt %>% group_by(L_ident, R_Gene) %>% summarise_all(sum)

New_L_Gene <- NULL
for (i in 1:nrow(papapa)){
  sss <- merged[intersect(grep(papapa$L_ident[i], merged$L_ident), grep(papapa$R_Gene[i], merged$R_Gene)),]
  New_L_Gene <- append(New_L_Gene, paste(gsub('\\w*_[RL]_\\d*_(.*)', '\\1', sss$L_Gene), collapse='_'))
}

papapa$L_Gene <- paste(papapa$L_ident, New_L_Gene, sep='_')
#write.table(papapa, 'interaction_cyto_merged_new_Lmerged2.txt', sep='\t', quote=F, row.names=F)
#papapa <- read.csv('interaction_cyto_merged_new_Lmerged2.txt', sep='\t')
papapa$R_Ident <- gsub("(\\w_[LR]_\\d*)_.*", "\\1", papapa$R_Gene)
papapa <- papapa[sort(colnames(papapa))]
write.table(papapa, './results_pp_dEpi_orgEpi/02_interaction_cyto_merged_new_Lmerged_0.20.txt', sep='\t', quote=F, row.names=F)#Ligand merger


papapa2 <- data.frame(cbind(as.character(merged$L_ident), as.character(merged$R_ident)), stringsAsFactors=F)
papapa2$soap <- 1
colnames(papapa2) <- c("L_ident", "R_ident", "soap")
papapa2$L_ident <- gsub("_L_", "_R_", papapa2$L_ident)
papapa2 <- papapa2 %>% group_by(L_ident, R_ident) %>% summarise_all(funs(sum))
papapa2$soap_log2 <- log2(papapa2$soap)

write.table(papapa2, './results_pp_dEpi_orgEpi/03_interaction_cyto_merged_consol.txt', sep='\t', quote=F, row.names=F)

papapa2[grep("pp_",papapa2$R_ident),]
