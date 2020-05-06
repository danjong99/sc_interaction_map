rm(list=ls())

library(dplyr)
lr_pair <- read.delim("~/Desktop/Data/All/interaction_map/ligand_receptor_mtx/LR_Pair_mus_musculus.txt")
epithelial_ligand_mtx <- read.delim("~/Desktop/Data/All/interaction_map/ligand_receptor_mtx/epithelial_ligand_mtx.txt")
epithelial_receptor_mtx <- read.delim("~/Desktop/Data/All/interaction_map/ligand_receptor_mtx/epithelial_receptor_mtx.txt")
pp_apc_ligand_mtx <- read.delim("~/Desktop/Data/All/interaction_map/ligand_receptor_mtx/pp_apc_ligand_mtx.txt")
pp_apc_receptor_mtx <- read.delim("~/Desktop/Data/All/interaction_map/ligand_receptor_mtx/pp_apc_receptor_mtx.txt")

# pp_apcs : new cluster
# (2,4,6): ILC13 --> 2
# (1): MPs --> 1
# (0,3): DCs --> 0
# (5,9): B cell --> exclude
# (7): pDCs --> 4
# (8): ILC2 --> 3

pp_apc_ligand_mtx$ident <- pp_apc_receptor_mtx$ident %>% unname()
pp_apc_ligand_mtx <- pp_apc_ligand_mtx[pp_apc_ligand_mtx$ident != 5 & 
                                         pp_apc_ligand_mtx$ident != 9,]
pp_apc_receptor_mtx$ident <- pp_apc_receptor_mtx$ident %>% unname()
pp_apc_receptor_mtx <- pp_apc_receptor_mtx[pp_apc_receptor_mtx$ident != 5 &
                                             pp_apc_receptor_mtx$ident != 9,]

pp_apc_ligand_mtx$ident <- gsub("^4$","2",pp_apc_ligand_mtx$ident)
pp_apc_ligand_mtx$ident <- gsub("^6$","2",pp_apc_ligand_mtx$ident)
pp_apc_ligand_mtx$ident <- gsub("^3$","0",pp_apc_ligand_mtx$ident)
pp_apc_ligand_mtx$ident <- gsub("^7$","4",pp_apc_ligand_mtx$ident)
pp_apc_ligand_mtx$ident <- gsub("^8$","3",pp_apc_ligand_mtx$ident)
pp_apc_ligand_mtx = droplevels(pp_apc_ligand_mtx)
pp_apc_ligand_mtx$ident <- as.numeric(pp_apc_ligand_mtx$ident)
pp_apc_ligand_mean_mtx <- pp_apc_ligand_mtx %>% group_by(ident) %>% summarise_all(funs(mean))
pp_apc_ligand_mtx2 <- pp_apc_ligand_mtx[,2:ncol(pp_apc_ligand_mtx)]
pp_apc_ligand_mtx2[pp_apc_ligand_mtx2>0] <- 1
pp_apc_ligand_mtx2 <- cbind(pp_apc_ligand_mtx[,1], pp_apc_ligand_mtx2)
colnames(pp_apc_ligand_mtx2)[1] <- "ident"
pp_apc_ligand_PA_mtx <- pp_apc_ligand_mtx2 %>% group_by(ident) %>% summarise_all(funs(mean))


pp_apc_receptor_mtx$ident <- gsub("^4$","2",pp_apc_receptor_mtx$ident)
pp_apc_receptor_mtx$ident <- gsub("^6$","2",pp_apc_receptor_mtx$ident)
pp_apc_receptor_mtx$ident <- gsub("^3$","0",pp_apc_receptor_mtx$ident)
pp_apc_receptor_mtx$ident <- gsub("^7$","4",pp_apc_receptor_mtx$ident)
pp_apc_receptor_mtx$ident <- gsub("^8$","3",pp_apc_receptor_mtx$ident)
pp_apc_receptor_mtx = droplevels(pp_apc_receptor_mtx)
pp_apc_receptor_mtx$ident <- as.numeric(pp_apc_receptor_mtx$ident)
pp_apc_receptor_mean_mtx <- pp_apc_receptor_mtx %>% group_by(ident) %>% summarise_all(funs(mean))
pp_apc_receptor_mtx2 <- pp_apc_receptor_mtx[,2:ncol(pp_apc_receptor_mtx)]
pp_apc_receptor_mtx2[pp_apc_receptor_mtx2 > 0] <- 1
pp_apc_receptor_mtx2 <- cbind(pp_apc_receptor_mtx[,1], pp_apc_receptor_mtx2)
colnames(pp_apc_receptor_mtx2)[1] <- "ident"
pp_apc_receptor_PA_mtx <- pp_apc_receptor_mtx2 %>% group_by(ident) %>% summarise_all(funs(mean))


epithelial_ligand_mean_mtx <- epithelial_ligand_mtx %>% group_by(ident) %>% summarise_all(funs(mean))
epithelial_ligand_mtx2 <- epithelial_ligand_mtx[,2:ncol(epithelial_ligand_mtx)]
epithelial_ligand_mtx2[epithelial_ligand_mtx2 > 0] <- 1
epithelial_ligand_mtx2 <- cbind(epithelial_ligand_mtx[,1], epithelial_ligand_mtx2)
colnames(epithelial_ligand_mtx2)[1] <- "ident"
epithelial_ligand_PA_mtx <- epithelial_ligand_mtx2 %>% group_by(ident) %>% summarise_all(funs(mean))

epithelial_receptor_mean_mtx <- epithelial_receptor_mtx %>% group_by(ident) %>% summarise_all(funs(mean))
epithelial_receptor_mtx2 <- epithelial_receptor_mtx[,2:ncol(epithelial_receptor_mtx)]
epithelial_receptor_mtx2[epithelial_receptor_mtx2 > 0] <- 1
epithelial_receptor_mtx2 <- cbind(epithelial_receptor_mtx[,1], epithelial_receptor_mtx2)
colnames(epithelial_receptor_mtx2)[1] <- "ident"
epithelial_receptor_PA_mtx <- epithelial_receptor_mtx2 %>% group_by(ident) %>% summarise_all(funs(mean))
