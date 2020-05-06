#### Defined Function 1
library(dplyr)

Gene_selector <- function(Expr_df, MeanUMI_df){
  #Expr_df <- Stroma_L_table
  #MeanUMI_df <- Stroma_L_Mean_table
  allsum <- data.frame(Expr_df[,1], stringsAsFactors=F)
  colnames(allsum)[1] <- 'ident'
  allsum$soap <- 1
  allsum <- allsum %>% group_by(ident) %>% summarise_all(sum)
  cc <- NULL
  #i <- grep("Csf1", colnames(Expr_df))
  
  for (i in 1:ncol(Expr_df)){
    positive_list <- data.frame(Expr_df[Expr_df[,i] > 0.1,1], stringsAsFactors=F)
    if (nrow(positive_list)==0){
      cc <- cbind(cc, as.logical(rep('FALSE', nrow(allsum))))
    } else {
      colnames(positive_list)[1] <- 'ident'
      positive_list$soap <- 1
      positive_sum <- positive_list %>% group_by(ident) %>% summarise_all(sum)
      
      if (nrow(positive_sum)!=nrow(allsum)){
        positive_sum2 <- merge(positive_sum, allsum, by='ident', all.y=T)
        positive_sum2[is.na(positive_sum2)] <- 0
        positive_sum2$soap.y <- NULL
        colnames(positive_sum2)[2] <- 'soap'
        positive_sum <- positive_sum2
      } else {}
      positive_sum$ratio <- positive_sum$soap/allsum$soap
      positive_sum$norm_ratio <- (positive_sum$soap/allsum$soap)*(sum(allsum$soap)/sum(positive_sum$soap))
      
      cc <- cbind(cc, positive_sum$ratio > 0.2)
    }
  }
  
  candidate_list <- which(cc==T, arr.ind=T)
  
  ret_df <- apply(candidate_list, 1, function(x){return(unlist(c(colnames(MeanUMI_df)[x[2]], MeanUMI_df$ident[x[1]],
                                                                 MeanUMI_df[x[1], x[2]])))})
  ret_df <- data.frame(t(ret_df), stringsAsFactors=F)
  ret_df <- ret_df[-grep("ident", ret_df[,1]),]
  colnames(ret_df) <- c('Gene', 'ident', 'MeanUMI')
  return(ret_df)
}

## Defined Function 2
Pair_finder <- function(LCandidates, RCandidates, LR_pair){
  pair_exist <- NULL
  for (i in (1:nrow(LR_pair))){
    if (((length(grep(LR_pair[i,1], LCandidates$Gene, ignore.case=T))!=0) & 
         (length(grep(LR_pair[i,2], RCandidates$Gene, ignore.case=T))!=0))==TRUE){
      pair_exist <- rbind(pair_exist, LR_pair[i,])
    } 
  }
  bobobo <- NULL
  for (k in 1:nrow(pair_exist)){
    if ((length(grep(pair_exist[k,1], LCandidates$Gene, ignore.case=T))==length(levels(factor(LCandidates$ident)))) &
        (length(grep(pair_exist[k,2], RCandidates$Gene, ignore.case=T))==length(levels(factor(RCandidates$ident))))){
    } else {
      for (LL in grep(pair_exist[k,1], LCandidates$Gene, ignore.case=T)){
        for (RR in grep(pair_exist[k,2], RCandidates$Gene, ignore.case=T)){
          bobobo <- rbind(bobobo, cbind(LCandidates[LL,], RCandidates[RR,]))
        }
      }
    }
  }
  colnames(bobobo) <- c("L_Gene", "L_ident", "L_MeanUMI",
                        "R_Gene", "R_ident", "R_MeanUMI")
  return(bobobo)
}
