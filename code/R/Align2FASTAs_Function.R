align2FASTAs <- function(FASTA_target_directory, FASTA_source_directory) {
  
  ##https://www.rdocumentation.org/packages/qdapTools/versions/1.3.3/topics/list2df
  list2df <- function(x) 
  { 
    MAX.LEN <- max(sapply(x, length), na.rm = TRUE) 
    DF <- data.frame(lapply(x, function(x) c(x, rep(NA, MAX.LEN - length(x))))) 
    colnames(DF) <- paste("V", seq(ncol(DF)), sep = "")   
    DF 
  }
  
  FASTA_source <- seqinr::read.fasta(file = FASTA_source_directory, seqtype = "AA", as.string = T, whole.header = T, seqonly = F)
  df_FASTA_source <-  list2df(FASTA_source)
  colnames(df_FASTA_source) <-  names(FASTA_source)
  character_FASTA_source <- apply(X = df_FASTA_source[1,], MARGIN = 1, FUN = as.character)
  
  FASTA_target <- seqinr::read.fasta(file = FASTA_target_directory, seqtype = "AA", as.string = T, whole.header = T, seqonly = F)
  df_FASTA_target <- list2df(FASTA_target)
  colnames(df_FASTA_target) <- names(FASTA_target)
  character_FASTA_target <- apply(X = df_FASTA_target[1,], MARGIN = 1, FUN = as.character)
  
  match <- c()
  max_scores <- c()
  for (i in 1:dim(character_FASTA_target)[1]) {
    
    alignment <- Biostrings::pairwiseAlignment(pattern = character_FASTA_source[,1],
                                               subject = character_FASTA_target[i,1])
    scores <- alignment@score
    max_scores[i] <- max(alignment@score)
    match[i] <- names(FASTA_source[which(scores == max(alignment@score))])
    
  }
  return(cbind(match, max_scores))
}
