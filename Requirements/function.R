rankEnrich <- function(exp, genelist){
  if (!is.matrix(exp) || is.null(rownames(exp)) || length(genelist) == 0) {
    stop("Invalid input. 'exp' must be a non-empty matrix and 'genelist' must be a non-empty list.")
  }
  row_names = rownames(exp)
  
  num_genes = nrow(exp)
  
  num_samples = ncol(exp)
  
  R = matrixStats::colRanks(exp, preserveShape = T, ties.method = 'average')
  
  result_matrix = matrix(0, nrow = length(genelist), ncol = num_samples)
  rownames(result_matrix) <- names(genelist)
  colnames(result_matrix) <- colnames(exp)
  
  for (i in names(genelist)) {
    gene_set = genelist[[i]]
    valid_genes = intersect(row_names, gene_set)
    
    if (length(valid_genes) == 0) {
      warning(paste("No valid genes for gene set", i, ". Setting enrichment score to 0."))
     es <-  rep(0,num_samples)
     
    }else if(length(valid_genes) == 1){
      gene_set_idx <- match(valid_genes,row_names)
      es <- R[gene_set_idx,]
      es <- es / num_genes
      
    }else{
        gene_set_idx <- match(valid_genes,row_names)
        es <- colSums(R[gene_set_idx,])
        es <- es / length(valid_genes)
        es <- es / num_genes
    }
    result_matrix[i, ] = es
  }
  return(result_matrix)
}
countToTpm <- function(counts, effLen, batch){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
TpmNorm <- function(data_exp,exp_len) {
  data_TPM <-apply(data_exp, 2, function(x) {
    countToTpm(x, exp_len$Length)
  })
  data_TPM <- as.data.frame(data_TPM)
  return(data_TPM)
}
