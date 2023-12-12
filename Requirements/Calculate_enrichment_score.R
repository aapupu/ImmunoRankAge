args = commandArgs(trailingOnly=TRUE)
source('function.R')
library(data.table)

if (length(args)!=2){
  
  print("Usage:file_path and norwayname must be provided")
  
}else{
  
  file_path <- args[1]
  norwayname <- args[2]
  
  rnaseq <- read.table(file_path,header = T,row.names = 1)
  
  selected_gene_length <- read.csv('Requirements/selected_gene_length.csv',header = T,row.names = 1)
  stopifnot(all(selected_gene_length$SYMBOL %in% rownames(rnaseq)))
  
  rnaseq <- rnaseq[selected_gene_length$SYMBOL,]
  if (norwayname=='count'){
    rnaseq <- TpmNorm(data_exp = rnaseq,exp_len = selected_gene_length)
    rnaseq <- log2(rnaseq + 1)
  }else if(norwayname=='tpm'){
    rnaseq <- log2(rnaseq + 1)
  }else{
    print('norwayname must be count or tpm')
  }

  enrich_select_df <- fread('Requirements/selected_geneset.csv',header = T)[,2:237]
  enrich_genelist <- list()
  for (i in colnames(enrich_select_df)){
    enrich_genelist[[i]] <- as.character(na.omit(enrich_select_df[[i]]))
  }
  
  es_TPM <- rankEnrich(as.matrix(rnaseq),genelist = enrich_genelist)
  write.csv(es_TPM,file = 'enrich_score.csv')
  
  print('Enrichment scoring is complete.')
}