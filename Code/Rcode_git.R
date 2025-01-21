# lm model
library(car)
results <- list()
for (gene in rownames(combat_data)) {
  model_formula <- as.formula(paste("`", gene, "`", " ~ Age + Sex2",sep = ''))
  model <- lm(model_formula, data = combined_data)
  
  # b Slope
  coefficients <- coef(model)
  # p
  anova_result <- Anova(model, type = "II")
  
  # 
  diff_expr_info <- data.frame(Gene = gene,
                               Intercept = coefficients[1],
                               Slope_age = coefficients[2],
                               Slope_sex = coefficients[3],
                               p_value = anova_result['Age','Pr(>F)'],
                               p_value_sex = anova_result['Sex2','Pr(>F)'])
  
  results[[gene]] <- diff_expr_info
}

diff_expr_results <- do.call(rbind, results)
diff_expr_results$p_adjusted <- p.adjust(diff_expr_results$p_value, method = "BH")
diff_expr_results$p_adjusted_sex <- p.adjust(diff_expr_results$p_value_sex, method = "BH")

# losses 
run_loess <- function(gene, data) {
  model_formula <- as.formula(paste("`", gene, "`", " ~ Age",sep = ''))
  
  # Generate LOESS model and save to lo_data
  lo <- loess(model_formula, data, span = 0.75)
  
  # Generate predicted values
  lo_predict <- predict(lo, data.frame(Age = seq(age_min, age_max, 1))) %>%
    # lo_predict <- predict(lo, data.frame(age = seq(age_min, age_max, 0.0833))) %>%
    as.data.frame()
  colnames(lo_predict) <- gene
  
  return(lo_predict)
}

bulkRNAseq_run_loess <- function(expression_matrix, data) {
  results <- list()
  for (gene in rownames(expression_matrix)) {
    lo_predict <- run_loess(gene, data)
    results[[gene]] <- lo_predict
  }
  
  lo_predict_all <- do.call(cbind, results)
  rownames(lo_predict_all) <- seq(age_min, age_max, 1)
  return(lo_predict_all)
}

# DEswan
library(DEswan)
res.DEswan_10=DEswan(data.df =  t(combat_data),
                     qt = HC_inform$Age,
                     window.center = seq(10,90,1),
                     covariates = HC_inform$Sex2,
                     buckets.size = 20)

# nichenet
ligand_target_matrix = readRDS("/data5/laiwp/nichenet/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]# target genes in rows, ligands in columns

lr_network = readRDS("/data5/laiwp/nichenet/lr_network.rds")
head(lr_network)

weighted_networks = readRDS("/data5/laiwp/nichenet/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

head(weighted_networks$lr_sig)# interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
#########
CD4_Treg_sub <- CreateSeuratObject(counts = t(CD4_Treg_sub),min.cells = 0,min.features = 0)
CD4_Treg_sub <- AddMetaData(object = CD4_Treg_sub,metadata = CD4_Treg_sub_meta)

geneset_oi = subset(CD4_Treg_deg,Frail_logfoldchanges>0.75) %>% pull(Frail_names)
# geneset_oi = subset(CD4_Treg_deg,Frail_logfoldchanges<(-0.75)) %>% pull(Frail_names)
# geneset_oi = get_expressed_genes(c('Advanced Aged','Frail'), CD4_Treg_sub, pct = 0.50,assay_oi = 'RNA')
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

CD4_Treg_sub@active.ident <- as.factor(CD4_Treg_sub$Group)
expressed_genes_receiver = get_expressed_genes('Frail', CD4_Treg_sub, pct = 0.10,assay_oi = 'RNA')
# expressed_genes_receiver = get_expressed_genes('Advanced Aged', CD4_Treg_sub, pct = 0.10,assay_oi = 'RNA')
# expressed_genes_receiver = get_expressed_genes(c('Advanced Aged','Frail'), CD4_Treg_sub, pct = 0.10,assay_oi = 'RNA')
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# geneset_oi = rownames(subset(cDC_degs[cDC_degs$cluster %in% 'migDC_LAMP3',],p_val_adj<0.05 & avg_log2FC>0.25)) %>% .[. %in% rownames(ligand_target_matrix)]


expressed_ligands = ligands
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(15, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
# best_upstream_ligands <- c(best_upstream_ligands,'IFNG')
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 300) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.15)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes of CD4_Treg (Frail)", color = "#B2182B",legend_position = "right", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic"),axis.text=element_text(color='black')) + scale_fill_gradient2(low = "whitesmoke",  high = "#B2182B", breaks = c(0,0.0035,0.0070))
p_ligand_target_network
#
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = stats::dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", legend_position = "right",
                                                                                        x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "right", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson

# Cellchat
cellchat_Frail <- createCellChat(object = t(Frail), meta = Frail_meta, group.by = "Celltype")
cellchat_Advanced_Aged <- createCellChat(object = t(Aged), meta = Aged_meta, group.by = "Celltype")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling",'Cell-Cell Contact','ECM-Receptor'))
cellchat_Frail@DB <- CellChatDB.use
cellchat_Advanced_Aged@DB <- CellChatDB.use

cellchat_Frail <- subsetData(cellchat_Frail) 
cellchat_Advanced_Aged <- subsetData(cellchat_Advanced_Aged) 


future::plan("multiprocess", workers = 30)
#
cellchat_Frail <- identifyOverExpressedGenes(cellchat_Frail)
cellchat_Frail <- identifyOverExpressedInteractions(cellchat_Frail)

cellchat_Advanced_Aged <- identifyOverExpressedGenes(cellchat_Advanced_Aged)
cellchat_Advanced_Aged <- identifyOverExpressedInteractions(cellchat_Advanced_Aged)
#
cellchat_Frail <- computeCommunProb(cellchat_Frail)
cellchat_Frail <- filterCommunication(cellchat_Frail, min.cells = 2)
cellchat_Advanced_Aged <- computeCommunProb(cellchat_Advanced_Aged)
cellchat_Advanced_Aged <- filterCommunication(cellchat_Advanced_Aged, min.cells = 2)

#
cellchat_Frail <- computeCommunProbPathway(cellchat_Frail)
cellchat_Advanced_Aged <- computeCommunProbPathway(cellchat_Advanced_Aged)
#
cellchat_Frail <- aggregateNet(cellchat_Frail)
cellchat_Advanced_Aged <- aggregateNet(cellchat_Advanced_Aged)

cellchat_Frail <- netAnalysis_computeCentrality(cellchat_Frail, slot.name = "netP")
cellchat_Advanced_Aged <- netAnalysis_computeCentrality(cellchat_Advanced_Aged, slot.name = "netP")

# # Gini Clonality Shannon_index
# gini_coef <- function(TCR_df) {
#   sorted_df <- TCR_df[order(TCR_df$Fraction, decreasing = TRUE), ]
  
#   sorted_df$Fraction_2 <- cumsum(sorted_df$Fraction)
  
#   total_freq <- sum(sorted_df$Fraction)
  
#   gini_coef <- 1 - (sum((sorted_df$Fraction_2 / total_freq) * (2 * (1:(nrow(sorted_df)) - 1) - nrow(sorted_df))) / nrow(sorted_df))
  
#   return(gini_coef)
# }

# Shannon_index <- function(list_p){
#   sum=0
#   for (p in list_p){
#     sum <- p*log(p) + sum
#   }
#   H <- -sum
#   return(H)
# }

# clonality <- function(diversity,richness){
#   E<- diversity/log(richness)
#   C <- 1-E
#   return(C)
# }
