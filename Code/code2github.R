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

# Gini Clonality Shannon_index
gini_coef <- function(TCR_df) {
  sorted_df <- TCR_df[order(TCR_df$Fraction, decreasing = TRUE), ]
  
  sorted_df$Fraction_2 <- cumsum(sorted_df$Fraction)
  
  total_freq <- sum(sorted_df$Fraction)
  
  gini_coef <- 1 - (sum((sorted_df$Fraction_2 / total_freq) * (2 * (1:(nrow(sorted_df)) - 1) - nrow(sorted_df))) / nrow(sorted_df))
  
  return(gini_coef)
}

Shannon_index <- function(list_p){
  sum=0
  for (p in list_p){
    sum <- p*log(p) + sum
  }
  H <- -sum
  return(H)
}

clonality <- function(diversity,richness){
  E<- diversity/log(richness)
  C <- 1-E
  return(C)
}