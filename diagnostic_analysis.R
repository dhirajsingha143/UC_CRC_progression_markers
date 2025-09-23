#install
install.packages("pROC")

library(pROC)
library(ggplot2)
library(dplyr)

# Define hub genes
genes <- c("IL1B", "CXCL8", "CXCL1", "LCN2", "CXCL5", "CXCL11", "MMP1")

# All pairwise comparisons
comparisons <- combn(unique(integrated_meta_data$condition), 2, simplify=FALSE)

roc_results <- list()
auc_table <- data.frame(Gene=character(), Comparison=character(), 
                        AUC=numeric(), CI_low=numeric(), CI_high=numeric())

# Loop over comparisons
for (comp in comparisons) {
  group1 <- comp[1]
  group2 <- comp[2]
  comp_name <- paste(group1, "vs", group2, sep="_")
  
  # Subset metadata
  keep_samples <- integrated_meta_data$condition %in% c(group1, group2)
  sub_meta <- integrated_meta_data[keep_samples, ]
  sub_expr <- integrated_expr_data[, keep_samples]
  
  # Binary labels: group1 = 0, group2 = 1
  labels <- factor(ifelse(sub_meta$condition == group1, 0, 1))
  
  for (gene in genes) {
    if (!(gene %in% rownames(sub_expr))) {
      message("Skipping ", gene, " (not found in expression data)")
      next
    }
    
    gene_expr <- as.numeric(sub_expr[gene, ])
    
    if (length(gene_expr) != length(labels)) {
      message("Skipping ", gene, " (length mismatch)")
      next
    }
    
    roc_obj <- roc(labels, gene_expr, ci=TRUE)
    roc_results[[paste(gene, comp_name, sep="_")]] <- roc_obj
    
    auc_table <- rbind(auc_table, data.frame(
      Gene = gene,
      Comparison = comp_name,
      AUC = as.numeric(auc(roc_obj)),
      CI_low = roc_obj$ci[1],
      CI_high = roc_obj$ci[3]
    ))
  }
}

# Inspect results
auc_table <- auc_table %>% arrange(desc(AUC))
print(auc_table)

# ---- ROC PLOTTING for one comparison ----
plot_comparison <- "HC_vs_PC"  # change as needed

roc_df <- data.frame()
for (gene in genes) {
  key <- paste(gene, plot_comparison, sep="_")
  if (!is.null(roc_results[[key]])) {
    roc_obj <- roc_results[[key]]
    coords_df <- coords(roc_obj, "all", ret=c("specificity","sensitivity"))
    df <- data.frame(
      Specificity = coords_df$specificity,
      Sensitivity = coords_df$sensitivity,
      Gene = gene,
      AUC = round(auc(roc_obj), 2)
    )
    roc_df <- rbind(roc_df, df)
  }
}

ggplot(roc_df, aes(x=1-Specificity, y=Sensitivity, color=Gene)) +
  geom_line(size=1) +
  geom_abline(linetype="dashed", color="gray") +
  theme_minimal(base_size=14) +
  labs(
    title = paste("ROC Curves for Hub Genes:", plot_comparison),
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  scale_color_brewer(palette="Set1",
                     labels = paste0(unique(roc_df$Gene),
                                     " (AUC=", unique(roc_df$AUC), ")"))


# Pairwise comparisons ROC 
## binary between-group performance ( HC vs disease groups )

library(pROC)
library(data.table)  # for dcast

genes <- c("IL1B", "CXCL8", "CXCL1", "LCN2", "CXCL5", "CXCL11", "MMP1")
genes <- genes[genes %in% rownames(integrated_expr_data)]  # keep only valid ones

comparisons <- list(
  "HC_vs_CRC" = c("HC", "LSC"),
  "HC_vs_LSC" = c("HC", "PC"),
  "HC_vs_PC"  = c("HC", "CRC")
)

roc_results <- list()
auc_table <- data.frame(Gene=character(), Comparison=character(),
                        AUC=numeric(), CI_low=numeric(), CI_high=numeric())

for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  subset_idx <- integrated_meta_data$condition %in% groups
  labels <- factor(integrated_meta_data$condition[subset_idx])
  
  for (gene in genes) {
    gene_expr <- as.numeric(integrated_expr_data[gene, subset_idx])
    
    roc_obj <- roc(labels, gene_expr, ci=TRUE)
    
    auc_table <- rbind(auc_table, data.frame(
      Gene = gene,
      Comparison = comp_name,
      AUC = as.numeric(auc(roc_obj)),
      CI_low = roc_obj$ci[1],
      CI_high = roc_obj$ci[3]
    ))
    
    roc_results[[paste0(gene, "_", comp_name)]] <- roc_obj
  }
}


# Convert auc_table into a data.table before casting
setDT(auc_table)

auc_matrix <- dcast(auc_table, Gene ~ Comparison, value.var="AUC")

print(auc_matrix)

# Plot heatmap
library(pheatmap)

# Convert to matrix (removing gene column)
auc_mat <- as.matrix(as.data.frame(auc_matrix[ , -1]))
rownames(auc_mat) <- auc_matrix$Gene

# Heatmap
pheatmap(auc_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "AUC Heatmap of Hub Genes Across Comparisons",
         display_numbers = TRUE,  # shows AUC values inside cells
         number_format = "%.2f")

# Gene Panel ROC curve AUC value
library(pROC)

# Candidate panel: IL1B, MS4A12, CXCL11, CXCL5
panel_genes <- c("IL1B", "CXCL8", "CXCL1", "LCN2", "CXCL5", "CXCL11", "MMP1")

# Extract expression matrix for panel genes (samples x genes)
expr_panel <- as.data.frame(t(integrated_expr_data[panel_genes, ]))

# Add labels
expr_panel$condition <- factor(ifelse(integrated_meta_data$condition == "LSC", 1, 0))  
# 1 = CRC, 0 = non-CRC (can also do pairwise comparisons separately)

# Logistic regression model
panel_model <- glm(condition ~ ., data=expr_panel, family=binomial)

# Predict probabilities
probabilities <- predict(panel_model, type="response")

# ROC curve for combined panel
roc_panel <- roc(expr_panel$condition, probabilities, ci=TRUE)

# Plot
plot(roc_panel, col="darkblue", lwd=2, main="ROC for Multi-Gene Panel")
legend("bottomright", legend=paste0("Panel AUC = ", round(auc(roc_panel), 2)), 
       col="darkblue", lwd=2)


