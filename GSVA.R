# =========================
# 1. Load packages
# =========================
if(!requireNamespace("GSVA", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("GSVA")
if(!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
library(GSVA)
library(pheatmap)
library(dplyr)

# =========================
# 2. Prepare expression matrix
# =========================
# expr_data: normalized expression matrix (genes x samples)
# hub_genes: vector of your hub gene symbols

# Subset expression to hub genes (ensure all hub genes exist in your data)
hub_expr <- expr_data[rownames(expr_data) %in% hub_genes, ]

# Convert to numeric matrix (GSVA requires numeric)
hub_expr <- as.matrix(hub_expr)

# =========================
# 3. Prepare gene sets
# =========================
# Since these are hub genes, we can treat them as a single gene set
gene_sets <- list(HubGenes = rownames(hub_expr))

# =========================
# 4. Run GSVA
# =========================
gsva_results <- gsva(hub_expr, gene_sets, method="gsva", kcdf="Gaussian", mx.diff=TRUE)

# gsva_results: matrix of enrichment scores (gene set x samples)
head(gsva_results)

# =========================
# 5. Visualize results
# =========================
pheatmap(gsva_results, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         scale = "row",
         color = colorRampPalette(c("blue","white","red"))(50),
         main = "GSVA enrichment for hub genes")

# =========================
# 6. Optional: relate to traits
# =========================
# trait_data: samples x traits (e.g., DiseaseStatus: HC, UC, CRC)
# Correlation of GSVA scores with traits
if(exists("trait_data")) {
  cor_res <- cor(t(gsva_results), as.numeric(factor(trait_data$DiseaseStatus)))
  print(cor_res)
}
