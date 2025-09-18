if(!requireNamespace("BiocManager", quitely = TRUE)){
  install.packages("BiocManager")
  library(BiocManager)
}
  
# install required packages
BiocManager::install("GEOquery")
BiocManager::install("dplyr")
BiocManager::install("biomaRt")
BiocManager::install("reshape2")
BiocManager::install("ggplot2")
BiocManager::install("limma")
BiocManager::install("ggforce")
BiocManager::install("ggrepel")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

# load packages
library(GEOquery)
library(dplyr)
library(biomaRt)
library(tibble)
library(reshape2)
library(ggplot2)
library(limma)
library(ggforce)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)


# GEO data set load local downloaded files

gset1 <- getGEO(filename = "GSE47908_series_matrix.txt.gz", getGPL = FALSE) # RMA
gset2 <- getGEO(filename = "GSE20916_series_matrix.txt.gz", getGPL = FALSE) # RMA
gset3 <- getGEO(filename = "GSE37283_series_matrix.txt.gz", getGPL = FALSE) # RMA

# meta data

meta1 <- pData(gset1)
meta2 <- pData(gset2)
meta3 <- pData(gset3)

# expression data

exp1 <- exprs(gset1)
exp2 <- exprs(gset2)
exp3 <- exprs(gset3)

# cleaning meta data

meta1$`disease state:ch1` <- make.names(meta1$`disease state:ch1`)
meta1$title <- make.names(meta1$title)

meta1_clean <- meta1 %>% 
  select(1,32) %>% 
  rename(condition = `disease state:ch1`) %>%
  mutate(condition = gsub("left.sided.coltis","LSC",condition)) %>% 
  mutate(condition = gsub("pancolitis","PC",condition)) %>% 
  mutate(condition = gsub("ulcerative.colitis.associated.dysplasia","UCD",condition)) %>% 
  mutate(condition = gsub("control","HC",condition)) %>% 
  mutate(title = gsub(".\\d+$","",title)) %>% 
  mutate(batch = "1")

group1 <- factor(meta1_clean$condition, levels = c("HC","LSC","PC","UCD"))

#meta2
meta2$`tissue:ch1` <- make.names(meta2$`tissue:ch1`)

meta2_clean <- meta2 %>% 
  select(1,41) %>%
  filter(!grepl("normal.colon",`tissue:ch1`)) %>%
  filter(!grepl("colon.tumor",`tissue:ch1`)) %>%
  rename(condition = `tissue:ch1`) %>% 
  mutate(condition = gsub("adenoma","AD",condition)) %>% 
  mutate(condition = gsub("adenocarcinoma","CRC",condition)) %>% 
  mutate(title = gsub("_\\d+$","",title)) %>% 
  mutate(batch = "2")

group2 <- factor(meta2_clean$condition, levels = c("AD","CRC"))

#meta3
meta3$source_name_ch1 <- make.names(meta3$source_name_ch1)

meta3_clean <- meta3 %>% 
  select(1,8) %>% 
  rename(condition = source_name_ch1) %>%
  mutate(condition = gsub("quiescent.ulcerative.colitis","qUC",condition)) %>%
  mutate(condition = gsub("normal.control","HC",condition)) %>%
  mutate(condition = gsub("ulcerative.colitis.with.neoplasia","UCD",condition)) %>%
  mutate(batch = "3")

group3 <- factor(meta3_clean$condition, levels = c("HC","qUC","UCD"))

# exp samples kept based on the metadata selection of samples

meta1_clean <- meta1_clean[order(group1),]
exp1 <- exp1[,rownames(meta1_clean)]

meta2_clean <- meta2_clean[order(group2),]
exp2 <- exp2[,rownames(meta2_clean)]

meta3_clean <- meta3_clean[order(group3),]
exp3 <- exp3[,rownames(meta3_clean)]


#boxplot(exp1,outline=FALSE,las=2,main="GSE47908")
#boxplot(exp2,outline=FALSE,las=2,main="GSE20916")
#oxplot(exp3,outline=FALSE,las=2,main="GSE37283")


#Probe annotation using bio mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)

grep("affy|U133_Plus_PM", attributes$name, value = TRUE)

#Annotate expr1 data
probe_ids1 <- rownames(exp1) #[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

annot_expr1 <-getBM(
  attributes = c("affy_hg_u133_plus_2","external_gene_name"),
  filters = "affy_hg_u133_plus_2",
  values = probe_ids1,
  mart = ensembl
)

colnames(annot_expr1) <- c("probe_id","gene_symbol")

#Annotate expr2 data
probe_ids2 <- rownames(exp2) #[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array

annot_expr2 <- getBM(
  attribute = c("affy_hg_u133_plus_2", "external_gene_name"),
  filters = "affy_hg_u133_plus_2",
  values = probe_ids2,
  mart = ensembl
)

colnames(annot_expr2) <- c("probe_id","gene_symbol")

#Annotate expr3 data
probe_ids3 <- rownames(exp3) #[HT_HG-U133_Plus_PM] Affymetrix HT HG-U133+ PM Array Plate

annot_expr3 <- getBM(
  attribute = c("affy_ht_hg_u133_plus_pm", "external_gene_name"),
  filters = "affy_ht_hg_u133_plus_pm",
  values = probe_ids3,
  mart = ensembl
)

colnames(annot_expr3) <- c("probe_id","gene_symbol")

# Merge annotation with expression data
#expr1
expr1 <- as.data.frame(exp1)
expr1$probe_id <- rownames(expr1)

expr1_annotated <- inner_join(annot_expr1,expr1, by = "probe_id") |>
  filter(gene_symbol != "") |>
  group_by(gene_symbol) |>
  summarise(across(where(is.numeric),mean)) |>
  column_to_rownames(var = "gene_symbol")

#expr2
expr2 <- as.data.frame(exp2)
expr2$probe_id <- rownames(expr2)

expr2_annotated <- inner_join(annot_expr2,expr2,by = "probe_id") |>
  filter(gene_symbol != "") |>
  group_by(gene_symbol) |>
  summarize(across(where(is.numeric),mean)) |>
  filter(!if_any(where(is.numeric), is.na)) |>
  column_to_rownames(var = "gene_symbol")

#expr3
expr3 <- as.data.frame(exp3)
expr3$probe_id <- rownames(expr3)

expr3_annotated <- inner_join(annot_expr3,expr3,by = "probe_id") |>
  filter(gene_symbol != "") |>
  group_by(gene_symbol) |>
  summarize(across(where(is.numeric),mean)) |>
  filter(!if_any(where(is.numeric), is.na)) |>
  column_to_rownames(var = "gene_symbol")

#Merge by common gene and prepare final integration expr data

## Find common genes

com_platform <- intersect(rownames(expr1_annotated),rownames(expr3_annotated))

common_genes <- intersect(com_platform,rownames(expr2_annotated))

## Subset both data sets to include only common genes
integrated_expr_data <- cbind(
  expr1_annotated[common_genes, ],
  expr2_annotated[common_genes, ],
  expr3_annotated[common_genes, ]
)


## integrated meta data
integrated_meta_data <- rbind(meta1_clean, meta2_clean, meta3_clean)

group_meta_data <- factor(integrated_meta_data$condition, levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC"))

integrated_meta_data <- integrated_meta_data[order(group_meta_data),]
integrated_expr_data <- integrated_expr_data[,rownames(integrated_meta_data)]

# Reorder the columns of the expression data to match the row order of the metadata
integrated_expr_data <- integrated_expr_data[,rownames(integrated_meta_data)]

# verify checks
all(colnames(integrated_expr_data) == rownames(integrated_meta_data))

# Check common genes in integrated expression data & consistency with meta data

if(!all(rownames(integrated_expr_data) %in% common_genes)){
  stop("Not all genes in expression data are common between the two datasets.")
}else{
  message("Expression data successfully combined with common probe ID's")
}

if (!all(colnames(integrated_expr_data) %in% rownames(integrated_meta_data))) {
  stop("Sample names in expression data do not match meta data")
}else{
  message("Samples of integrated_expr_data match integrated_meta_data sample IDs")
}

# Box
boxplot(integrated_expr_data, outline=FALSE, las=2, main="Integrated Expression Data")

# Visualization of RAW integrated data sets

integrated_expr_data <- as.data.frame(integrated_expr_data)

#-----------Visualizations after pre-processing and integration--------------

integrated_expr_data$gene <- rownames(integrated_expr_data)

# Melt the expression data for ggplot
melted_data <- melt(integrated_expr_data, id.vars = "gene")

colnames(melted_data) <- c("gene", "sample_id", "expression")

# Restore row names (optional)
rownames(integrated_expr_data) <- integrated_expr_data$gene
integrated_expr_data$gene <- NULL

# Check sample ordering ## TESTING PURPOSE ONLY

integrated_meta_data$sample_id <- rownames(integrated_meta_data)

final_melted_data <- dplyr::left_join(melted_data, integrated_meta_data, by = "sample_id")

# Check sample ordering ## TESTING PURPOSE ONLY
col_order <- rownames(integrated_meta_data)

final_melted_data$sample_id <- factor(final_melted_data$sample_id, levels = col_order)

# Define condition order
final_melted_data$condition <- factor(
  final_melted_data$condition,
  levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC")
)


# Create the box plot
raw_exp_boxplot <- ggplot(final_melted_data, 
                          aes(x = sample_id, y = expression, fill = condition)) +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(limits = col_order) +  # sample order
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6)) +
  labs(x = "Sample", y = "Expression", title = "Raw expression data")

raw_exp_density_plot <- ggplot(final_melted_data, aes(x = expression, fill = condition)) +
  geom_density(alpha = 0.4) +
  labs(title = "Expression Distribution Raw expression data", x = "Expression", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  )

# Generate PCA plot

pca <- prcomp(t(integrated_expr_data), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC"))
pca_df$batch <- integrated_meta_data$batch

raw_exp_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = batch)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group, color = Group), type = "norm", linetype = 2) +
  theme_minimal() +
  ggtitle("PCA of Raw Expression Data") +
  theme(plot.title = element_text(hjust = 0.5))

#---------------------------------------------------------------------------------------------
# Normalization of expression data
normalized_expr_data <- normalizeBetweenArrays(integrated_expr_data, method = "quantile")

normalized_expr_data <- as.data.frame(normalized_expr_data)

#---Visualizations after pre-processing and integration------------------------------------------------

normalized_expr_data$gene <- rownames(normalized_expr_data)

# Melt the expression data for ggplot
melted_data <- melt(normalized_expr_data, id.vars = "gene")

colnames(melted_data) <- c("gene", "sample_id", "expression")

# Restore row names (optional)
rownames(normalized_expr_data) <- normalized_expr_data$gene
normalized_expr_data$gene <- NULL

# Check sample ordering ## TESTING PURPOSE ONLY

integrated_meta_data$sample_id <- rownames(integrated_meta_data)

final_melted_data <- dplyr::left_join(melted_data, integrated_meta_data, by = "sample_id")

# Check sample ordering ## TESTING PURPOSE ONLY
col_order <- rownames(integrated_meta_data)

final_melted_data$sample_id <- factor(final_melted_data$sample_id, levels = col_order)

# Define condition order
final_melted_data$condition <- factor(
  final_melted_data$condition,
  levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC")
)

# Create the box plot
normalized_exp_boxplot <- ggplot(final_melted_data, 
                                 aes(x = sample_id, y = expression, fill = condition)) +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(limits = col_order) +  # sample order
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6)) +
  labs(x = "Sample", y = "Expression", title = "Raw expression data")

normalized_exp_density_plot <- ggplot(final_melted_data, aes(x = expression, fill = condition)) +
  geom_density(alpha = 0.4) +
  labs(title = "Expression Distribution of Noramlized expression data", x = "Expression", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  )

# Generate PCA plot
pca <- prcomp(t(normalized_expr_data), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC"))
pca_df$batch <- integrated_meta_data$batch

normalized_exp_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = batch)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group, color = Group), type = "norm", linetype = 2) +
  theme_minimal() +
  ggtitle("PCA of Normalized Expression Data") +
  theme(plot.title = element_text(hjust = 0.5))


#---Processing: Batch correction for integrated expr data------------------------------------------------------------

batch_corrected_expr <- limma::removeBatchEffect(normalized_expr_data, batch = integrated_meta_data$batch)

#---Visualization after batch effect correction----

batch_corrected_expr <- as.data.frame(batch_corrected_expr)

#---Visualizations after pre-processing and integration------------------------------------------------

batch_corrected_expr$gene <- rownames(batch_corrected_expr)

# Melt the expression data for ggplot
melted_data <- melt(batch_corrected_expr)

colnames(melted_data) <- c("gene", "sample_id", "expression")

# Check sample ordering ## TESTING PURPOSE ONLY
batch_corrected_expr$gene <- NULL

final_melted_data <- dplyr::left_join(melted_data, integrated_meta_data, by = "sample_id")

# Check sample ordering ## TESTING PURPOSE ONLY
col_order <- rownames(integrated_meta_data)

final_melted_data$sample_id <- factor(final_melted_data$sample_id, levels = col_order)

# Define condition order
final_melted_data$condition <- factor(
  final_melted_data$condition,
  levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC")
)

# Create the box plot
batch_corrected_exp_boxplot <- ggplot(final_melted_data, 
                                      aes(x = sample_id, y = expression, fill = condition)) +
  geom_boxplot(outliers = FALSE) +
  scale_x_discrete(limits = col_order) +  # sample order
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 6)) +
  labs(x = "Sample", y = "Expression", title = "Raw expression data")

batch_corrected_exp_density_plot <- ggplot(final_melted_data, aes(x = expression, fill = condition)) +
  geom_density(alpha = 0.4) +
  labs(title = "Expression Distribution of Noramlized expression data", x = "Expression", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(
    values = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "orange",
      "qUC" = "yellow",
      "UCD" = "red",
      "AD"  = "blue",
      "CRC" = "purple"
    )
  )

# Generate PCA plot after batch correction

pca <- prcomp(t(batch_corrected_expr), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC", "PC", "qUC", "UCD", "AD", "CRC"))
pca_df$batch <- integrated_meta_data$batch


batch_corrected_exp_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = batch)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group, color = Group), type = "norm", linetype = 2) +
  theme_minimal() +
  ggtitle("PCA of Batch Corrected Expression Data") +
  theme(plot.title = element_text(hjust = 0.5))


#---Differential expression identification----
# Load the DEG functions
source("scripts/deg_functions.R")

# Contrast selection for differential expression analysis
custom_contrasts <- c(
  "LSC_vs_HC = LSC - HC", # inflammatory signatures
  "PC_vs_LSC = PC - LSC",
  "qUC_vs_PC = qUC - PC",
  "UCD_vs_qUC = UCD - qUC",
  "AD_vs_UCD = AD - UCD",
  "CRC_vs_AD = CRC - AD",
  "CRC_vs_HC = CRC - HC"
)

# Run DEG analysis with custom contrasts
deg_out <- run_deg_analysis(
  expr_data = batch_corrected_expr,
  meta_data = integrated_meta_data,
  condition_col = "condition",
  logfc_threshold = 1.5,
  pval_threshold = 0.05,
  contrasts_list = custom_contrasts
)

#---Retrieve full DEG tables into vectors----

# for full results
deg_LSC <- deg_out$`LSC_vs_HC = LSC - HC`$full_result
deg_PC <- deg_out$`PC_vs_LSC = PC - LSC`$full_result
deg_qUC <- deg_out$`qUC_vs_PC = qUC - PC`$full_result
deg_UCD <- deg_out$`UCD_vs_qUC = UCD - qUC`$full_result
deg_AD <- deg_out$`AD_vs_UCD = AD - UCD`$full_result
deg_CRC <- deg_out$`CRC_vs_AD = CRC - AD`$full_result
deg_CRC_HC <- deg_out$`CRC_vs_HC = CRC - HC`$full_result

# for Up-reg

up_LSC <- deg_out$`LSC_vs_HC = LSC - HC`$up
up_PC <- 
  
  
# for Down-reg
  
down_LSC <- deg_out$`LSC_vs_HC = LSC - HC`$down


# Get Gene Symbols

# for full results
genes_deg_LSC <- deg_LSC$Gene

# for up reg
genes_up_LSC <- up_LSC$Gene

# for down reg
genes_down_LSC <- down_LSC$Gene

# Fetch annotations via BioMart Function ---

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")


get_entrez <- function(gene_list) {
  getBM(
    attributes = c("external_gene_name", "entrezgene_id"),
    filters = "external_gene_name",
    values = gene_list,
    mart = ensembl
  ) %>% distinct()
}


#---Fetch Entrez IDs for each disease----

# for full results
annot_LSC <- get_entrez(genes_deg_LSC)

# for up-reg
annot_up_LSC <- get_entrez(genes_up_LSC)

#for down-reg
annot_down_LSC <- get_entrez(genes_down_LSC)


#---Merge annotations with DEG results----
# for full results
deg_out$`LSC_vs_HC = LSC - HC`$full_result <- left_join(deg_LSC,annot_LSC,by = c("Gene" = "external_gene_name"))

# for up-reg
deg_out$`LSC_vs_HC = LSC - HC`$up <- left_join(up_LSC,annot_up_LSC,by = c("Gene" = "external_gene_name"))

# for down-reg


#---Enrichment Analysis (Gene Ontology)----
# Perform GO enrichment analysis

# Prepare your named gene list (Entrez IDs per contrast)

gene_list_Up <- list(
  LSC_vs_HC <- deg_out$`LSC_vs_HC = LSC - HC`$up
)

##---Ontology enrichment results for all diseases----
gene_list_Up_entrez <- lapply(gene_list_Up, function(x) {
  mapped <- bitr(x, fromType = "Gene", toType = "entrezgene_id", OrgDb = org.Hs.eg.db)
  unique(mapped$ENTREZID)
})

#---GO > BP > UP----
go_compare_BP_up <- compareCluster(
  geneCluster = gene_list_Up,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP", # or "MF", "CC", or "ALL"
  pvalueCutoff = 0.05,
  readable = TRUE
)




