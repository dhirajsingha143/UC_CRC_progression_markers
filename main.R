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
BiocManager::install("ReactomePA")
BiocManager::install("ComplexHeatmap")
BiocManager::install("colorRamp2")
BiocManager::install("ggVennDiagram")

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
library(stringr)
library(ReactomePA)
library(ComplexHeatmap)
library(colorRamp2)
library(ggVennDiagram)


# GEO data set load local downloaded files

gset1 <- getGEO(filename = "GEO_Datasets/GSE47908_series_matrix.txt.gz", getGPL = FALSE) # RMA
gset2 <- getGEO(filename = "GEO_Datasets/GSE20916_series_matrix.txt.gz", getGPL = FALSE) # RMA
gset3 <- getGEO(filename = "GEO_Datasets/GSE37283_series_matrix.txt.gz", getGPL = FALSE) # RMA

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
  dplyr::select(1,32) %>% 
  dplyr::rename(condition = colnames(meta1)[32]) %>%
  dplyr::mutate(condition = gsub("left.sided.coltis","LSC",condition)) %>% 
  dplyr::mutate(condition = gsub("pancolitis","PC",condition)) %>% 
  dplyr::mutate(condition = gsub("ulcerative.colitis.associated.dysplasia","UCD",condition)) %>% 
  dplyr::mutate(condition = gsub("control","HC",condition)) %>% 
  dplyr::mutate(title = gsub(".\\d+$","",title)) %>% 
  dplyr::mutate(batch = "1")

group1 <- factor(meta1_clean$condition, levels = c("HC","LSC","PC","UCD"))

#meta2
meta2$`tissue:ch1` <- make.names(meta2$`tissue:ch1`)

meta2_clean <- meta2 %>% 
  dplyr::select(1,41) %>%
  dplyr::filter(!grepl("normal.colon",`tissue:ch1`)) %>%
  dplyr::filter(!grepl("colon.tumor",`tissue:ch1`)) %>%
  dplyr::rename(condition = `tissue:ch1`) %>% 
  dplyr::mutate(condition = gsub("adenoma","AD",condition)) %>% 
  dplyr::mutate(condition = gsub("adenocarcinoma","CRC",condition)) %>% 
  dplyr::mutate(title = gsub("_\\d+$","",title)) %>% 
  dplyr::mutate(batch = "2")

group2 <- factor(meta2_clean$condition, levels = c("AD","CRC"))

#meta3
meta3$source_name_ch1 <- make.names(meta3$source_name_ch1)

meta3_clean <- meta3 %>% 
  dplyr::select(1,8) %>% 
  dplyr::rename(condition = source_name_ch1) %>%
  dplyr::mutate(condition = gsub("quiescent.ulcerative.colitis","qUC",condition)) %>%
  dplyr::mutate(condition = gsub("normal.control","HC",condition)) %>%
  dplyr::mutate(condition = gsub("ulcerative.colitis.with.neoplasia","UCD",condition)) %>%
  dplyr::mutate(batch = "3")

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
custom_contrasts <- c(        # Step-wise Comparison for (local shift at each stage progression)
  "LSC_vs_HC = LSC - HC",     # best for diagnostic markers
  "PC_vs_LSC = PC - LSC",
  "qUC_vs_PC = qUC - PC",
  "UCD_vs_qUC = UCD - qUC",
  "AD_vs_UCD = AD - UCD",
  "CRC_vs_AD = CRC - AD",
  "CRC_vs_HC = CRC - HC"
)

custom_contrasts <- c(
  "LSC_vs_HC   = LSC - HC",   # Cumulative Contrast (exp deviation from HC as disease accumulates) 
  "PC_vs_HC    = PC  - HC",   # helps identify early vs late markers
  "UCD_vs_HC   = UCD - HC",
  "AD_vs_HC    = AD  - HC",
  "CRC_vs_HC   = CRC - HC"
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


--------------------------------------------------------------------------------
# Get Gene Symbols of DEGs-regulations

# HC contrast
# for up reg
genes_up_LSC <- deg_out$`LSC_vs_HC   = LSC - HC`$up$Gene
genes_up_PC <-  deg_out$`PC_vs_HC    = PC  - HC`$up$Gene
genes_up_UCD <- deg_out$`UCD_vs_HC   = UCD - HC`$up$Gene
genes_up_AD <- deg_out$`AD_vs_HC    = AD  - HC`$up$Gene
genes_up_CRC <- deg_out$`CRC_vs_HC   = CRC - HC`$up$Gene

# for down reg
genes_down_LSC <- deg_out$`LSC_vs_HC   = LSC - HC`$down$Gene
genes_down_PC <-  deg_out$`PC_vs_HC    = PC  - HC`$down$Gene
genes_down_UCD <- deg_out$`UCD_vs_HC   = UCD - HC`$down$Gene
genes_down_AD <- deg_out$`AD_vs_HC    = AD  - HC`$down$Gene
genes_down_CRC <- deg_out$`CRC_vs_HC   = CRC - HC`$down$Gene

--------------------------------------------------------------------------------
# Fetch annotations via BiomaRt Function ---

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

# for Up-reg
annot_up_LSC <- get_entrez(genes_up_LSC)
annot_up_PC <- get_entrez(genes_up_PC)
annot_up_UCD <- get_entrez(genes_up_UCD)
annot_up_AD <- get_entrez(genes_up_AD)
annot_up_CRC <- get_entrez(genes_up_CRC)

# for Down-reg
annot_down_LSC <- get_entrez(genes_down_LSC)
annot_down_PC <- get_entrez(genes_down_PC)
annot_down_UCD <- get_entrez(genes_down_UCD)
annot_down_AD <- get_entrez(genes_down_AD)
annot_down_CRC <- get_entrez(genes_down_CRC)

#---Enrichment Analysis (GO/KEGG/Reactome)----
# UP DEG table entreze id annotated
annot_up_LSC <- inner_join(annot_up_LSC, deg_out$`LSC_vs_HC   = LSC - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_PC <- inner_join(annot_up_PC, deg_out$`PC_vs_HC    = PC  - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_UCD <- inner_join(annot_up_UCD, deg_out$`UCD_vs_HC   = UCD - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_AD <- inner_join(annot_up_AD, deg_out$`AD_vs_HC    = AD  - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_CRC <- inner_join(annot_up_CRC, deg_out$`CRC_vs_HC   = CRC - HC`$up, by = c("external_gene_name" = "Gene"))

# DOWN DEG table entreze id annotated
annot_down_LSC <- inner_join(annot_down_LSC, deg_out$`LSC_vs_HC   = LSC - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_PC <- inner_join(annot_down_PC, deg_out$`PC_vs_HC    = PC  - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_UCD <- inner_join(annot_down_UCD, deg_out$`UCD_vs_HC   = UCD - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_AD <- inner_join(annot_down_AD, deg_out$`AD_vs_HC    = AD  - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_CRC <- inner_join(annot_down_CRC, deg_out$`CRC_vs_HC   = CRC - HC`$down, by = c("external_gene_name" = "Gene"))


source("scripts/enrichment.R")
#-------------------------------------------------------------------------------
# UP genes LSC
fc_values <- annot_up_LSC$logFC
names(fc_values) <- annot_up_LSC$external_gene_name   # make sure IDs match enrichment

res <- run_enrichment(
  gene_list  = annot_up_LSC$entrezgene_id,  # Entrez IDs of DEGs table
  fc_values  = fc_values,                   # gene vector of logFC values
  prefix     = "LSC_vs_HC_up",              # File naming prefix
  outdir     = "results/LSC_vs_HC/Enrichment/UP",      # dir name
  showCategory = 15,
  save_tables = TRUE
)

# DOWN genes LSC
fc_values <- annot_down_LSC$logFC
names(fc_values) <- annot_down_LSC$external_gene_name 

res <- run_enrichment(
  gene_list  = annot_down_LSC$entrezgene_id,
  fc_values  = fc_values,                     
  prefix     = "LSC_vs_HC_down",              
  outdir     = "results/LSC_vs_HC/Enrichment/DOWN",      
  showCategory = 15,
  save_tables = TRUE
)

# UP genes PC
fc_values <- annot_up_PC$logFC
names(fc_values) <- annot_up_PC$external_gene_name

res <- run_enrichment(
  gene_list  = annot_up_PC$entrezgene_id,  
  fc_values  = fc_values,                   
  prefix     = "PC_vs_HC_up",              
  outdir     = "results/PC_vs_HC/Enrichment/UP",      
  showCategory = 15,
  save_tables = TRUE
)

# DOWN genes PC
fc_values <- annot_down_PC$logFC
names(fc_values) <- annot_down_PC$external_gene_name

res <- run_enrichment(
  gene_list  = annot_down_PC$entrezgene_id,
  fc_values  = fc_values,                     
  prefix     = "PC_vs_HC_down",              
  outdir     = "results/PC_vs_HC/Enrichment/DOWN",      
  showCategory = 15,
  save_tables = TRUE
)

# UP genes UCD
fc_values <- annot_up_UCD$logFC
names(fc_values) <- annot_up_UCD$external_gene_name

res <- run_enrichment(
  gene_list  = annot_up_UCD$entrezgene_id,  
  fc_values  = fc_values,                   
  prefix     = "UCD_vs_HC_up",              
  outdir     = "results/UCD_vs_HC/Enrichment/UP",      
  showCategory = 15,
  save_tables = TRUE
)

# DOWN genes UCD
fc_values <- annot_down_UCD$logFC
names(fc_values) <- annot_down_UCD$external_gene_name

res <- run_enrichment(
  gene_list  = annot_down_UCD$entrezgene_id,
  fc_values  = fc_values,                     
  prefix     = "UCD_vs_HC_down",              
  outdir     = "results/UCD_vs_HC/Enrichment/DOWN",      
  showCategory = 15,
  save_tables = TRUE
)

# UP genes AD
fc_values <- annot_up_AD$logFC
names(fc_values) <- annot_up_AD$external_gene_name

res <- run_enrichment(
  gene_list  = annot_up_AD$entrezgene_id,  
  fc_values  = fc_values,                   
  prefix     = "AD_vs_HC_up",              
  outdir     = "results/AD_vs_HC/Enrichment/UP",      
  showCategory = 15,
  save_tables = TRUE
)

# DOWN genes AD
fc_values <- annot_down_AD$logFC
names(fc_values) <- annot_down_AD$external_gene_name

res <- run_enrichment(
  gene_list  = annot_down_AD$entrezgene_id,
  fc_values  = fc_values,                     
  prefix     = "AD_vs_HC_down",              
  outdir     = "results/AD_vs_HC/Enrichment/DOWN",      
  showCategory = 15,
  save_tables = TRUE
)

# UP genes CRC
fc_values <- annot_up_CRC$logFC
names(fc_values) <- annot_up_CRC$external_gene_name

res <- run_enrichment(
  gene_list  = annot_up_CRC$entrezgene_id,  
  fc_values  = fc_values,                   
  prefix     = "CRC_vs_HC_up",              
  outdir     = "results/CRC_vs_HC/Enrichment/UP",      
  showCategory = 15,
  save_tables = TRUE
)

# DOWN genes CRC
fc_values <- annot_down_CRC$logFC
names(fc_values) <- annot_down_CRC$external_gene_name

res <- run_enrichment(
  gene_list  = annot_down_CRC$entrezgene_id,
  fc_values  = fc_values,                     
  prefix     = "CRC_vs_HC_down",              
  outdir     = "results/CRC_vs_HC/Enrichment/DOWN",      
  showCategory = 15,
  save_tables = TRUE
)

#--------------------------------------------------------------------------------






#---Heat Map Generation for DEG Expressions----

top_genes <- unique(c(
  annot_up_LSC$external_gene_name,
  annot_down_LSC$external_gene_name,
  annot_up_PC$external_gene_name,
  annot_down_PC$external_gene_name,
  annot_up_UCD$external_gene_name,
  annot_down_UCD$external_gene_name,
  annot_up_AD$external_gene_name,
  annot_down_AD$external_gene_name,
  annot_up_CRC$external_gene_name,
  annot_down_CRC$external_gene_name
))

heat_expr_mat <- integrated_expr_data[top_genes, ]

#rownames(gene_row_group) <- gene_row_group$gene

# Arrange samples by disease then batch
sample_order <- integrated_meta_data %>%
  arrange(condition, batch) %>%
  rownames()

# Reorder expression matrix
heat_expr_mat <- heat_expr_mat[, integrated_meta_data$sample_id]

# Reorder sample annotation
annotation_col <- integrated_meta_data[integrated_meta_data$sample_id,]

# Column annotations
col_ha <- HeatmapAnnotation(
  Disease = annotation_col$condition,
  Batch = annotation_col$batch,
  col = list(
    Disease = c("HC"  = "gray",
                "LSC" = "green",
                "PC"  = "orange",
                "qUC" = "yellow",
                "UCD" = "red",
                "AD"  = "blue",
                "CRC" = "purple"
    ),
    Batch = c("1" = "salmon", "2" = "gold", "3" = "lightblue")
  ),
  annotation_height = unit(6, "mm")
)

# Row scaling (z-score)
heat_z <- t(scale(t(as.matrix(heat_expr_mat))))

Heatmap(
  heat_z,
  name = "Z-score",
  top_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_split = annotation_col$disease,  # optional: facet columns by disease
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  heatmap_legend_param = list(title = "Expression"),
  row_names_gp = gpar(fontsize = 4.5),        # Gene names (rows)
  column_names_gp = gpar(fontsize = 5)
)

# Box Plot of Mean DEG Expression per Sample by Disease

# Explicitly reorder disease levels
integrated_meta_data$condition <- factor(integrated_meta_data$condition,
                                       levels = c("HC", "LSC", "PC", "qUC", "UCD","AD","CRC"))

# Recalculate sample-wise means (if not already)
sample_means <- colMeans(heat_expr_mat)

# Plot
boxplot(sample_means ~ integrated_meta_data$condition,
        ylab = "Mean Z-score Expression",
        xlab = "Disease",
        main = "Mean DEG Expression per Sample by Disease",
        col = c("orange", "green", "purple", "cyan", "gray","lightblue","pink"),)

#--------------------------------------------------------------------------
# PPI Network Analysis
# Load the functions for PPI network analysis
source("scripts/PPI_functions.R")

#===PPI > LSC_vs_HC====
# Process Multiple Regulation per contrasts matrix

## up-regulation

# Building PPI network
network_LSC_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_LSC_vs_HC_up <- identify_hub_genes(
  g = network_LSC_vs_HC_up$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "up",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_LSC_vs_HC_up <- detect_ppi_modules(
  g = network_LSC_vs_HC_up$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "up", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation

# Building PPI network
network_LSC_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_LSC_vs_HC_down <- identify_hub_genes(
  g = network_LSC_vs_HC_down$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "down",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_LSC_vs_HC_down <- detect_ppi_modules(
  g = network_LSC_vs_HC_down$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "down", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/LSC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)


#===PPI > PC_vs_HC====
# Process Multiple Regulation per contrasts matrix

## up-regulation

# Building PPI network
network_PC_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_PC_vs_HC_up <- identify_hub_genes(
  g = network_PC_vs_HC_up$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "up",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_PC_vs_HC_up <- detect_ppi_modules(
  g = network_PC_vs_HC_up$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "up", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation _pending processing

# Building PPI network
network_PC_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_PC_vs_HC_down <- identify_hub_genes(
  g = network_PC_vs_HC_down$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "down",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_PC_vs_HC_down <- detect_ppi_modules(
  g = network_PC_vs_HC_down$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "down", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)


#===PPI > UCD_vs_HC====
# Process Multiple Regulation per contrasts matrix

## up-regulation

# Building PPI network
network_UCD_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_UCD_vs_HC_up <- identify_hub_genes(
  g = network_UCD_vs_HC_up$graph,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "up",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_UCD_vs_HC_up <- detect_ppi_modules(
  g = network_UCD_vs_HC_up$graph,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "up", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation

# Building PPI network
network_UCD_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_UCD_vs_HC_down <- identify_hub_genes(
  g = network_UCD_vs_HC_down$graph,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "down",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_UCD_vs_HC_down <- detect_ppi_modules(
  g = network_UCD_vs_HC_down$graph,
  condition = "UCD_vs_HC   = UCD - HC",
  regulation = "down", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)


#===PPI > AD_vs_HC====
# Process Multiple Regulation per contrasts matrix

## up-regulation

# Building PPI network
network_AD_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_AD_vs_HC_up <- identify_hub_genes(
  g = network_AD_vs_HC_up$graph,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "up",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_AD_vs_HC_up <- detect_ppi_modules(
  g = network_AD_vs_HC_up$graph,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "up", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation

# Building PPI network
network_AD_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_AD_vs_HC_down <- identify_hub_genes(
  g = network_AD_vs_HC_down$graph,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "down",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_AD_vs_HC_down <- detect_ppi_modules(
  g = network_AD_vs_HC_down$graph,
  condition = "AD_vs_HC    = AD  - HC",
  regulation = "down", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/AD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

#===PPI > CRC_vs_HC====

## up-regulation

# Building PPI network
network_CRC_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_CRC_vs_HC_up <- identify_hub_genes(
  g = network_CRC_vs_HC_up$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "up",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_CRC_vs_HC_up <- detect_ppi_modules(
  g = network_CRC_vs_HC_up$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "up", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation

# Building PPI network
network_CRC_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_CRC_vs_HC_down <- identify_hub_genes(
  g = network_AD_vs_HC_down$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "down",  # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_CRC_vs_HC_down <- detect_ppi_modules(
  g = network_AD_vs_HC_down$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "down", # Options: "up" , "down"
  save = TRUE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Define hub gene list

hub_genes <- list(
  LSC <- unique(c(hubs_LSC_vs_HC_up$hub_df$Gene_Symbol,hubs_LSC_vs_HC_down$hub_df$Gene_Symbol)),
  PC <- unique(c(hubs_PC_vs_HC_up$hub_df$Gene_Symbol,hubs_PC_vs_HC_down$hub_df$Gene_Symbol)),
  UCD <- unique(c(hubs_UCD_vs_HC_up$hub_df$Gene_Symbol,hubs_UCD_vs_HC_down$hub_df$Gene_Symbol)),
  AD <- unique(c(hubs_AD_vs_HC_up$hub_df$Gene_Symbol,hubs_AD_vs_HC_down$hub_df$Gene_Symbol)),
  CRC <- unique(c(hubs_CRC_vs_HC_up$hub_df$Gene_Symbol,hubs_CRC_vs_HC_down$hub_df$Gene_Symbol))
)

# hub genes venn diagram

# Prepare your list
venn_list <- list(
  LSC = LSC,PC = PC,UCD = UCD,AD = AD,CRC = CRC
)

# Plot
ggVennDiagram(venn_list, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "#FDE725FF", high = "#440154FF") +
  theme_void() +
  ggtitle("Hub Genes Across UC to CRC progression stages") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

#---Heat Map Generation for DEG Expressions----

top_genes <- unique(c(
  LSC,
  PC,
  UCD,
  AD,
  CRC
))

heat_expr_mat <- integrated_expr_data[top_genes, ]

# Reorder expression matrix
heat_expr_mat <- heat_expr_mat[, integrated_meta_data$sample_id]

# Reorder sample annotation
annotation_col <- integrated_meta_data[integrated_meta_data$sample_id,]

# Column annotations
col_ha <- HeatmapAnnotation(
  Disease = annotation_col$condition,
  Batch = annotation_col$batch,
  col = list(
    Disease = c("HC"  = "gray",
                "LSC" = "green",
                "PC"  = "orange",
                "qUC" = "yellow",
                "UCD" = "red",
                "AD"  = "blue",
                "CRC" = "purple"
    ),
    Batch = c("1" = "salmon", "2" = "gold", "3" = "lightblue")
  ),
  annotation_height = unit(6, "mm")
)

# Row scaling (z-score)
heat_z <- t(scale(t(as.matrix(heat_expr_mat))))

Heatmap(
  heat_z,
  name = "Z-score",
  top_annotation = col_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_split = annotation_col$disease,  # optional: facet columns by disease
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  heatmap_legend_param = list(title = "Expression"),
  row_names_gp = gpar(fontsize = 4.5),        # Gene names (rows)
  column_names_gp = gpar(fontsize = 5)
)



