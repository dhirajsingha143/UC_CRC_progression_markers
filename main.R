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
BiocManager::install("hpar")
BiocManager::install("pheatmap")

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
library(enrichplot)
library(ggupset)       # check for installation
library(stringr)
library(ReactomePA)
library(ComplexHeatmap)
library(colorRamp2)
library(ggVennDiagram)
library(hpar)
library(pheatmap)


# GEO data set load local downloaded files

gset1 <- getGEO(filename = "GEO_Datasets/GSE47908_series_matrix.txt.gz", getGPL = FALSE)  # RMA UC  - data set
gset2 <- getGEO(filename = "GEO_Datasets/GSE110224_series_matrix.txt.gz", getGPL = FALSE) # RMA CRC - data set

# meta data

meta1 <- pData(gset1)
meta2 <- pData(gset2)

# expression data

exp1 <- exprs(gset1)
exp2 <- exprs(gset2)

# cleaning meta data
#meta1
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

meta2_clean <- meta2 |> 
  dplyr::select(1,62) |> 
  dplyr::rename(condition = colnames(meta2)[62]) |> 
  dplyr::mutate(condition = gsub("normal.adjacent","HC",condition)) |> 
  dplyr::mutate(condition = gsub("primary.colorectal.adenocarcinoma","CRC",condition)) |> 
  dplyr::mutate(batch = "2")

group2 <- factor(meta2_clean$condition, levels = c("HC","CRC"))


# exp samples kept based on the metadata selection of samples

meta1_clean <- meta1_clean[order(group1),]
exp1 <- exp1[,rownames(meta1_clean)]

meta2_clean <- meta2_clean[order(group2),]
exp2 <- exp2[,rownames(meta2_clean)]


#Probe annotation using bio mart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)

grep("affy|U133_Plus_2", attributes$name, value = TRUE)

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


# Merge annotation with expression data
#expr1
expr1 <- as.data.frame(exp1)
expr1$probe_id <- rownames(expr1)

expr1_annotated <- inner_join(annot_expr1,expr1, by = "probe_id") |>
  filter(gene_symbol != "") |>
  group_by(gene_symbol) |>
  summarise(across(where(is.numeric),mean)) |>
  tibble::column_to_rownames(var = "gene_symbol")

#expr2
expr2 <- as.data.frame(exp2)
expr2$probe_id <- rownames(expr2)

expr2_annotated <- inner_join(annot_expr2,expr2,by = "probe_id") |>
  filter(gene_symbol != "") |>
  group_by(gene_symbol) |>
  summarize(across(where(is.numeric),mean)) |>
  filter(!if_any(where(is.numeric), is.na)) |>
  tibble::column_to_rownames(var = "gene_symbol")

#Merge by common gene and prepare final integration expr data

## Find common genes

common_genes <- intersect(rownames(expr1_annotated),rownames(expr2_annotated))

## Subset both data sets to include only common genes

integrated_expr_data <- cbind(
  expr1_annotated[common_genes, ],
  expr2_annotated[common_genes, ]
)

## integrated meta data
integrated_meta_data <- rbind(meta1_clean, meta2_clean)

group_meta_data <- factor(integrated_meta_data$condition, levels = c("HC", "LSC","PC","UCD","CRC"))

integrated_meta_data <- integrated_meta_data[order(group_meta_data),]

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
  levels = c("HC", "LSC","PC","UCD","CRC")
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "blue"
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "blue"
    )
  )

# Generate PCA plot

pca <- prcomp(t(integrated_expr_data), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC","PC","UCD","CRC"))
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
  levels = c("HC", "LSC","PC","UCD","CRC")
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
    )
  )

# Generate PCA plot
pca <- prcomp(t(normalized_expr_data), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC","PC","UCD","CRC"))
pca_df$batch <- integrated_meta_data$batch

normalized_exp_pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = batch)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = Group, color = Group), type = "norm", linetype = 2) +
  theme_minimal() +
  ggtitle("PCA of Normalized Expression Data") +
  theme(plot.title = element_text(hjust = 0.5))


#---Processing: Batch correction for integrated expr data------------------------------------------------------------

design <- model.matrix(~0 + condition, data = integrated_meta_data)

colnames(design) <- levels(factor(integrated_meta_data$condition))

batch_corrected_expr <- limma::removeBatchEffect(normalized_expr_data, batch = integrated_meta_data$batch, design = design)


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
  levels = c("HC", "LSC","PC","UCD","CRC")
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
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
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
    )
  )

# Generate PCA plot after batch correction

pca <- prcomp(t(batch_corrected_expr), scale. = TRUE)

pca_df <- as.data.frame(pca$x)
pca_df$Group <- factor(integrated_meta_data$condition, levels = c("HC", "LSC","PC","UCD","CRC"))
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
  "LSC_vs_HC   = LSC - HC",       # Cumulative Contrast (exp deviation from HC as disease accumulates) 
  "PC_vs_HC    = PC  - HC",       # helps identify early vs late markers
  "UCD_vs_HC    = UCD  - HC",   
  "CRC_vs_HC   = CRC - HC"
)


# Run DEG analysis with custom contrasts
deg_out <- run_deg_analysis(
  expr_data = batch_corrected_expr,
  meta_data = integrated_meta_data,
  condition_col = "condition",
  logfc_threshold = 1,
  pval_threshold = 0.05,
  contrasts_list = custom_contrasts
)


--------------------------------------------------------------------------------
# Get Gene Symbols of DEGs-regulations

# HC contrast
# for up reg
genes_up_LSC <- deg_out$`LSC_vs_HC   = LSC - HC`$up$Gene
genes_up_PC <-  deg_out$`PC_vs_HC    = PC  - HC`$up$Gene
genes_up_UCD <- deg_out$`UCD_vs_HC    = UCD  - HC`$up$Gene
genes_up_CRC <- deg_out$`CRC_vs_HC   = CRC - HC`$up$Gene

# for down reg
genes_down_LSC <- deg_out$`LSC_vs_HC   = LSC - HC`$down$Gene
genes_down_PC <-  deg_out$`PC_vs_HC    = PC  - HC`$down$Gene
genes_down_UCD <- deg_out$`UCD_vs_HC    = UCD  - HC`$down$Gene
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
annot_up_CRC <- get_entrez(genes_up_CRC)

# for Down-reg
annot_down_LSC <- get_entrez(genes_down_LSC)
annot_down_PC <- get_entrez(genes_down_PC)
annot_down_UCD <- get_entrez(genes_down_UCD)
annot_down_CRC <- get_entrez(genes_down_CRC)

#---Enrichment Analysis (GO/KEGG/Reactome)----
# UP DEG table entreze id annotated
annot_up_LSC <- inner_join(annot_up_LSC, deg_out$`LSC_vs_HC   = LSC - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_PC <- inner_join(annot_up_PC, deg_out$`PC_vs_HC    = PC  - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_UCD <- inner_join(annot_up_UCD, deg_out$`UCD_vs_HC    = UCD  - HC`$up, by = c("external_gene_name" = "Gene"))
annot_up_CRC <- inner_join(annot_up_CRC, deg_out$`CRC_vs_HC   = CRC - HC`$up, by = c("external_gene_name" = "Gene"))

# DOWN DEG table entreze id annotated
annot_down_LSC <- inner_join(annot_down_LSC, deg_out$`LSC_vs_HC   = LSC - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_PC <- inner_join(annot_down_PC, deg_out$`PC_vs_HC    = PC  - HC`$down, by = c("external_gene_name" = "Gene"))
annot_down_UCD <- inner_join(annot_down_UCD, deg_out$`UCD_vs_HC    = UCD  - HC`$down, by = c("external_gene_name" = "Gene"))
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
    Disease = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
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
  column_split = annotation_col$condition,  # optional: facet columns by disease
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  heatmap_legend_param = list(title = "Expression"),
  row_names_gp = gpar(fontsize = 4.5),        # Gene names (rows)
  column_names_gp = gpar(fontsize = 5)
)

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
  save = FALSE,
  out_dir = "results/LSC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_LSC_vs_HC_up <- identify_hub_genes(
  g = network_LSC_vs_HC_up$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "up",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/LSC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_LSC_vs_HC_up <- detect_ppi_modules(
  g = network_LSC_vs_HC_up$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "up", # Options: "up" , "down"
  save = FALSE,
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
  save = FALSE,
  out_dir = "results/LSC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_LSC_vs_HC_down <- identify_hub_genes(
  g = network_LSC_vs_HC_down$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "down",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/LSC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_LSC_vs_HC_down <- detect_ppi_modules(
  g = network_LSC_vs_HC_down$graph,
  condition = "LSC_vs_HC   = LSC - HC",
  regulation = "down", # Options: "up" , "down"
  save = FALSE,
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
  save = FALSE,
  out_dir = "results/PC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_PC_vs_HC_up <- identify_hub_genes(
  g = network_PC_vs_HC_up$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "up",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/PC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_PC_vs_HC_up <- detect_ppi_modules(
  g = network_PC_vs_HC_up$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "up", # Options: "up" , "down"
  save = FALSE,
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
  save = FALSE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_PC_vs_HC_down <- identify_hub_genes(
  g = network_PC_vs_HC_down$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "down",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_PC_vs_HC_down <- detect_ppi_modules(
  g = network_PC_vs_HC_down$graph,
  condition = "PC_vs_HC    = PC  - HC",
  regulation = "down", # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/PC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)


#===PPI > UCD_vs_HC====
# Process Multiple Regulation per contrasts matrix

## up-regulation

# Building PPI network
network_UCD_vs_HC_up <- build_ppi_network(
  deg_out = deg_out,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "up",  # Options: "up" , "down"
  top_n = 100,
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_UCD_vs_HC_up <- identify_hub_genes(
  g = network_UCD_vs_HC_up$graph,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "up",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_UCD_vs_HC_up <- detect_ppi_modules(
  g = network_UCD_vs_HC_up$graph,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "up", # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/UP",
  return_plot = TRUE
)

## down-regulation

# Building PPI network
network_UCD_vs_HC_down <- build_ppi_network(
  deg_out = deg_out,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "down",  # Options: "up" , "down"
  top_n = 100,
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_UCD_vs_HC_down <- identify_hub_genes(
  g = network_UCD_vs_HC_down$graph,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "down",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_UCD_vs_HC_down <- detect_ppi_modules(
  g = network_UCD_vs_HC_down$graph,
  condition = "UCD_vs_HC    = UCD  - HC",
  regulation = "down", # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/UCD_vs_HC/PPI/DOWN",
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
  save = FALSE,
  out_dir = "results/CRC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Identify hub genes
hubs_CRC_vs_HC_up <- identify_hub_genes(
  g = network_CRC_vs_HC_up$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "up",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/CRC_vs_HC/PPI/UP",
  return_plot = TRUE
)

# Detect modules
modules_CRC_vs_HC_up <- detect_ppi_modules(
  g = network_CRC_vs_HC_up$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "up", # Options: "up" , "down"
  save = FALSE,
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
  save = FALSE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Identify hub genes
hubs_CRC_vs_HC_down <- identify_hub_genes(
  g = network_CRC_vs_HC_down$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "down",  # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Detect modules
modules_CRC_vs_HC_down <- detect_ppi_modules(
  g = network_CRC_vs_HC_down$graph,
  condition = "CRC_vs_HC   = CRC - HC",
  regulation = "down", # Options: "up" , "down"
  save = FALSE,
  out_dir = "results/CRC_vs_HC/PPI/DOWN",
  return_plot = TRUE
)

# Define hub gene list

hub_genes <- list(
  LSC <- unique(c(hubs_LSC_vs_HC_up$hub_df$Gene_Symbol,hubs_LSC_vs_HC_down$hub_df$Gene_Symbol)),
  PC <- unique(c(hubs_PC_vs_HC_up$hub_df$Gene_Symbol,hubs_PC_vs_HC_down$hub_df$Gene_Symbol)),
  UCD <- unique(c(hubs_UCD_vs_HC_up$hub_df$Gene_Symbol,hubs_UCD_vs_HC_down$hub_df$Gene_Symbol)),
  CRC <- unique(c(hubs_CRC_vs_HC_up$hub_df$Gene_Symbol,hubs_CRC_vs_HC_down$hub_df$Gene_Symbol))
)

hub_genes <- list(
  LSC <- unique(c(hubs_LSC_vs_HC_up$hub_df$Gene_Symbol)),
  PC <- unique(c(hubs_PC_vs_HC_up$hub_df$Gene_Symbol)),
  UCD <- unique(c(hubs_UCD_vs_HC_up$hub_df$Gene_Symbol)),
  CRC <- unique(c(hubs_CRC_vs_HC_up$hub_df$Gene_Symbol))
)

hub_genes <- list(
  LSC <- unique(c(hubs_LSC_vs_HC_down$hub_df$Gene_Symbol)),
  PC <- unique(c(hubs_PC_vs_HC_down$hub_df$Gene_Symbol)),
  UCD <- unique(c(hubs_UCD_vs_HC_down$hub_df$Gene_Symbol)),
  CRC <- unique(c(hubs_CRC_vs_HC_down$hub_df$Gene_Symbol))
)


#---Heat Map Generation for DEG Expressions----

top_genes <- unique(c(
  LSC,
  PC,
  UCD,
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
    Disease = c(
      "HC"  = "gray",
      "LSC" = "green",
      "PC"  = "purple",
      "UCD"  = "orange",
      "CRC" = "yellow"
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
  column_split = annotation_col$condition,  # optional: facet columns by disease
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  heatmap_legend_param = list(title = "Expression"),
  row_names_gp = gpar(fontsize = 4.5),        # Gene names (rows)
  column_names_gp = gpar(fontsize = 5)
)


# Prepare your list
venn_list <- list(
  LSC = LSC,
  PC = PC,
  UCD = UCD,
  CRC = CRC
)

# Plot
ggVennDiagram(venn_list, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "#FDE725FF", high = "#440154FF") +
  theme_void() +
  ggtitle("Hub Genes Across UC to CRC progression") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# common genes between Stages

LSC_PC <- intersect(LSC, PC)
PC_CRC <- intersect(PC, CRC)

UC <- intersect(LSC_PC, PC_CRC)

UC_CRC <- intersect(UC, CRC) # early to CRC common significantly hub markers



# > LSC_PC
# [1] "IL1B"    "CXCL8"   "MMP9"    "CXCL1"   "LCN2"    "CXCL2"   "CXCR2"   "CXCL13"  "CXCL5"   "S100A12" "CXCL11"  "MMP1"    "CCL18"  
# [14] "CXCL3" 
# > LSC_CRC
# [1] "IL1B"   "CXCL8"  "CXCL1"  "LCN2"   "CXCL5"  "CXCL9"  "CXCL11" "MMP1" 
#> PC_CRC
# [1] "IL1B"   "CXCL8"  "LCN2"   "CXCL1"  "CXCL5"  "MMP3"   "MMP1"   "CXCL11"
----------------------------------------------------------------------------
# > UC_CRC
# [1] "IL1B"   "CXCL8"  "CXCL1"  "LCN2"   "CXCL5"  "CXCL11" "MMP1" 

# Enrichment of Common hub genes between Stages for shared pathway identification

# annotations of common genes to entrenz id

get_entrez <- function(gene_list) {
  getBM(
    attributes = c("external_gene_name", "entrezgene_id"),
    filters = "external_gene_name",
    values = gene_list,
    mart = ensembl
  ) %>% distinct()
}

UC_CRC_annot <- get_entrez(UC_CRC)

source("scripts/enrichment.R")
#-------------------------------------------------------------------------------
# UP HUB genes Common UC to CRC
fc_values <- annot_up_LSC$logFC
names(fc_values) <- annot_up_LSC$external_gene_name   # make sure IDs match enrichment

res <- run_enrichment(
  gene_list  = UC_CRC_annot$entrezgene_id,  # Entrez IDs of DEGs table
  fc_values  = NULL,                   # gene vector of logFC values
  prefix     = "UC_vs_CRC_HUB_up",              # File naming prefix
  outdir     = "results/UC_vs_CRC_HUB/Enrichment/UP",      # dir name
  showCategory = 15,
  save_tables = TRUE
)

#---Box Plot of genes Expression Across Disease Groups----

# Load

# Define pairwise comparisons (optional)

my_comparisons <- list(
  c("HC", "LSC"), c("HC", "PC"), c("HC", "UCD"), c("HC", "CRC")
)

gene_id <- "CCL18"
title <- NULL


expr_vec <- as.numeric(batch_corrected_expr[gene_id, ])
if (length(expr_vec) != nrow(integrated_meta_data)) stop("Expression and metadata mismatch")

# Create data frame for plotting
plot_df <- data.frame(
  expr = expr_vec,
  group = factor(integrated_meta_data$condition, levels = c("HC", "LSC","PC","UCD","CRC"))
)

# Generate boxplot using ggpubr
p <- ggboxplot(plot_df, x = "group", y = "expr", fill = "group",
              add =  "jitter", notch = TRUE) +
  labs(
    title = ifelse(is.null(title), paste("Expression of", gene_id), title),
    x = "Disease", y = "Batch Corrected Expression (log2)"
  ) +
  theme_bw(base_size = 13) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif", hide.ns = TRUE)


# Validation with HPA db

# Install required packages if not already installed


# Example hub gene list
hub_genes <- c("IL1B", "CXCL8", "CXCL1", "LCN2", "CXCL5", "CXCL11", "MMP1")  # replace with your genes

# Load all HPA pathology data
hpa_data <- allHparData()    # this contain all the data from HPA

# Load cancer data of HPA version 16.1
cancer_types <- hpaCancer16.1()

#Filter cancer
colorectal_cancer <- cancer_types %>% 
  filter(str_detect(Tumor, "colon|rectum|colorectal"))

# filter hub genes to HPA database
hub_in_crc <- colorectal_cancer %>%
  filter(Gene.name %in% hub_genes)

# summary

hub_summary <- hub_in_crc %>%
  dplyr::select(Gene.name, Level, Count.patients, Total.patients) %>%
  tidyr::pivot_wider(names_from = Level, values_from = Count.patients, values_fill = 0)

hub_summary

--------------------------------------------------------------------------------
# Visualization
ggplot(hub_in_crc, aes(x = Gene.name, y = Count.patients, fill = Level)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Not detected"="grey80",
                               "Low"="skyblue",
                               "Medium"="orange",
                               "High"="red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Gene Expression in Colorectal Cancer (HPA)",
       y = "Proportion of Patients", x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

heat_mat <- hub_summary %>%
  tibble::column_to_rownames("Gene.name") %>%
  dplyr::select(High, Medium, Low, `Not detected`) %>%
  as.matrix()

pheatmap(heat_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("grey95","skyblue","orange","red"))(50),
         main = "HPA Colorectal Cancer: Hub Genes (Patient Counts)",
         fontsize = 12)

ggplot(hub_in_crc, aes(x = Gene.name, y = Count.patients, fill = Level)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("Not detected"="grey80",
                               "Low"="skyblue",
                               "Medium"="orange",
                               "High"="red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Genes: Patient Counts by Expression Level (HPA CRC)",
       y = "Number of Patients", x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

install.packages("ggalluvial")
library(ggalluvial)

ggplot(hub_in_crc,
       aes(axis1 = Gene.name, axis2 = Level, y = Count.patients)) +
  geom_alluvium(aes(fill = Level), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("Not detected"="grey80",
                               "Low"="skyblue",
                               "Medium"="orange",
                               "High"="red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Genes – Patient Distribution (HPA CRC)",
       y = "Patient Counts", x = "")
--------------------------------------------------------------------------------

# Expression in normal tissues

# Load cancer data of HPA version 16.1
normal_types <- hpaNormalTissue16.1()

#Filter cancer
common_hub_genes <- normal_types %>% 
filter(str_detect(Tissue, "colon|rectum|colorectal"))

# filter hub genes to HPA database
normal_tissue <- common_hub_genes %>%
filter(Gene.name %in% hub_genes)

ggplot(normal_tissue, aes(x = Gene.name, fill = Level)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Not detected"="grey80",
                               "Low"="skyblue",
                               "Medium"="orange",
                               "High"="red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Gene Expression in Normal Colon/Rectum (HPA v16.1)",
       y = "Proportion of Cell Types", x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# summary

hub_summary <- hub_in_crc %>%
  dplyr::select(Gene.name, Level, Count.patients, Total.patients) %>%
  tidyr::pivot_wider(names_from = Level, values_from = Count.patients, values_fill = 0)

hub_summary


# rnaGtexTissue Normal RNA-seq, healthy colon/rectum samples

normal_rna_types <- rnaGtexTissue()

#Filter cancer
common_hub_genes_rna <- normal_rna_types %>% 
  filter(str_detect(Tissue, "colon|rectum|colorectal"))

# filter hub genes to HPA database
normal_rna_tissue <- common_hub_genes_rna %>%
  filter(Gene.name %in% hub_genes)


ggplot(normal_rna_tissue, aes(x = Gene.name, y = TPM, fill = Gene.name)) +
  geom_col(width = 0.7, alpha = 0.8, show.legend = FALSE) +
  geom_text(aes(label = round(TPM,1)), vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Gene Expression in Normal Colon (RNA-seq, TPM)",
       x = "Genes", y = "TPM")

ggplot(normal_rna_tissue, aes(x = reorder(Gene.name, TPM), y = TPM)) +
  geom_segment(aes(xend = Gene.name, y = 0, yend = TPM), color = "grey60") +
  geom_point(size = 6, color = "firebrick") +
  geom_text(aes(label = round(TPM,1)), hjust = -0.3, size = 4) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Genes Ranked by Expression in Colon (TPM)",
       x = "Genes", y = "TPM")

#comparison
--------------------------------------------------------------------------------
# Prepare cancer data (assuming hub_in_crc already exists)
  cancer_plot <- hub_in_crc %>%
  mutate(Source = "Cancer") %>%
  select(Gene.name, Level, Source, Count.patients)

# Prepare normal data
normal_plot <- normal_tissue %>%
  mutate(Source = "Normal",
         Count.patients = 1) %>%   # each row = one sample (cell type), so weight = 1
  select(Gene.name, Level, Source, Count.patients)

# Combine
combined_plot <- bind_rows(cancer_plot, normal_plot)

ggplot(combined_plot, aes(x = Gene.name, fill = Level, weight = Count.patients)) +
  geom_bar(position = "fill") +
  facet_wrap(~Source) +
  scale_fill_manual(values = c("Not detected"="grey80",
                               "Low"="skyblue",
                               "Medium"="orange",
                               "High"="red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Hub Genes Expression: Cancer vs Normal (HPA v16.1)",
       y = "Proportion", x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Construct the full path to the lock directory
lock_dir <- "G:/My Drive/AUUP/L319/Dhiraj/R_Projects_sync/UC_CRC_progression_markers/renv/library/windows/R-4.5/x86_64-w64-mingw32/00LOCK"

# Check if the directory exists before trying to delete it
if (dir.exists(lock_dir)) {
  # Use unlink() to remove the directory and its contents
  unlink(lock_dir, recursive = TRUE)
  print("Lock directory removed. You can now try to install the package again.")
} else {
  print("No lock directory found.")
}


