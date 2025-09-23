# Install if not already installed
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")

library(WGCNA)
library(limma)
options(stringsAsFactors = FALSE)

allowWGCNAThreads() 

#Prepare expression data
datExpr <- t(integrated_expr_data) 

#Check for missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if(!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Sample clustering (detect outliers)
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


# Remove outliers if needed:
clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
keepSamples <- (clust==1)
datExpr <- datExpr[keepSamples, ]

#Choose soft-thresholding power
powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free topology fit index
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")

# Choose the power where R^2 > 0.8 or where it plateaus.

# Construct network and detect modules

softPower <- 10     # choose based on previous step
adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering on TOM-based dissimilarity", sub="", xlab="")

# Module identification
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)

# Convert numeric labels to colors
dynamicColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE)

#Relate modules to traits
#trait_data: samples × traits matrix (e.g., Disease status)
trait_data <- integrated_meta_data %>%
  dplyr::select(condition, batch)

# Convert categorical condition into numeric factors
trait_data$condition_num <- as.numeric(factor(trait_data$condition, 
                                              levels = c("HC", "LSC", "PC", "CRC")))

# Set rownames = sample IDs (must match datExpr rows)
rownames(trait_data) <- rownames(integrated_meta_data)

head(trait_data)


# Get module eigengenes
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

trait_data_sub <- trait_data[rownames(MEs), ]

# Double-check alignment
all(rownames(MEs) == rownames(trait_data_sub))

trait_numeric <- trait_data_sub %>% dplyr::select(where(is.numeric))

# Make sure rownames are aligned
trait_data_sub <- trait_data[rownames(MEs), ]

# Keep numeric traits only
trait_numeric <- trait_data_sub %>% dplyr::select(where(is.numeric))

# Correlation
moduleTraitCor <- cor(MEs, trait_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(MEs))

# Check first few correlations
head(moduleTraitCor)

# Heatmap of module-trait relationships

library(pheatmap)

pheatmap(moduleTraitCor,
         display_numbers = TRUE,
         cluster_cols = FALSE,   # only rows will be clustered
         main = "Module–Trait Relationships")

# Identify hub genes within key modules

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
geneTraitSignificance <- as.data.frame(cor(datExpr, trait_numeric, use = "p"))

# Assuming integrated_expr_data has genes as rownames
datExpr <- t(integrated_expr_data)

# Assign proper column names to datExpr
colnames(datExpr) <- rownames(integrated_expr_data)

# Now try again
blue_genes <- colnames(datExpr)[dynamicColors == "blue"]
head(blue_genes)
length(blue_genes)


#for the blue module:

module <- "blue"
moduleGenes <- colnames(datExpr)[dynamicColors == module]

kME <- geneModuleMembership[moduleGenes, paste0("ME", module)]
GS <- geneTraitSignificance[moduleGenes, 1]

hub_candidates <- data.frame(
  Gene = moduleGenes,
  kME = kME,
  GS = GS
)

# Apply cutoffs
hubs <- hub_candidates[abs(hub_candidates$kME) > 0.8 & abs(hub_candidates$GS) > 0.3, ]
head(hubs[order(-abs(hubs$kME)), ])
