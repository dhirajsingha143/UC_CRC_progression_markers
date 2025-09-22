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
# trait_data: samples × traits matrix (e.g., Disease status)
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes
moduleTraitCor <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Heatmap of module-trait relationships
library(pheatmap)
pheatmap(moduleTraitCor, display_numbers = TRUE, main="Module-Trait Relationships")

# Identify hub genes within key modules

# Module of interest
module <- "blue"  # example
moduleGenes <- colnames(datExpr)[dynamicColors == module]

# Module membership vs trait significance
geneModuleMembership <- cor(datExpr[,moduleGenes], MEs[,module], use = "p")
geneTraitSignificance <- cor(datExpr[,moduleGenes], trait_data$Disease, use = "p")



