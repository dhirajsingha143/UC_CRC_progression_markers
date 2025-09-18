if(!requireNamespace("BiocManager", quitely = TRUE)){
  install.packages("BiocManager")
  library(BiocManager)
}
  
# install required packages
BiocManager::install("GEOquery")

# load packages
library(GEOquery)


# GEO data set load local downloaded files

gset1 <- getGEO(filename = "GSE47908_series_matrix.txt.gz", getGPL = FALSE)
gset2 <- getGEO(filename = "GSE20916_series_matrix.txt.gz", getGPL = FALSE)
gset3 <- getGEO(filename = "GSE10714_series_matrix.txt.gz", getGPL = FALSE)
gset4 <- getGEO(filename = "GSE37283_series_matrix.txt.gz", getGPL = FALSE)

# meta data

meta1 <- pData(gset1)
meta2 <- pData(gset2)
meta3 <- pData(gset3)
meta4 <- pData(gset4)

# expression data

exp1 <- exprs(gset1)
exp2 <- exprs(gset2)
exp3 <- exprs(gset3)
exp4 <- exprs(gset4)



