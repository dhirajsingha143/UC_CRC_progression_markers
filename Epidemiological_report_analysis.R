library(ggplot2)
library(rworldmap) # installing from Github insted

install.packages("devtools")
library(devtools)

install_github('AndySouth/rworldmap', build_vignettes=TRUE)

library(dplyr)

# Load your data (assuming your CSV has 'country' and 'prevalence' columns)

# Incidence data

IBD <- read.csv("GEO_Datasets/Incidence_IBD.csv")
CD <- read.csv("GEO_Datasets/Incidence_CD.csv")
UC <- read.csv("GEO_Datasets/Incidence_UC.csv")
CRC <- read.csv("GEO_Datasets/Incidence_CRC.csv")
CRC <- CRC[-186,]

avg_IBD <- mean(as.numeric(IBD$incidence_rate), na.rm = TRUE)
avg_CD <- mean(as.numeric(CD$incidence_rate), na.rm = TRUE)
avg_UC <- mean(as.numeric(UC$incidence_rate), na.rm = TRUE)
avg_CRC <- mean(as.numeric(CRC$ASR..World.), na.rm = TRUE)

# input data from cancer today website

avg_CRC <- 18.4

df <- data.frame(
  Condition = c("IBD", "CD", "UC", "CRC"),
  Average_Incidence = c(avg_IBD, avg_CD, avg_UC, avg_CRC)
)


