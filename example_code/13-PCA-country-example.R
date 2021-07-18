
library(tidyverse)
library(GGally)
library(cluster)
library(purrr)
library(dendextend)
library(ggcorrplot)
library(DataExplorer)
library(factoextra)
library(gridExtra)
library(vegan)


data <- read_csv("_GitBook/data/Country-data.csv")
head(data)
summary(data)

data_std <- decostand(data[, 2:10], method = "standardize")

data_pca <- rda(data_std)
data_pca

eigenvals(data_pca)

# species scores
scores(data_pca, choices = c(1, 2, 3), display = "sp")

# site scores
head(scores(data_pca, choices = c(1, 2, 3), display = "sites"))
