# Quantitative Ecology (BCB743)
# Principal Component Analysis
# Author: AJ Smit
# Date: 9 July 2021


# Set-up ------------------------------------------------------------------

library(tidyverse)
library(vegan)


# The Doubs River data ----------------------------------------------------

env <- read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)
head(env)


# Do the PCA --------------------------------------------------------------

env_pca <- rda(env, scale = TRUE)
str(env_pca)
env_pca

env_pca$CA$eig

round(env_pca$CA$eig[1], 3)

sum(env_pca$CA$eig)

round(env_pca$CA$eig[1] / sum(env_pca$CA$eig) * 100, 1) # result in %

summary(env_pca)


# Graphical represenations of ordinations ---------------------------------

par(mfrow = c(1, 2))
biplot(env_pca, scaling = 1, main = "PCA scaling 1", choices = c(1, 2))
biplot(env_pca, scaling = 2, main = "PCA scaling 2", choices = c(1, 2))

# we need to load the function first from its R file:
source("Num_Ecol_R_book_ed1/cleanplot.pca.R")
cleanplot.pca(env_pca, scaling = 1)
cleanplot.pca(env_pca, scaling = 2)

biplot(env_pca, type = c("text", "points"), col = c("black", "black"))
ordisurf(env_pca ~ bod, env, add = TRUE, col = "turquoise", knots = 1)
ordisurf(env_pca ~ alt, env, add = TRUE, col = "salmon", knots = 1)
