# Quantitative Ecology (BCB743)
# Principal Coordinate Analysis
# Author: AJ Smit
# Date: 15 July 2021


# Set-up ------------------------------------------------------------------

library(tidyverse)
library(vegan)


# The Doubs River species data --------------------------------------------

spe <- read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv")
spe <- dplyr::select(spe, -1)
spe <- dplyr::slice(spe, -8)
head(spe, 8)


# Do the PCoA -------------------------------------------------------------

# optional: calculate the Bray-Curtis dissimilarity, which is
# appropriate for abundance data

spe_bray <- vegdist(spe)
spe_bray

# we can only specify the input via a formula interface
spe_pcoa <- capscale(spe_bray ~ 1)
spe_pcoa

summary(spe_pcoa)

# I advocate providing a raw species table to capscale() to
# retain the species information
spe_pcoa <- capscale(spe ~ 1, distance = "bray")
spe_pcoa

summary(spe_pcoa)

# the percentage inertia explained by the first three axes is
round(sum(spe_pcoa$CA$eig[1:3]) / sum(spe_pcoa$CA$eig) * 100, 2)


# Ordination diagrams -----------------------------------------------------

par(mfrow = c(2, 2))
plot(spe_pcoa, scaling = 1, main = "PCoA fish abundances - biplot scaling 1")
plot(spe_pcoa, scaling = 2, main = "PCoA fish abundances - biplot scaling 2")
plot(spe_pcoa, choices = c(1, 3), scaling = 1, main = "PCoA fish abundances - biplot scaling 1")
plot(spe_pcoa, choices = c(1, 3), scaling = 2, main = "PCoA fish abundances - biplot scaling 2")

# improve the plots
pl1 <- ordiplot(spe_pcoa, type = "none", scaling = 1, main = "PCoA fish abundances - biplot scaling 1")
points(pl1, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl1, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl1, "species", col = "blue4", cex = 0.9)
text(pl1, "sites", col = "red4", cex = 0.9)

pl2 <- ordiplot(spe_pcoa, type = "none", scaling = 2, main = "PCoA fish abundances - biplot scaling 2")
points(pl2, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl2, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl2, "species", col = "blue4", cex = 0.9)
text(pl2, "sites", col = "red4", cex = 0.9)

# we can also fit response surfaces using ordisurf()
par(mar = c(4, 4, 0.9, 0.5) + .1, mfrow = c(2, 2))
with(spe, tmp <- ordisurf(spe_pcoa ~ Satr, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Satr"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_pcoa ~ Scer, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Scer"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_pcoa ~ Teso, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Teso"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_pcoa ~ Cogo, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Cogo"))
abline(h = 0, v = 0, lty = 3)

env <- read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)
env <- dplyr::slice(env, -8)

(spe_pcoa_env <- envfit(spe_pcoa, env, scaling = 2))
plot(spe_pcoa_env, col = "grey40")
plot(spe_pcoa_env, p.max = 0.05, col = "red")
