# Quantitative Ecology (BCB743)
# Correspondence Analysis
# Author: AJ Smit
# Date: 14 July 2021


# Set-up ------------------------------------------------------------------

library(tidyverse)
library(vegan)


# The Doubs River species data --------------------------------------------

spe <- read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv")
spe <- dplyr::select(spe, -1)
head(spe, 8)


# Do the CA ---------------------------------------------------------------

spe_ca <- cca(spe) # !!! Why?

# one of the rows sums to 0
# which?
apply(spe, 1, sum)

# remove offending row
spe <- spe[rowSums(spe) > 0, ]
head(spe, 8)

# do the CA again
spe_ca <- cca(spe)
spe_ca

summary(spe_ca)

# sum of all eigenvalues
# is the total inertia
round(sum(spe_ca$CA$eig), 5)

# the first eigenvalue
round(spe_ca$CA$eig[1], 5)

# sum of CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]), 5)

# proportion explained by CA1 and CA2
round(sum(spe_ca$CA$eig[1:2]) / sum(spe_ca$CA$eig) * 100, 2)


# Ordination diagrams -----------------------------------------------------

par(mfrow = c(1, 2))
plot(spe_ca, scaling = 1, main = "CA fish abundances - biplot scaling 1")
# choices = c(1, 2)
plot(spe_ca, scaling = 2, main = "CA fish abundances - biplot scaling 2")

# Scaling 1: This scaling emphasises relationships between rows accurately in low-dimensional ordination space. Distances among objects (samples or sites) in the biplot are approximations of their 𝜒2  distances in multidimensional space. Objects found near a point representing a species are likely to contain a high contribution of that species.

# Scaling 2: This scaling emphasises relationships between columns accurately in low-dimensional ordination space. Distances among objects (samples or sites) in the biplot are not approximations of their 𝜒2  distances in multidimensional space, but the distances among species are. Species positioned close to the point representing an object (a sample or site) are more likely to be found in that object or to have higher frequency there.

# advanced plots
require('viridis')
palette(viridis(8))
par(mar = c(4, 4, 0.9, 0.5) + .1, mfrow = c(2, 2))
with(spe, tmp <- ordisurf(spe_ca ~ Satr, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Satr"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Scer, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Scer"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Teso, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Teso"))
abline(h = 0, v = 0, lty = 3)
with(spe, tmp <- ordisurf(spe_ca ~ Cogo, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Cogo"))
abline(h = 0, v = 0, lty = 3)

# a posteriori projection of environmental variables in a CA
env <- read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv")
env <- dplyr::select(env, -1)

# we removed the 8th row in spe, so do it here too
env <- dplyr::slice(env, -8)

# the last plot produced (CA scaling 2) must be active
# scaling 2 is default
(spe_ca_env <- envfit(spe_ca, env, scaling = 2))
plot(spe_ca_env, col = "grey40")
plot(spe_ca_env, p.max = 0.05, col = "red") # plot significant variables with a different colour
