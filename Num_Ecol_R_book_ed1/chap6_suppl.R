################################################################################
### CHAPTER 6 - CANONICAL ORDINATION
###
### Updated by AJS -- focus on CCA
###
### Online supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 3 February 2017
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(MASS)
# library(ellipse)
# library(FactoMineR)
library(tidyverse)

source("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/Num_Ecol_R_book_ed1/text_mult.R")

# Import the data from CSV files
# (files must be in the working directory)
spe <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1))
env <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
spa <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpa.csv", row.names = 1))
# Remove empty site 8
env <- env[-8, ]
spe <- spe[-8, ]
spa <- spa[-8, ]

# Set aside the variable 'dfs' (distance from the source) for later use
dfs <- env[, 1]

# Remove the 'dfs' variable from the env dataset
env <- decostand(env[, -1], method = "standardize")

# Create two subsets of explanatory variables
# Physiography (upstream-downstream gradient)
envtopo <- env[, c(1:2)]
names(envtopo)
# Water quality
envchem <- env[, c(4:9)]
names(envchem)

# Hellinger-transform the species dataset
# Suited to species abundance data. It gives low weights to species
# with low counts and many zeros. The transformation comprises dividing
# each value in a data matrix by its row sum, and taking the square root
# of the quotient.
spe.hel <- as_tibble(decostand(spe, "hellinger"))
spe.hel

# RDA of Hellinger-transformed fish data, constrained
# by all env. vars. in the environmental dataset (sans d.f.s.)
spe.rda <- rda(spe.hel ~ ., env)

# Is the fit significant?
anova(spe.rda) # ... yes!

# The adjusted R2 --- the variance explained by the constrained axes:
RsquareAdj(spe.rda)$adj.r.squared

# Variance explained by full model:
sum(spe.rda$CCA$eig) / spe.rda$tot.chi * 100

# Which axes are significant?
anova(spe.rda, by = "axis")

# Which terms are significant?
anova(spe.rda, by = "terms")

# Make some plots:
plot(spe.rda, scaling = 1, display = c("sp", "lc", "cn"),
     main = "Site scaling (scaling 1)")
plot(spe.rda, scaling = 2, display = c("sp", "lc", "cn"),
     main = "Species scaling (scaling 2)")

# Another way - build from scratch:
## better control -- remember to set scaling etc. identically
plot(spe.rda, type = "n", scaling = "sites", choices = 1:2)
text(spe.rda, dis = "cn", scaling = "sites", choices = 1:2)
points(spe.rda, pch = 21, col = "salmon", bg = "grey90", cex = 1.2, scaling = "sites", choices = 1:2)
text(spe.rda, "species", col = "blue", cex = 0.8, scaling = "sites", choices = 1:2)
text(spe.rda, "sites", col = "red", cex = 0.4, scaling = "sites", choices = 1:2)

# Make a ggplot2 ordination plot:
# The significant terms:
spe.rda.axis.test <- anova(spe.rda, by = "terms", parallel = 4)
spe.rda.ax <- which(spe.rda.axis.test[, 4] < 0.05)
spe.rda.sign.ax <- colnames(env[,spe.rda.ax])

# The biplot scores for constraining variables:
scores(spe.rda, display = "bp", choices = c(1:2))

spe.rda.scrs <- scores(spe.rda, display = c("sp","wa","lc","bp","cn"))
spe.rda.df_sites <- data.frame(spe.rda.scrs$constraints)

multiplier <- ordiArrowMul(spe.rda.scrs$biplot, fill = 0.45)
spe.rda.bp <- spe.rda.scrs$biplot * multiplier
spe.rda.bp <- as.data.frame(spe.rda.bp)
spe.rda.bp$labels <- rownames(spe.rda.bp)
colnames(spe.rda.bp) <- c("x", "y", "labels")
spe.rda.bp.sign <- spe.rda.bp[spe.rda.bp$labels %in% spe.rda.sign.ax,]

spe.rda.text <- text.mult(spe.rda.scrs$biplot, fill = 0.45)
spe.rda.text <- as.data.frame(spe.rda.text)
spe.rda.text$labels <- rownames(spe.rda.text)
colnames(spe.rda.text) <- c("x", "y", "labels")
spe.rda.text.sign <- spe.rda.text[spe.rda.text$labels %in% spe.rda.sign.ax,]

spe.rda.p <- ggplot(data = spe.rda.df_sites, aes(RDA1, RDA2)) +
  geom_point(size = 2.0, shape = 4) +
  geom_segment(data = spe.rda.bp.sign,
               aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "red", alpha = 1, size = 0.4) +
  geom_text(data = as.data.frame(spe.rda.text.sign),
            aes(x, y, label = rownames(spe.rda.text.sign)),
            color = "black") +
  xlab("RDA1") + ylab("RDA2") +
  ggtitle("Significant environmental drivers") +
  theme_grey() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(.80, .75),
        aspect.ratio = 0.8)


# Do a canonical correspondence analysis (CCA)
(spe.cca <- cca(spe.hel ~ ., env))
# spe.cca <- cca(spe.hel ~ alt + slo + flo + pH + har + pho + nit + amm, env)
summary(spe.cca)

# The adjusted R2 --- the variance explained by the constrained axes
(spe.cca.R2a <- RsquareAdj(spe.cca)$adj.r.squared)

# Variance explained by full model
sum(spe.cca$CCA$eig) / spe.cca$tot.chi * 100

# Global test of the CCA result..., i.e.,
# is the fit significant?
# ... 'parallel = 4' for multicore processing
# ... this needs to be adjusted for you own computer
anova(spe.cca, parallel = 4) # ... yes!

# Tests of all canonical axes
anova(spe.cca, by = "axis", parallel = 4)
(spe.cca.axis.test <- anova(spe.cca, by = "term", parallel = 4))

# Select only sign. terms
(spe.cca.ax <- which(spe.cca.axis.test[, 4] < 0.05))
env2 <- env[, spe.cca.ax]

# Make a reduced model with only the significant terms
(spe.cca2 <- cca(spe.hel ~., env2))
summary(spe.cca2)
anova(spe.cca2, parallel = 4)
anova(spe.cca2, by = "axis", parallel = 4)
anova(spe.cca2, by = "term", parallel = 4)

# Check for collinearity...
# ...analyse linear dependencies among constraints and conditions
# "variance inflation factors" is a diagnostic tool to identify useless constraints;
# a common rule is that values over 10 indicate redundant constraints

vif.cca(spe.cca) # all seems fine

# Make some ordiplot graphs
par(mfrow = c(1, 2))

# Scaling 1 (sites): distance triplot with wa scores
plot(spe.cca2, scaling = "sites", choices = 1:2, type = "none",
     main = "CCA (sites scaling; wa scores)")
# plot the sites
text(spe.cca2, dis = "cn", col = "red")
points(spe.cca2, display = "sites", scaling = "sites", pch = 21, col = "pink", bg = "pink", cex = 1.2)
text(spe.cca2, display = "sites", scaling = "sites", cex = 0.8, col = "grey30")
# plot the species
# points(spe.cca2, display = "species", scaling = "sites", pch = 21, col = "blue", bg = NA, cex = 1.2)
text(spe.cca2, display = "species", scaling = "sites", cex = 0.8, col = "blue")
spe.sc1 <- scores(spe.cca2, choices = 1:2, scaling = "sites", display = "sp")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "grey80")

# Scaling 2 (species; default): correlation triplot with wa scores
plot(spe.cca2, scaling = "species", choices = 1:2, type = "none",
     main = "CCA (sites scaling; wa scores)")
# plot the sites
text(spe.cca2, dis = "cn", col = "red")
points(spe.cca2, display = "sites", scaling = "species", pch = 21, col = "pink", bg = "pink", cex = 1.2)
text(spe.cca2, display = "sites", scaling = "species", cex = 0.8, col = "grey30")
# plot the species
# points(spe.cca2, display = "species", scaling = "sites", pch = 21, col = "blue", bg = NA, cex = 1.2)
text(spe.cca2, display = "species", scaling = "species", cex = 0.8, col = "blue")
spe.sc1 <- scores(spe.cca2, choices = 1:2, scaling = "species", display = "sp")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "grey80")


