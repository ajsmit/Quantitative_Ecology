################################################################################
### CHAPTER 3: EUCLIDIAN DISTANCE
###
### Additional supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Prepared by : AJ Smit
### Prepared on : 31 July 2018
################################################################################


# Set-up ------------------------------------------------------------------

library(vegan)
library(ggplot2)
library(tibble)
library(gclus)
source("Num_Ecol_R_book_ed1/coldiss.R")


# Make a table with Euclidian distance ------------------------------------

ex.xy <- data.frame(row.names = letters[1:7],
                x = c(4, 5, 6, 1, 2, 8, 9),
                y = c(1, 5, 7, 4, 3, 3, 1))


# Make a plot -------------------------------------------------------------

# in base graphics
plot(x = ex.xy$x, y = ex.xy$y, type = "p", xlab = "X", ylab = "Y")

# or in ggplot
ggplot(data = ex.xy, aes(x = x, y = y, label = rownames(ex.xy))) +
  geom_point() +
  geom_label(nudge_x = 0.15, nudge_y = 0.15) +
  scale_x_continuous(breaks = c(1:9)) + scale_y_continuous(breaks = c(1:9)) +
  labs(x = "X", y = "Y")


# Calculate Euclidian distances -------------------------------------------

# example with contrived position coordinates (two dimensions)
ex.xy.euc <- vegdist(ex.xy, method = "euclidian")


# An example of higher-dimensional environmental data ---------------------

ex.env <- data.frame(row.names = letters[1:7],
                     pH = c(7.1, 7.5, 7.6, 7.0, 7.1, 7.2, 6.9),
                     O2 = c(6.5, 5.5, 5.7, 5.4, 6.3, 6.3, 6.1),
                     temp = c(12.1, 12.3, 11.9, 11.8, 12.0, 12.1, 12.2),
                     depth = c(1.1, 1.3, 1.5, 1.6, 1.8, 1.9, 2.2))

# the data are measured along different scales, so standardise
ex.env.std <- decostand(ex.env, method = "standardize")

# now we can calculate Euclidian distances
ex.env.euc <- vegdist(ex.env.std, method = "euclidian")


# Species data ------------------------------------------------------------

spe <- as_tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names=1))
spe <- spe[-8,] # once we work with the env and xy data, we need to remove site 8 too
dim(spe)
spe.bc <- round(vegdist(spe, method = "bray", diag = TRUE, upper = TRUE), 2)
spe.bc

spp.pa <- decostand(spe, method = "pa")

spe.pa <- round(vegdist(spp.pa, method = "jaccard", diag = TRUE, upper = TRUE), 2)


# but because differences of the same size are not influenced by the magnitude
# of the abundance, it makes sense to log-transform the data first
spe.log.bc <- round(vegdist(log1p(spe)), 2)

# Same but on log-transformed data
dev.new(title = "Percentage difference (Bray-Curtis), ln(y+1) data",
        width = 10, height = 5)
coldiss(spe.log.bc, byrank = FALSE, diag = TRUE)


# Correlation matrices ----------------------------------------------------

# First we need to transpose the matrix
spe[1:10, 1:10]
spe.t <- t(spe)
spe.t[1:10, 1:10]

# then Chi-square pre-transformation followed by Euclidean distance
spe.t.chi <- decostand(spe.t, "chi.square")
round(spe.t.chi[1:10, 1:10], 2)
spe.t.D16 <- dist(spe.t.chi)
round(as.matrix(spe.t.D16)[1:10, 1:10], 2)
dev.new(title="D16 on fish species (R-mode)", width=10, height=5)
coldiss(spe.t.D16, diag = TRUE)

# Jaccard index on fish presence-absence
spe.t.S7 <- vegdist(spe.t, "jaccard", binary = TRUE)
round(as.matrix(spe.t.S7)[1:10, 1:10], 2)
dev.new(title="S7 on fish species (R-mode)", width=10, height=5)
coldiss(spe.t.S7, diag=TRUE)

# Pearson r linear correlation among environmental variables
env <- as_tibble(read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
env <- env[-8,]
env.pearson <- cor(env)	# default method = "pearson"
round(env.pearson, 2)

# Reorder the variables prior to plotting
env.o <- order.single(env.pearson)

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
pairs(USJudgeRatings[1:5], panel = panel.smooth,
      cex = 1.5, pch = 24, bg = "light blue",
      diag.panel = panel.hist, cex.labels = 2, font.labels = 2)

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

dev.new(title = "Linear correlation matrix", width = 9, height = 9)
op <- par(mfrow = c(1, 1), pty = "s")
pairs(env[,env.o], lower.panel = panel.smooth, upper.panel = panel.cor,
      diag.panel = panel.hist, main = "Pearson Correlation Matrix")
par(op)
