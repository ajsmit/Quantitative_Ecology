################################################################################
### CHAPTER 2: EXPLORATORY DATA ANALYSIS
###
### Online supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 30 January 2017
################################################################################

# Load required packages
library(vegan)
library(labdsv)
library(tidyverse)

# Load additionnal functions
# (files must be in the working directory)
source("Num_Ecol_R_book_ed1/panelutils.R")

# Import the data from CSV files
# (files must be in the working directory)
# Species (community) data frame (fish abundances)
spe <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1))
spe <- spe[-8,]
spe
# Environmental data frame
env <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
env <- env[-8,]
env
# Spatial data frame – cartesian coordinates
spa <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpa.csv", row.names = 1))
spa <- spa[-8,]

(env_std <- as.tibble(round(decostand(env, method = "standardize"), 2)))
env_dist <- round(vegdist(env_std[,2], method = "euclidian", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist))

spe_dist <- round(vegdist(spe, method = "bray", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist))

plot.dat <- tibble(x = seq(1:30),
                   h = as.matrix(alt.dist)[,1])

library(ggplot2)
ggplot(data = plot.dat, aes(x = x, y = h)) +
  geom_point()


# Look at the moss/mite data
# Species (community) data frame
spe.m <- as.tibble(read.delim("Num_Ecol_R_book_ed1/mite.txt"))
# Environmental data frame
env.m <- as.tibble(read.delim("Num_Ecol_R_book_ed1/mite_env.txt"))
# Spatial data frame – cartesian coordinates
spa.m <- as.tibble(read.csv("Num_Ecol_R_book_ed1/mite_xy.txt", row.names = 1))

# ------------------------------------------------------------------------------
# IMPORTANT NOTE
# We originally extracted these data from the "doubs" dataset in package ade4
# v. 1.7-4, and we restored the original units of the environmental variables.
# Names of data frames and variables have changed in the package ade4 since
# the publication of the book.
# Furthermore, a mistake in the environmental dataset (as compared to the
# original publication by J. Verneaux) has not yet been corrected in the data
# currently available from ade4 version 1.7-4 (doubs$env[7,1] = 36.8, not 26.8).
# One of us (F. Gillet) returned to Verneaux's thesis to retrieve the original
# coordinates of the sites, which have also been corrected in our script.
# Therefore, we assume in our R scripts that you use the CSV files provided
# here and not the original dataset extracted from ade4!
# ------------------------------------------------------------------------------



# Basic functions
# ***************

spe			# Display the whole data frame in the console
			# Not recommended for large datasets!
spe[1:5,1:10]		# Display only 5 lines and 10 columns
head(spe)			# Display only the first 6 lines
nrow(spe)			# Number of rows (sites)
ncol(spe)			# Number of columns (species)
dim(spe)			# Dimensions of the data frame (rows, columns)
colnames(spe)		# Column labels (descriptors = species)
rownames(spe)		# Row labels (objects = sites)
summary(spe)		# Descriptive statistics for columns


# Overall distribution of abundances (dominance codes)
# ****************************************************

# Minimum and maximum of abundance values in the whole data set
range(spe)
# Count cases for each abundance class
(ab <- table(unlist(spe)))
# Create a graphic window with title
dev.new(title = "Distribution of abundance classes")
# Barplot of the distribution, all species confounded
barplot(ab, las=1, xlab="Abundance class", ylab="Frequency", col=gray(5:0/5))
# Number of absences
sum(spe==0)
# Proportion of zeros in the community data set
sum(spe==0)/(nrow(spe)*ncol(spe))


# Map of the locations of the sites
# *********************************

# New graphic window
dev.new(title="Site Locations")
# Create an empty frame (proportional axes 1:1, with titles)
# Geographic coordinates x and y from the spa data frame
plot(spa, asp=1, type="n", main="Site Locations", xlab="x 	coordinate (km)", ylab="y coordinate (km)")
# Add a blue line connecting the sites (Doubs river)
lines(spa, col="light blue")
# Add site labels
text(spa, row.names(spa), cex=0.8, col="red")
# Add text blocks
text(68, 20, "Upstream", cex=1.2, col="red")
text(15, 35, "Downstream", cex=1.2, col="red")


# Maps of some fish species
# *************************

# New graphic window (size 9x9 inches)
dev.new(title="Species Locations", width=9, height=9)
# Divide the plot window into 4 frames, 2 per row
par(mfrow=c(2,2))
# Plot four species
plot(spa, asp=1, cex.axis=0.8, col="brown", cex=spe$Satr, main="Brown trout",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, col="brown", cex=spe$Thth, main="Grayling",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, col="brown", cex=spe$Baba, main="Barbel",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, col="brown", cex=spe$Abbr, main="Common bream",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")


# Compare species: number of occurrences
# **************************************

# Compute the number of sites where each species is present
# To sum by columns, the second argument of apply(), MARGIN, is set to 2
spe.pres <- apply(spe > 0, 2, sum)
# Sort the results in increasing order
sort(spe.pres)
# Compute percentage frequencies
spe.relf <- 100*spe.pres/nrow(spe)
# Round the sorted output to 1 digit
round(sort(spe.relf), 1)
# Plot the histograms
dev.new(title="Frequency Histograms", width=8, height=5)
# Divide the window horizontally
par(mfrow=c(1,2))
hist(spe.pres, main="Species Occurrences", right=FALSE, las=1,
	xlab="Number of occurrences", ylab="Number of species",
	breaks=seq(0,30,by=5), col="bisque")
hist(spe.relf, main="Species Relative Frequencies", right=FALSE, las=1,
	xlab="Frequency of occurrences (%)", ylab="Number of species",
		breaks=seq(0, 100, by=10), col="bisque")


# Compare sites: species richness
# *******************************

# Compute the number of species at each site
# To sum by rows, the second argument of apply(), MARGIN, is set to 1
sit.pres <- apply(spe > 0, 1, sum)
# Sort the results in increasing order
sort(sit.pres)
dev.new(title="Species Richness", width=10, height=5)
par(mfrow=c(1,2))
# Plot species richness vs. position of the sites along the river
plot(sit.pres,type="s", las=1, col="gray",
	main="Species Richness vs. \n Upstream-Downstream Gradient",
	xlab="Positions of sites along the river", ylab="Species richness")
text(sit.pres, row.names(spe), cex=.8, col="red")
# Use geographic coordinates to plot a bubble map
plot(spa, asp=1, main="Map of Species Richness", pch=21, col="white",
	bg="brown", cex=5*sit.pres/max(sit.pres), xlab="x coordinate (km)",
	ylab="y coordinate (km)")
lines(spa, col="light blue")


# Compute alpha diversity indices of the fish communities
# *******************************************************

# Get help on the diversity() function
?diversity

N0 <- rowSums(spe > 0)         # Species richness
H <- diversity(spe)            # Shannon entropy
N1 <- exp(H)                   # Shannon diversity (number of abundant species)
N2 <- diversity(spe, "inv")    # Simpson diversity (number of dominant species)
J <- H/log(N0)                 # Pielou evenness
E10 <- N1/N0                   # Shannon evenness (Hill's ratio)
E20 <- N2/N0                   # Simpson evenness (Hill's ratio)
(div <- data.frame(N0, H, N1, N2, E10, E20, J))


# Transformation and standardization of the species data
# ******************************************************

# Get help on the decostand() function
?decostand

## Simple transformations

# Partial view of the raw data (abundance codes)
spe[1:5, 2:4]
# Transform abundances to presence-absence (1-0)
spe.pa <- decostand(spe, method="pa")
spe.pa[1:5, 2:4]


## Species profiles: standardization by columns

# Scale abundances by dividing them by the maximum value for each species
# Note: MARGIN=2 (column, default value) for argument "max"
spe.scal <- decostand(spe, "max")
spe.scal[1:5,2:4]
# Display the maximum by column
apply(spe.scal, 2, max)

# Scale abundances by dividing them by the species totals
# (relative abundance by species)
# Note: here, override the default MARGIN=1 argument of "total"
spe.relsp <- decostand(spe, "total", MARGIN=2)
spe.relsp[1:5,2:4]
# Display the sum by column
# Classical: apply(spe.relsp, 2, sum)
colSums(spe.relsp)

## Site profiles: standardization by rows

# Scale abundances by dividing them by the site totals
# (relative abundance by site)
spe.rel <- decostand(spe, "total") # default MARGIN=1
spe.rel[1:5,2:4]
# Display the sum of row vectors to determine if the scaling worked properly
# Classical: apply(spe.rel, 1, sum)
rowSums(spe.rel)

# Give a length of 1 to each row vector (Euclidean norm)
# This is called the chord transformation
spe.norm <- decostand(spe, "normalize") # default MARGIN=1
spe.norm[1:5,2:4]
# Verify the norm of row vectors
norm <- function(x) sqrt(x%*%x)
apply(spe.norm, 1, norm)

# Compute square root of relative abundances by site
# This is called the Hellinger transformation
spe.hel <- decostand(spe, "hellinger")
spe.hel[1:5,2:4]
# Check the norm of row vectors
apply(spe.hel, 1, norm)


## Double profiles: standardization by columns and rows

# Chi-square transformation
spe.chi <- decostand(spe, "chi.square")
spe.chi[1:5,2:4]
# Check what happened to site 8 where no species was found
spe.chi[7:9,]

# Wisconsin standardization
# Abundances are first ranged by species maxima and then by site totals
spe.wis <- wisconsin(spe)
spe.wis[1:5,2:4]


# Boxplots of transformed abundances of a common species (stone loach)
dev.new(title="Loach")
par(mfrow=c(2,2))
boxplot(spe$Babl, sqrt(spe$Babl), log1p(spe$Babl),
	las=1, main="Simple transformation",
	names=c("raw data", "sqrt", "log"), col="bisque")
boxplot(spe.scal$Babl, spe.relsp$Babl,
	las=1, main="Standardization by species",
	names=c("max", "total"), col="lightgreen")
boxplot(spe.hel$Babl, spe.rel$Babl, spe.norm$Babl,
	las=1, main="Standardization by sites",
	names=c("Hellinger", "total", "norm"), col="lightblue")
boxplot(spe.chi$Babl, spe.wis$Babl,
	las=1, main="Double standardization",
	names=c("Chi-square", "Wisconsin"), col="orange")

# Plot profiles along the upstream-downstream gradient
dev.new(title="Species profiles", width=9, height=9)
par(mfrow=c(2,2))
plot(env$dfs, spe$Satr, type="l", col=4, main="Raw data",
	xlab="Distance from the source [km]", ylab="Raw abundance code")
lines(env$dfs, spe$Thth, col=3)
lines(env$dfs, spe$Baba, col="orange")
lines(env$dfs, spe$Abbr, col=2)
lines(env$dfs, spe$Babl, col=1, lty="dotted")

plot(env$dfs, spe.scal$Satr, type="l", col=4, main="Species profiles (max)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$dfs, spe.scal$Thth, col=3)
lines(env$dfs, spe.scal$Baba, col="orange")
lines(env$dfs, spe.scal$Abbr, col=2)
lines(env$dfs, spe.scal$Babl, col=1, lty="dotted")

plot(env$dfs, spe.hel$Satr, type="l", col=4,
	main="Site profiles (Hellinger)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$dfs, spe.hel$Thth, col=3)
lines(env$dfs, spe.hel$Baba, col="orange")
lines(env$dfs, spe.hel$Abbr, col=2)
lines(env$dfs, spe.hel$Babl, col=1, lty="dotted")

plot(env$dfs, spe.chi$Satr, type="l", col=4,
	main="Double profiles (Chi-square)",
	xlab="Distance from the source [km]", ylab="Standardized abundance")
lines(env$dfs, spe.chi$Thth, col=3)
lines(env$dfs, spe.chi$Baba, col="orange")
lines(env$dfs, spe.chi$Abbr, col=2)
lines(env$dfs, spe.chi$Babl, col=1, lty="dotted")
legend("topright", c("Brown trout", "Grayling", "Barbel", "Common bream",
	"Stone loach"), col=c(4,3,"orange",2,1), lty=c(rep(1,4),3))


# Bubble maps of some environmental variables
# *******************************************

dev.new(title="Bubble maps", width=9, height=9)
par(mfrow=c(2,2))
plot(spa, asp=1, cex.axis=0.8, main="Altitude", pch=21,
	col="white", bg="red", cex=5*env$alt/max(env$alt),
	xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, main="Flow rate", pch=21,
	col="white", bg="blue", cex=5*env$flo/max(env$flo),
	xlab="x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, main="Oxygen", pch=21,
	col="white", bg="green3",cex=5*env$oxy/max(env$oxy),
	xlab= "x", ylab="y")
lines(spa, col="light blue")
plot(spa, asp=1, cex.axis=0.8, main="Nitrate", pch=21,
	col="white", bg="brown", cex=5*env$nit/max(env$nit),
	xlab="x", ylab="y")
lines(spa, col="light blue")


# Line plots
# **********

dev.new(title="Descriptor line plots")
par(mfrow=c(2,2))
plot(env$dfs, env$alt, type="l", xlab="Distance from the source (km)",
	ylab="Altitude (m)", col="red", main="Altitude")
plot(env$dfs, env$flo, type="l", xlab="Distance from the source (km)",
	ylab="Flow rate (m3/s)", col="blue", main="Flow rate")
plot(env$dfs, env$oxy, type="l", xlab="Distance from the source (km)",
	ylab="Oxygen (mg/L)", col="green3", main="Oxygen")
plot(env$dfs, env$nit, type="l", xlab="Distance from the source (km)",
	ylab="Nitrate (mg/L)", col="brown", main="Nitrate")


# Scatter plots for all pairs of environmental variables
# ******************************************************

# Bivariate plots with histograms on the diagonal and smooth fitted curves
dev.new(title="Bivariate descriptor plots")
op <- par(mfrow=c(1,1), pty="s")
pairs(env, panel=panel.smooth, diag.panel=panel.hist,
	main="Bivariate Plots with Histograms and Smooth Curves")
par(op)


# Simple transformation of an environmental variable
# **************************************************

range(env$slo)
# Log-transformation of the slope variable (y = ln(x))
# Compare histograms and boxplots of raw and transformed values
dev.new(title="Transformation and standardization of variable slope")
par(mfrow=c(2,2))
hist(env$slo, col="bisque", right=FALSE)
hist(log(env$slo), col="light green", right=F, main="Histogram of ln(env$slo)")
boxplot(env$slo, col="bisque", main="Boxplot of env$slo", ylab="env$slo")
boxplot(log(env$slo), col="light green", main="Boxplot of ln(env$slo)",
	ylab="log(env$slo)")


# Standardization of all environmental variables
# **********************************************

# Center and scale = standardize variables (z-scores)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)	# means = 0
apply(env.z, 2, sd)		# standard deviations = 1

# Same standardization using the scale() function (which returns a matrix)
env.z <- as.data.frame(scale(env))
