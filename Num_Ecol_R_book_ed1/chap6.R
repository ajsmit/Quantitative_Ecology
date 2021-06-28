################################################################################
### CHAPTER 6 - CANONICAL ORDINATION
###
### Online supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 3 February 2017
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(tibble)

# Load additionnal functions
# (files must be in the working directory)
source("Num_Ecol_R_book_ed1/evplot.R")
source("Num_Ecol_R_book_ed1/hcoplot.R")
source("Num_Ecol_R_book_ed1/triplot.rda.R")
source("Num_Ecol_R_book_ed1/plot.lda.R")

# Import the data from CSV files
# (files must be in the working directory)
spe <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1))
env <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
spa <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpa.csv", row.names = 1))
# Remove empty site 8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]

# Set aside the variable 'dfs' (distance from the source) for later use
dfs <- env[, 1]

# Remove the 'dfs' variable from the env dataset
env <- env[, -1]

# Recode the slope variable (slo) into a factor (qualitative)
# variable (to show how these are handled in the ordinations)
slo2 <- rep("very_steep", nrow(env))
slo2[env$slo <= quantile(env$slo)[4]] <- "steep"
slo2[env$slo <= quantile(env$slo)[3]] <- "moderate"
slo2[env$slo <= quantile(env$slo)[2]] <- "low"
slo2 <- factor(slo2, levels = c("low", "moderate", "steep", "very_steep"))
table(slo2)

# Create an env2 data frame with slope as a qualitative variable
env2 <- env
env2$slo <- slo2

# Create two subsets of explanatory variables
# Physiography (upstream-downstream gradient)
envtopo <- env[, c(1:2)]
names(envtopo)
# Water quality
envchem <- env[, c(3:9)]
names(envchem)

# Hellinger-transform the species dataset
# Suited to species abundance data. It gives low weights to species
# with low counts and many zeros. The transformation comprises dividing
# each value in a data matrix by its row sum, and taking the square root
# of the quotient.
spe.hel <- as.tibble(decostand(spe, "hellinger"))
spe.hel


# Redundancy analysis (RDA)
# *************************

# RDA of the Hellinger-transformed fish species data, constrained
# by all the environmental variables contained in env2
(spe.rda <- rda(spe.hel ~ ., env2)) # Observe the shortcut
                                    # formula
summary(spe.rda)	# Scaling 2 (default)

# Sum of constrained eigenvalues is the constrained inertia
sum(spe.rda$CCA$eig)

# Sum of unconstrained eigenvalues is the unconstrained inertia
sum(spe.rda$CA$eig)

# Proportion of variation explained by RDA1 is therefore
spe.rda$CCA$eig[1] / (sum(spe.rda$CCA$eig) + sum(spe.rda$CA$eig))

# Canonical coefficients from the rda object
coef(spe.rda)

# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(spe.rda)$r.squared)

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

## IMPORTANT NOTE: a bug in the scaling 1 biplot scores of explanatory variables
## has been corrected in vegan 2.4-2. The arrows of the explanatory variables
## in scaling 1 biplots in the book are incorrect.

## Triplots of the rda results (wa scores)
## Scaling for species and site scores. Either species (2) or site (1) scores are scaled
## Site scores as weighted averages (vegan's default)
# Scaling 1: distance triplot
# dev.new(title = "RDA scaling 1 + wa")
#
# RDA 1 vs 2
par(mfrow = c(2, 2))
plot(spe.rda, scaling = 1,
     main = "Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "red")

# Scaling 2 (default): correlation triplot
# dev.new(title = "RDA scaling 2 + wa")
plot(spe.rda, main = "Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices = 1:2, display = "sp")
arrows(0, 0, spe.sc2[, 1] * 0.92, spe.sc2[, 2] * 0.92, length = 0, lty = 1, col = "red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(spe.rda, scaling = 1, display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(spe.rda, display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1] * 0.92, spe.sc2[,2] * 0.92, length = 0, lty = 1, col = "red")

# RDA 2 vs 3
par(mfrow = c(2, 2))
plot(spe.rda, scaling = 1, choices = 2:3,
     main = "Triplot RDA spe.hel ~ env2 - scaling 1 - wa scores")
spe.sc1 <- scores(spe.rda, choices = 2:3, scaling = 1, display = "sp")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "red")

# Scaling 2 (default): correlation triplot
# dev.new(title = "RDA scaling 2 + wa")
plot(spe.rda, choices = 2:3, main = "Triplot RDA spe.hel ~ env2 - scaling 2 - wa scores")
spe.sc2 <- scores(spe.rda, choices = 2:3, display = "sp")
arrows(0, 0, spe.sc2[, 1] * 0.92, spe.sc2[, 2] * 0.92, length = 0, lty = 1, col = "red")

## Triplots of the rda results (lc scores)
## Site scores as linear combinations of the environmental variables
# Scaling 1
# dev.new(title="RDA scaling 1 + lc")
plot(spe.rda, scaling = 1, choices = 2:3, display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc1[, 1] * 0.92, spe.sc1[, 2] * 0.92, length = 0, lty = 1, col = "red")

# Scaling 2
# dev.new(title="RDA scaling 2 + lc")
plot(spe.rda, choices = 2:3, display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env2 - scaling 2 - lc scores")
arrows(0, 0, spe.sc2[,1] * 0.92, spe.sc2[,2] * 0.92, length = 0, lty = 1, col = "red")


# Global test of the RDA result
anova(spe.rda, permutations = how(nperm = 999))

# Tests of all canonical axes
anova(spe.rda, by = "axis", permutations = how(nperm = 999))
anova(spe.rda, by = "term", permutations = how(nperm = 999))

# Variance inflation factors (VIF)
vif.cca(spe.rda)

# Apply Kaiser-Guttman criterion to residual axes
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]


# Partial RDA: effect of water chemistry, holding physiography
# constant
# ************************************************************

# Simple interface; X and W may be separate tables of quantitative variables
(spechem.physio <- rda(spe.hel, envchem, envtopo))
summary(spechem.physio)

# Formula interface; X and W variables must be in the same data frame
(spechem.physio2 <- rda(spe.hel ~ pH + har + pho + nit + amm + oxy + bod
	+ Condition(alt + slo + flo), data=env))

# Test of the partial RDA (using the results with the formula
# interface to allow the tests of the axes to be run)
anova(spechem.physio2, permutations=how(nperm=999))
anova(spechem.physio2, permutations=how(nperm=999), by="axis")


# Partial RDA triplots (with fitted site scores) - function plot.cca
# Scaling 1
dev.new(title="Partial RDA scaling 1")
plot(spechem.physio, scaling=1, display=c("sp", "lc", "cn"),
	main="Triplot RDA spe.hel ~ chem | Topo - scaling 1 - lc scores")
spe3.sc <- scores(spechem.physio, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe3.sc[, 1]*0.92, spe3.sc[, 2]*0.92, length=0, lty=1, col="red")

# Scaling 2
dev.new(title="Partial RDA scaling 2")
plot(spechem.physio, display=c("sp", "lc", "cn"),
	main="Triplot RDA spe.hel ~ chem | Topo - scaling 2 - lc scores")
spe4.sc <- scores(spechem.physio, choices=1:2, display="sp")
arrows(0, 0, spe4.sc[,1]*0.88, spe4.sc[,2]*0.88, length=0, lty=1, col="red")


# Variance inflation factors (VIF)
vif.cca(spe.rda)           # All explanatory variables
vif.cca(spechem.physio)    # Partial RDA


# Forward selection of explanatory variables
# ******************************************

# RDA with all explanatory variables
spe.rda.all <- rda(spe.hel ~ ., data=env)
# Global adjusted R^2
(R2a.all <- RsquareAdj(spe.rda.all)$adj.r.squared)

# Forward selection using forward.sel()
forward.sel(spe.hel, env, adjR2thresh=R2a.all)

# Forward selection using vegan's ordistep()
# This function allows the use of factors.
mod0 <- rda(spe.hel ~ 1, data=env)
step.forward <- ordistep(mod0,scope=formula(spe.rda.all),
	direction="forward", permutations=how(nperm=199))
RsquareAdj(step.forward)

# Compare adjusted R^2
RsquareAdj(rda(spe.hel ~ alt, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt + oxy, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt + oxy + bod, data=env))$adj.r.squared
RsquareAdj(rda(spe.hel ~ alt + oxy + bod + slo, data=env))$adj.r.squared

# Forward selection using vegan's ordiR2step()
# using a double stopping criterion (Blanchet et al. 2008a)
step.forward2 <- ordiR2step(rda(spe.hel ~ 1, data=env),
	scope=formula(spe.rda.all ), direction="forward",
	permutations=how(nperm=199))


# Parsimonious RDA
# ****************

(spe.rda.pars <- rda(spe.hel ~ alt + oxy + bod, data=env))
anova(spe.rda.pars, permutations=how(nperm=999))
anova(spe.rda.pars, permutations=how(nperm=999), by="axis")
(R2a.pars <- RsquareAdj(spe.rda.pars)$adj.r.squared)
# Compare variance inflation factors
vif.cca(spe.rda.all)
vif.cca(spe.rda.pars)

# Triplots of the parsimonious RDA (with fitted site scores)
dev.new(title="Parsimonious RDA scaling 1", width=14, height=7)
par(mfrow=c(1,2))
# Scaling 1
plot(spe.rda.pars, scaling=1, display=c("sp", "lc", "cn"),
	main="Triplot RDA spe.hel ~ alt+oxy+bod - scaling 1 - lc scores")
spe4.sc <- scores(spe.rda.pars, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe4.sc[, 1]*0.92, spe4.sc[, 2]*0.92, length=0, lty=1, col="red")
# Scaling 2
plot(spe.rda.pars, display=c("sp", "lc", "cn"),
	main="Triplot RDA spe.hel ~ alt+oxy+bod - scaling 2 - lc scores")
spe5.sc <- scores(spe.rda.pars, choices=1:2, display="sp")
arrows(0, 0, spe5.sc[,1]*0.9, spe5.sc[,2]*0.9, length=0, lty=1, col="red")


# Variation partitioning with two sets of explanatory variables
# *************************************************************

# Explanation of fraction labels (two, three and four explanatory matrices)
dev.new(title="Symbols of variation partitioning fractions", width=12, height=4)
par(mfrow=c(1,3))
showvarparts(2)
showvarparts(3)
showvarparts(4)

## 1. Variation partitioning with all explanatory variables
(spe.part.all <- varpart(spe.hel, envchem, envtopo))

# Plot of the partitioning results
dev.new(title="Variation partitioning - all variables")
plot(spe.part.all, digits=2)


## 2. Variation partitioning after forward selection of explanatory variables
# Separate forward selection in each subset of environmental variables
spe.chem <- rda(spe.hel, envchem)
R2a.all.chem <- RsquareAdj(spe.chem)$adj.r.squared
forward.sel(spe.hel, envchem, adjR2thresh=R2a.all.chem, nperm=9999)

spe.topo <- rda(spe.hel, envtopo)
R2a.all.topo <- RsquareAdj(spe.topo)$adj.r.squared
forward.sel(spe.hel, envtopo, adjR2thresh=R2a.all.topo, nperm=9999)

# Parsimonious subsets of explanatory variables (based on forward selections)
names(envchem)
envchem.pars <- envchem[, c(4,6,7)]
names(envtopo)
envtopo.pars <- envtopo[, c(1,2)]

# Variation partitioning
(spe.part <- varpart(spe.hel, envchem.pars, envtopo.pars))
dev.new(title="Variation partitioning - parsimonious subsets")
plot(spe.part, digits=2)

# Tests of all testable fractions
# Test of fractions [a+b]
anova(rda(spe.hel, envchem.pars), permutations=how(nperm=999))
# Test of fractions [b+c]
anova(rda(spe.hel, envtopo.pars), permutations=how(nperm=999))
# Test of fractions [a+b+c]
env.pars <- cbind(envchem.pars, envtopo.pars)
anova(rda(spe.hel, env.pars), permutations=how(nperm=999))
# Test of fraction [a]
anova(rda(spe.hel, envchem.pars, envtopo.pars), permutations=how(nperm=999))
# Test of fraction [c]
anova(rda(spe.hel, envtopo.pars, envchem.pars), permutations=how(nperm=999))


## 3. Variation partitioning without the 'nit' variable
envchem.pars2 <- envchem[, c(6,7)]
(spe.part2 <- varpart(spe.hel, envchem.pars2, envtopo.pars))
dev.new(title="Variation partitioning - parsimonious subset 2")
plot(spe.part2, digits=2)


# Two-way MANOVA by RDA
# *********************

# Creation of a factor 'altitude' (3 levels, 9 sites each)
alt.fac <- gl(3, 9,labels=c("high", "mid", "low"))
# Creation of a factor mimicking 'pH'
pH.fac <- as.factor(c(1, 2, 3, 2, 3, 1, 3, 2, 1, 2, 1, 3, 3, 2,
    1, 1, 2, 3, 2, 1, 2, 3, 2, 1, 1, 3, 3))
# Are the factors balanced?
table(alt.fac, pH.fac)

# Creation of Helmert contrasts for the factors and their interaction
alt.pH.helm <- model.matrix(~ alt.fac * pH.fac,
	contrasts=list(alt.fac="contr.helmert", pH.fac="contr.helmert"))[,-1]
alt.pH.helm

# Check property 1 of Helmert contrasts: all variables sum to 0
apply(alt.pH.helm, 2, sum)
# Check property 2 of Helmert contrasts: variables are uncorrelated
cor(alt.pH.helm)

# Verify multivariate homogeneity of within-group covariance
# matrices using the betadisper() function (vegan package)
# implementing Marti Anderson's testing method

spe.hel.d1 <- dist(spe.hel[1:27,])

# Factor "altitude"
(spe.hel.alt.MHV <- betadisper(spe.hel.d1, alt.fac))
anova(spe.hel.alt.MHV)     # Parametric test (not recommended here)
permutest(spe.hel.alt.MHV) # Permutation test
# Factor "pH"
(spe.hel.pH.MHV <- betadisper(spe.hel.d1, pH.fac))
anova(spe.hel.pH.MHV)     # Parametric test (not recommended here)
permutest(spe.hel.pH.MHV) # Permutation test

# Test the interaction first. The factors alt and pH form
# the matrix of covariables
interaction.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 5:8], alt.pH.helm[, 1:4])
anova(interaction.rda, permutations=how(nperm=999))

# Test the main factor alt. The factor pH and the interaction
# form the matrix of covariables.
factor.alt.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 1:2], alt.pH.helm[, 3:8])
anova(factor.alt.rda, permutations=how(nperm=999), strata=pH.fac)

# Test the main factor pH. The factor alt and the interaction
# form the matrix of covariables.
factor.pH.rda <- rda(spe.hel[1:27, ], alt.pH.helm[, 3:4],
	alt.pH.helm[, c(1:2, 5:8)])
anova(factor.pH.rda, permutations=how(nperm=999), strata=alt.fac)

# RDA and triplot for the significant factor alt
alt.rda.out <- rda(spe.hel[1:27,]~., as.data.frame(alt.fac))
dev.new(title="Multivariate ANOVA - altitude")
plot(alt.rda.out, scaling=1, display=c("sp", "wa", "cn"),
	main="Multivariate ANOVA, factor altitude - scaling 1 - wa scores")
spe.manova.sc <- scores(alt.rda.out, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.manova.sc[, 1], spe.manova.sc[, 2], length=0, col="red")


# Distance-based redundancy analysis (db-RDA)
# ******************************************

# 1. Explicit steps
spe.bray <- vegdist(spe[1:27, ], "bray")
spe.pcoa <- cmdscale(spe.bray, k=nrow(spe[1:27, ])-1, eig=TRUE, add=TRUE)
spe.scores <- spe.pcoa$points
# Test of the interaction. Factors alt and pH, Helmert-coded, form the
# matrix of covariables
interact.dbrda <- rda(spe.scores, alt.pH.helm[, 5:8],
	alt.pH.helm[, 1:4])
anova(interact.dbrda, permutations=how(nperm=999))

# 2. Direct way using vegan's function capscale. Runs with model
#    interface only.
#    Response matrix can be raw data:
interact.dbrda2 <- capscale(spe[1:27,] ~ alt.fac*pH.fac +
	Condition(alt.fac+pH.fac), distance="bray", add=TRUE)
anova(interact.dbrda2, permutations=how(nperm=999))
# Computation using Helmert contrasts:
# interact.dbrda2 <- capscale(spe[1:27, ] ~ alt.fac1:pH.fac1 + alt.fac2:pH.fac1
#  + alt.fac1:pH.fac2 + alt.fac2:pH.fac2 + Condition(alt.fac1 + alt.fac2
#  + pH.fac1 + pH.fac2), data=as.data.frame(alt.pH.helm),
#  distance="bray", add=TRUE)

#    Or response matrix can be a dissimilarity matrix:
interact.dbrda3 <- capscale(spe.bray ~ alt.fac*pH.fac
	+ Condition(alt.fac+pH.fac), add=TRUE)
anova(interact.dbrda3, permutations=how(nperm=999))
# Computation using Helmert contrasts:
# interact.dbrda3 <- capscale(spe.bray ~ alt.fac1:pH.fac1 + alt.fac2:pH.fac1
# + alt.fac1:pH.fac2 + alt.fac2:pH.fac2
# + Condition(alt.fac1 + alt.fac2 + pH.fac1 + pH.fac2),
# data=as.data.frame(alt.pH.helm), add=TRUE)


# RDA with a second degree explanatory variable
# *********************************************

# Create a matrix of dfs and its orthogonal second degree term using poly()
dfs.df <- poly(dfs, 2)
colnames(dfs.df) <- c("dfs", "dfs2")
# Verify if both variables are significant
forward.sel(spe, dfs.df)

# RDA and test
spe.dfs.rda <- rda(spe ~ ., as.data.frame(dfs.df))
anova(spe.dfs.rda, permutations=how(nperm=999))

# Triplot using lc (model) site scores and scaling 2
dev.new(title="RDA w. 2nd-order variables - scaling 2")
plot(spe.dfs.rda, scaling=2, display=c("sp", "lc", "cn"),
	main="Triplot RDA spe ~ dfs+dfs2 - scaling 2 - lc scores")
spe6.sc <- scores(spe.dfs.rda, choices=1:2, scaling=2, display="sp")
arrows(0, 0, spe6.sc[, 1], spe6.sc[, 2], length=0, lty=1, col="red")

# Maps of four fish species
dev.new(title="Four fish species", 9, 9)
par(mfrow=c(2,2))
plot(spa, asp=1, col="brown", cex=spe$Satr,
	xlab="x (km)", ylab="y (km)", main="Brown trout")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$Thth,
	xlab="x (km)", ylab="y (km)", main="Grayling")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$Alal,
	xlab="x (km)", ylab="y (km)", main="Bleak")
lines(spa, col="light blue")
plot(spa, asp=1, col="brown", cex=spe$Titi,
	xlab="x (km)", ylab="y (km)", main="Tench")
lines(spa, col="light blue")




### Code it yourself corner 4 ##################################################

myRDA <- function(Y,X)
{

	# 1. Preparation of the data
	# **************************

	Y.mat <- as.matrix(Y)
	Yc <- scale(Y.mat, scale=FALSE)

	X.mat <- as.matrix(X)
	Xcr <- scale(X.mat)

	# Dimensions
	n <- nrow(Y)
	p <- ncol(Y)
	m <- ncol(X)

	# 2. Computation of the multivariate linear regression
	# ****************************************************

	# Matrix of regression coefficients (eq. 11.9)
	B <- solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% Yc

	# Matrix of fitted values (eq. 11.10)
	Yhat <- Xcr %*% B

	# Matrix of residuals
	Yres <- Yc - Yhat


	# 3. PCA on fitted values
	# ***********************
	# Covariance matrix (eq. 11.12)
	S <- cov(Yhat)

	# Eigenvalue decomposition
	eigenS <- eigen(S)

	# How many canonical axes?
	kc <- length(which(eigenS$values > 0.00000001))

	# Eigenvalues of canonical axes
	ev <- eigenS$values[1:kc]
	# Total variance (inertia) of the centred matrix Yc
	trace = sum(diag(cov(Yc)))

	# Orthonormal eigenvectors (contributions of response variables,
	# scaling 1)
	U <- eigenS$vectors[,1:kc]
	row.names(U) <- colnames(Y)

	# Site scores (vegan's 'wa' scores, scaling 1; eq. 11.17)
	F <- Yc %*% U
	row.names(F) <- row.names(Y)

	# Site constraints (vegan's 'lc' scores, scaling 1; eq. 11.18)
	Z <- Yhat %*% U
	row.names(Z) <- row.names(Y)

	# Canonical coefficients (eq. 11.19)
	CC <- B %*% U
	row.names(CC) <- colnames(X)

	# Explanatory variables
	# Species-environment correlations
	corXZ <- cor(X,Z)

	# Diagonal matrix of weights
	D <- diag(sqrt(ev/trace))

	# Biplot scores of explanatory variables
	coordX <- corXZ %*% D    # Scaling 1
	coordX2 <- corXZ         # Scaling 2
	row.names(coordX) <- colnames(X)
	row.names(coordX2) <- colnames(X)

	# Scaling to sqrt of the relative eigenvalue (for scaling 2)
	U2 <- U %*% diag(sqrt(ev))
	row.names(U2) <- colnames(Y)
	F2 <- F %*% diag(1/sqrt(ev))
	row.names(F2) <- row.names(Y)
	Z2 <- Z %*% diag(1/sqrt(ev))
	row.names(Z2) <- row.names(Y)

	# Unadjusted R2
	R2 <- sum(ev/trace)
	# Adjusted R2
	R2a <- 1-((n-1)/(n-m-1))*(1-R2)


	# 4. PCA on residuals
	# *******************
	# Write your own code as in Chapter 5. It could begin with:
	#     eigenSres <- eigen(cov(Yres))
	#     evr <- eigenSres$values


	# 5. Output

	result <- list(trace, R2, R2a, ev, CC, U, F, Z, coordX, U2, F2, Z2, coordX2)
	names(result) <- c("Total_variance", "R2", "R2adj", "Can_ev", "Can_coeff",
    "Species_sc1", "wa_sc1", "lc_sc1", "Biplot_sc1", "Species_sc2",
    "wa_sc2", "lc_sc2", "Biplot_sc2")

	result
}


doubs.myRDA <- myRDA(spe.hel, env)
summary(doubs.myRDA)

################################################################################



# Canonical correspondence analysis (CCA)
# ***************************************

# CCA of the raw fish species data, constrained by all the
# environmental variables in env2
(spe.cca <- cca(spe ~ ., env2))
summary(spe.cca)	# Scaling 2 (default)

# Unadjusted R^2
RsquareAdj(spe.cca)

## CCA triplots (using lc site scores)
dev.new(title="CCA triplot - lc scores", width=14, height=7)
par(mfrow=c(1,2))
# Scaling 1: species scores scaled to relative eigenvalues,
# sites are weighted averages of the species
plot(spe.cca, scaling=1, display=c("sp","lc","cn"),
	main="Triplot CCA spe ~ env2 - scaling 1")
# Default scaling 2: site scores scaled to relative eigenvalues,
# species are weighted averages of the sites
plot(spe.cca, display=c("sp","lc","cn"),
	main="Triplot CCA spe ~ env2 - scaling 2")

dev.new(title="CCA biplot - without species", width=14, height=7)
par(mfrow=c(1,2))
# CCA scaling 1 biplot without species (using lc site scores)
plot(spe.cca, scaling=1, display=c("lc", "cn"),
	main="Biplot CCA spe ~ env2 - scaling 1")
# CCA scaling 2 biplot with species but without sites
plot(spe.cca, scaling=2, display=c("sp", "cn"),
	main="Biplot CCA spe ~ env2 - scaling 2")

# Permutation test of the overall analysis
anova(spe.cca, permutations=how(nperm=999))
# Permutation test of each axis
anova(spe.cca, by="axis", permutations=how(nperm=999))


# CCA-based forward selection using vegan's ordistep()
# ****************************************************

# This function allows the use of factors like 'slo' in env2
cca.step.forward <- ordistep(cca(spe ~ 1, data=env2), scope=formula(spe.cca),
	direction="forward", permutations=how(nperm=199))
cca.step.forward


# Parsimonious CCA using alt, oxy and bod
# ***************************************

spe.cca.pars <- cca.step.forward
# Or rerun the analysis explicitly:
# spe.cca.pars <- cca(spe ~ alt + oxy + bod, data=env2)
anova(spe.cca.pars, permutations=how(nperm=999))
anova(spe.cca.pars, permutations=how(nperm=999), by="axis")
RsquareAdj(spe.cca.pars)

# Compare variance inflation factors
vif.cca(spe.cca)
vif.cca(spe.cca.pars)


# Three-dimensional interactive ordination plots
# (requires vegan3d package)
# **********************************************

# Plot of the sites only (wa scores)
ordirgl(spe.cca.pars, type="t", scaling=1)

# Connect weighted average scores to linear combination scores
orglspider(spe.cca.pars, scaling=1, col="purple")

# Plot the sites (wa scores) with a clustering result
# Colour sites according to cluster membership
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward.D2"), 4)
ordirgl(spe.cca.pars, type="t", scaling=1, ax.col="black", col=gr+1)
# Connect sites to cluster centroids
orglspider(spe.cca.pars, gr, scaling=1)

# Complete CCA 3D triplot
ordirgl(spe.cca.pars, type="t", scaling=2)
orgltext(spe.cca.pars, display="species", type="t", scaling=2, col="cyan")

# Plot species groups (Jaccard similarity, useable in R mode)
gs <- cutree(hclust(vegdist(t(spe), method="jaccard"), "ward.D2"), 4)
ordirgl(spe.cca.pars, display="species", type="t", col=gs+1)


# Linear discriminant analysis (LDA)
# **********************************

# Ward clustering result of Hellinger-transformed species data,
# cut in 4 groups
gr <- cutree(hclust(vegdist(spe.hel, "euc"), "ward.D2"), 4)

# Environmental matrix with only 3 variables (alt, oxy and bod)
env.pars2 <- as.matrix(env[, c(1, 9, 10)])

# Verify multivariate homogeneity of within-group covariance matrices
# using the betadisper() function {vegan}
env.pars2.d1 <- dist(env.pars2)
(env.MHV <- betadisper(env.pars2.d1, gr))
anova(env.MHV)
permutest(env.MHV)	# Permutational test (recommended)

# Log transform alt and dbo
env.pars3 <- cbind(log(env$alt), env$oxy, log(env$bod))
colnames(env.pars3) <- c("alt.ln", "oxy", "bod.ln")
rownames(env.pars3) <- rownames(env)
env.pars3.d1 <- dist(env.pars3)
(env.MHV2 <- betadisper(env.pars3.d1, gr))
permutest(env.MHV2)


# METHODOLOGICAL UPDATE below to clarify difference between identification
# and discrimination functions

# Computation of LDA - identification functions (on unstandardized variables)
# ---------------------------------------------------------------------------

env.pars3.df <- as.data.frame(env.pars3)
(spe.lda <- lda(gr ~ alt.ln + oxy + bod.ln, data=env.pars3.df))
# Alternate coding without formula interface: spe.lda <- lda(env.pars3.df, gr)
# The result object contains the information necessary to interpret the LDA
summary(spe.lda)

# Display the group means for the 3 variables
spe.lda$means

# Extract the identification functions (matrix C, eq. 11.33)
(C <- spe.lda$scaling)

# Classification of two new objects (identification)
# A new object is created with two sites; (1) ln(alt)=6.8, oxygen=9 and ln(bod)=0.8
# and (2) ln(alt)=5.5, oxygen=10 and ln(bod)=1.0
newo <- data.frame(c(6.8,5.5), c(9,10), c(0.8,1))
colnames(newo) <- colnames(env.pars3)
newo
(predict.new <- predict(spe.lda, newdata=newo))

# Computation of LDA - discrimination functions (on standardized variables)
# -------------------------------------------------------------------------

env.pars3.sc <- as.data.frame(scale(env.pars3.df))
spe.lda2 <- lda(gr ~., data=env.pars3.sc)

# Display the group means for the 3 variables
spe.lda2$means

# Extract the classification functions
(C2 <- spe.lda2$scaling)

# Compute the canonical eigenvalues
spe.lda2$svd^2

# Position the objects in the space of the canonical variates
(Fp2 <- predict(spe.lda2)$x)
# alternative way: Fp2 <- as.matrix(env.pars3.sc) %*% C2

# Classification of the objects
(spe.class2 <- predict(spe.lda2)$class)

# Posterior probabilities of the objects to belong to the groups
# (rounded for easier interpretation)
(spe.post2 <- round(predict(spe.lda2)$posterior,2))

# Contingency table of prior versus predicted classifications
(spe.table2 <- table(gr, spe.class2))

# Proportion of correct classification (classification success)
diag(prop.table(spe.table2, 1))

# Plot the LDA results using homemade function plot.lda.R
dev.new(title="Discriminant analysis")
plot.lda(lda.out=spe.lda2, groups=gr, plot.sites=2, plot.centroids=1, mul.coef=2.35)

# LDA with jackknife-based classification (i.e., leave-one-out cross-validation)
(spe.lda.jac <- lda(gr ~ alt.ln + oxy + bod.ln, data=env.pars3.sc, CV=TRUE))
summary(spe.lda.jac)

# Numbers and proportions of correct classification
spe.jac.class <- spe.lda.jac$class
spe.jac.table <- table(gr, spe.jac.class)
# Classification success
diag(prop.table(spe.jac.table, 1))



# Canonical correlation analysis (CCorA)
# **************************************

# Preparation of data (transformations to make variable
# distributions approximately symmetrical)
envchem2 <- envchem
envchem2$pho <- log(envchem$pho)
envchem2$nit <- sqrt(envchem$nit)
envchem2$amm <- log1p(envchem$amm)
envchem2$bod <- log(envchem$bod)
envtopo2 <- envtopo
envtopo2$alt <- log(envtopo$alt)
envtopo2$slo <- log(envtopo$slo)
envtopo2$flo <- sqrt(envtopo$flo)

# CCorA (on standardized variables)
chem.topo.ccora <- CCorA(envchem2, envtopo2, stand.Y=TRUE, stand.X=TRUE,
	permutations=how(nperm=999))
chem.topo.ccora
dev.new(title="Canonical correlation analysis", width=12, height=8)
biplot(chem.topo.ccora, plot.type="biplot")



# Co-inertia analysis
# *******************

# PCA on both matrices using ade4 functions
dudi.chem <- dudi.pca(envchem2, scale=TRUE, scan=FALSE, nf=3)
dudi.topo <- dudi.pca(envtopo2, scale=TRUE, scan=FALSE, nf=2)
cumsum(dudi.chem$eig/sum(dudi.chem$eig)) # Cumulated relative variation of eigenvalues
cumsum(dudi.topo$eig/sum(dudi.topo$eig)) # Cumulated relative variation of eigenvalues
all.equal(dudi.chem$lw, dudi.topo$lw)    # Equal row weights in the 2 analyses?

# Co-inertia analysis
(coia.chem.topo <- coinertia(dudi.chem, dudi.topo, scan=FALSE, nf=2))
summary(coia.chem.topo)

# Relative variation on first eigenvalue
coia.chem.topo$eig[1]/sum(coia.chem.topo$eig)
# Permutation test
randtest(coia.chem.topo, nrepet=999)

# Plot results
dev.new(title="Co-inertia analysis")
plot(coia.chem.topo)


# Multiple factor analysis (MFA)
# ******************************

# MFA on 3 groups of variables:
# Regroup the 3 tables (Hellinger-transformed species,
# physiographical variables, chemical variables)
tab3 <- data.frame(spe.hel, envtopo, envchem)
dim(tab3)
# Number of variables in each group
(grn <- c(ncol(spe), ncol(envtopo), ncol(envchem)))

# Close the previous graphic windows
graphics.off()

# Compute the MFA with multiple plots
t3.mfa <- MFA(tab3, group=grn, type=c("c","s","s"), ncp=2,
	name.group=c("Fish community","Physiography","Water quality"))
t3.mfa

# Alternative plots
dev.new(title="My MFA plots")
par(mfrow=c(2,2))
plot(t3.mfa, choix="ind", habillage="none", cex=0.5)
plot(t3.mfa, choix="ind", habillage="none", partial="all", cex=0.5)
plot(t3.mfa, choix="var", habillage="group", cex=0.5)
plot(t3.mfa, choix="axes", cex=0.5)

# RV coefficients with tests (p-values above the diagonal of the matrix)
rvp <- t3.mfa$group$RV
rvp[1,2] <- coeffRV(spe.hel, scale(envtopo))$p.value
rvp[1,3] <- coeffRV(spe.hel, scale(envchem))$p.value
rvp[2,3] <- coeffRV(scale(envtopo), scale(envchem))$p.value
round(rvp[-4,-4], 6)

# Eigenvalues and % of variation
t3.mfa$eig
ev <- t3.mfa$eig[,1]
names(ev) <- 1:nrow(t3.mfa$eig)
dev.new(title="MFA eigenvalues")
evplot(ev)

# Alternative to the standard, automatic plots:
# Plot only the significant variables (correlations)

# Select the most characteristic variables
aa <- dimdesc(t3.mfa, axes=1:2, proba=0.0001)

# Plot
dev.new(title="MFA: correlations among significant variables")
varsig <- t3.mfa$quanti.var$cor[unique(c(rownames(aa$Dim.1$quanti),
	rownames(aa$Dim.2$quanti))),]
plot(varsig[,1:2], asp=1, type="n", xlim=c(-1,1), ylim=c(-1,1))
abline(h=0, lty=3)
abline(v=0, lty=3)
symbols(0, 0, circles=1, inches=FALSE, add=TRUE)
arrows(0, 0, varsig[,1], varsig[,2], length=0.08, angle=20)
for (v in 1:nrow(varsig))
{
	if (abs(varsig[v,1]) > abs(varsig[v,2]))
	{
		if (varsig[v,1] >= 0) pos <- 4
		else pos <- 2
	}
	else
	{
		if (varsig[v,2] >= 0) pos <- 3
		else pos <- 1
	}
	text(varsig[v,1], varsig[v,2], labels=rownames(varsig)[v], pos=pos)
}



