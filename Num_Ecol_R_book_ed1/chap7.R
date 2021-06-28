################################################################################
### CHAPTER 7: SPATIAL ANALYSIS
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 6 February 2017
################################################################################


# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ape)
library(spdep)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(AEM)


# Load additionnal functions
# (files must be in the working directory)
source("plot.links.R")
source("sr.value.R")
source("quickMEM.R")

# Import the data from CSV files
# (files must be in the working directory)
mite <- read.table("mite.txt")
mite.env <- read.table("mite_env.txt")
mite.xy <- read.table("mite_xy.txt")

# Transform the data
mite.h <- decostand (mite, "hellinger")
mite.xy.c <- scale(mite.xy, center=TRUE, scale=FALSE)



# Spatial correlogram (based on Moran's I)
# ****************************************

# Search for neighbours of all points within a radius of 0.7 m
# and multiples (i.e., 0 to 0.7 m, 0.7 to 1.4 m and so on). The points do not 
# form a connected graph at 0.7 m.
dev.new(title="Linkage map")
plot.links(mite.xy, thresh=0.7)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, 0.7)
summary(nb1)

# Correlogram of substrate density
subs.dens <- mite.env[,1]
subs.correlog <- sp.correlogram(nb1, subs.dens, order=14, method="I", 
	zero.policy=TRUE)
print(subs.correlog, p.adj.method="holm")
dev.new(title="Correlogram of substrate density")
plot(subs.correlog)


# Mantel correlogram of the oribatid mite data
# ********************************************

# The species data are first detrended; see Section 7.3
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
mite.h.D1 <- dist(mite.h.det)
(mite.correlog <- mantel.correlog(mite.h.D1, XY=mite.xy, nperm=999))
summary(mite.correlog)

# Number of classes
mite.correlog$n.class # or: mite.correlog[2]
# Break points
mite.correlog$break.pts # or: mite.correlog[3]

# Plot the Mantel correlogram
dev.new(title="Mantel correlogram of mite data", width=12, height=7)
plot(mite.correlog)


# Trend-surface analysis
# **********************

## Simple models on a square, regularly sampled surface

# Construct and plot a 10 by 10 grid
xygrid <- expand.grid(1:10, 1:10)
dev.new(title="Regular grid")
plot(xygrid)

# Centring
xygrid.c <- scale(xygrid, scale=FALSE)

# Create and plot some first, second and third-degree functions of X and Y
X <- xygrid.c[,1]
Y <- xygrid.c[,2]
XY <- X+Y
XY2 <- X^2 + Y^2
XY3 <- X^2 - X*Y - Y^2
XY4 <- X+Y + X^2 + X*Y + Y^2
XY5 <- X^3 + Y^3
XY6 <- X^3 + X^2*Y + X*Y^2 + Y^3
XY7 <- X + Y + X^2 + X*Y + Y^2 + X^3 + X^2*Y + X*Y^2 + Y^3
xy3deg <- cbind(X, Y, XY, XY2, XY3, XY4, XY5, XY6, XY7)

dev.new(title="Polynomials")
s.value(xygrid, xy3deg, symbol="circle")
	# Try other combinations, for instance with minus signs or with
	# coefficients not equal to 1.


## Trend-surface analysis of the mite data

# Computation of a raw (non-orthogonal) third-degree polynomial 
# function on the previously centred X-Y coordinates
mite.poly <- poly(as.matrix(mite.xy.c), degree=3, raw=TRUE)
colnames(mite.poly) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")

# RDA with all 9 polynomial terms
(mite.trend.rda <- rda(mite.h ~ ., data=as.data.frame(mite.poly)))

# Computation of the adjusted R^2
(R2adj.poly <- RsquareAdj(mite.trend.rda)$adj.r.squared)

# RDA using a third-degree orthogonal polynomial of the geographic coordinates
mite.poly.ortho <- poly(as.matrix(mite.xy), degree=3)
colnames(mite.poly.ortho) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")
(mite.trend.rda.ortho <- rda(mite.h ~ ., data=as.data.frame(mite.poly.ortho)))
(R2adj.poly2 <- RsquareAdj(mite.trend.rda.ortho)$adj.r.squared)

# Forward selection using Blanchet et al. (2008a) double stopping criterion
(mite.trend.fwd <- forward.sel(mite.h, mite.poly.ortho, adjR2thresh=R2adj.poly2))

# New RDA using the 6 terms retained
(mite.trend.rda2 <- rda(mite.h ~ ., 
	data=as.data.frame(mite.poly)[,mite.trend.fwd[,2]]))

# Overall test and test of the canonical axes
anova(mite.trend.rda2)
anova(mite.trend.rda2, by="axis")

# Plot of the three independent significant spatial structures
# (canonical axes) plus the fourth (p-value around 0.06). 
mite.trend.fit <- scores(mite.trend.rda2, choices=1:4, display="lc", 
	scaling=1)
dev.new(title="Mite Trend Surface Analysis")
s.value(mite.xy, mite.trend.fit, symbol="circle")

# Other plotting code with homemade function sr.value:
# dev.new(title="Mite Trend Surface Analysis")
# par(mfrow=c(1,3))
# sr.value(mite.xy,mite.trend.fit[,1])
# sr.value(mite.xy,mite.trend.fit[,2])
# sr.value(mite.xy,mite.trend.fit[,3])

# ------------------------------------------------------------------------------
	# If you want to construct a raw polynomial function directly
	# within the rda call, here is the syntax (2nd degree):
	# mite.trend.rda <- rda(mite.h ~ Xm + Ym + I(Xm^2) + I(Xm*Ym) + I(Ym^2))
	# Notice how squared variables and product variables are 
	# requested to be treated "as they are" by function I(). 
	# Otherwise R would consider them as ANOVA terms.
# ------------------------------------------------------------------------------


# METHODOLOGICAL UPDATE: first-generation PCNM replaced by dbMEM
# --------------------------------------------------------------

# dbMEM analysis (artificial data)
# ********************************

# 1. One-dimensional sampling: transect with 100 equispaced points.
#    The distance between adjacent points is 1. Function dbmem() automatically
#    computes the threshold value, 1 in this case.

# Generate transect points
tr100 <- 1:100

# Creation of the dbMEM eigenfunctions with Moran's I corresponding to positive
# spatial correlation: argument MEM.autocor = "positive", which is the default
tr100.dbmem.tmp <- dbmem(tr100, silent=FALSE)
tr100.dbmem <- as.data.frame(tr100.dbmem.tmp)
# Display the eigenvalues
attributes(tr100.dbmem.tmp)$values
# Number of (positive) eigenvalues
length(attributes(tr100.dbmem.tmp)$values)

# Plot some dbMEM variables modelling positive spatial correlation
dev.new(title="dbMEM variables (transect)")
par(mfrow=c(4,2))
somedbmem <- c(1, 2, 4, 8, 15, 20, 30, 40)
for(i in 1:length(somedbmem)){
	plot(tr100.dbmem[,somedbmem[i]], type="l", xlab="X coordinate", ylab=c("dbMEM", somedbmem[i]))
}

# 2. Two-dimensional sampling: equispaced grid with smallest distance between 
#    points equal to 1.

# Generate grid point coordinates
xygrid2 <- expand.grid(1:20, 1:20)

# Creation of the dbMEM eigenfunctions with positive Moran's I 
xygrid2.dbmem.tmp <- dbmem(xygrid2)
xygrid2.dbmem <- as.data.frame(xygrid2.dbmem.tmp)
# Count the eigenvalues
length(attributes(xygrid2.dbmem.tmp)$values)
# Plot some dbMEM variables using s.value {adegraphics}
dev.new(title="dbMEM variables (grid)", width=8, height=8)
somedbmem2 <- c(1, 2, 5, 10, 20, 40, 80, 120, 189)
s.value(xygrid2, xygrid2.dbmem[,somedbmem2], method="color", symbol="circle", ppoints.cex=0.5)


# dbMEM analysis of the oribatid mite data
# ****************************************

# Is there a linear trend in the mite data?
anova(rda(mite.h, mite.xy))	# Result: significant trend
# Computation of linearly detrended mite data
mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))

## 1. Construct the matrix of dbMEM variables
mite.dbmem.tmp <- dbmem(mite.xy, silent=FALSE)
mite.dbmem <- as.data.frame(mite.dbmem.tmp)
# Truncation distance used above:
(thr <- give.thresh(dist(mite.xy)))

# Display and count the eigenvalues
attributes(mite.dbmem.tmp)$values
length(attributes(mite.dbmem.tmp)$values)
 # Argument silent=FALSE allows the function to display the truncation level.

## 2. Run the global dbMEM analysis on the *detrended* Hellinger-transformed
##    mite data
(mite.dbmem.rda <- rda(mite.h.det, mite.dbmem))
anova(mite.dbmem.rda)

## 3. Since the analysis is significant, compute the adjusted R2
##    and run a forward selection of the dbmem variables
(mite.R2a <- RsquareAdj(mite.dbmem.rda)$adj.r.squared)
(mite.dbmem.fwd <- forward.sel(mite.h.det, as.matrix(mite.dbmem), 
	adjR2thresh=mite.R2a))
(nb.sig.dbmem <- nrow(mite.dbmem.fwd))	# Number of signif. dbMEM
# Identity of significant dbMEM increasing order
(dbmem.sign <- sort(mite.dbmem.fwd[,2]))
# Write the significant dbMEM to a new object
dbmem.red <- mite.dbmem[,c(dbmem.sign)]

## 4. New dbMEM analysis with 8 significant dbMEM variables
##    Adjusted R-square after forward selection: R2adj=0.2418
(mite.dbmem.rda2 <- rda(mite.h.det ~ ., data=dbmem.red))
(mite.fwd.R2a <- RsquareAdj(mite.dbmem.rda2)$adj.r.squared)
anova(mite.dbmem.rda2)
(axes.test <- anova(mite.dbmem.rda2, by="axis"))
# Number of significant axes
(nb.ax <- length(which(axes.test[,ncol(axes.test)] <= 0.05)))

## 5. Plot the significant canonical axes
mite.rda2.axes <- scores(mite.dbmem.rda2, choices=c(1:nb.ax), display="lc", 
	scaling=1)
dev.new(title="dbMEM analysis of mite data", width=8, height=6)
par(mfrow=c(1,nb.ax))
for(i in 1:nb.ax){
sr.value(mite.xy, mite.rda2.axes[,i], sub=paste("RDA",i), csub=2)
                 }

# Interpreting the spatial variation: regression of the significant
# canonical axes on the environmental variables
shapiro.test(resid(lm(mite.rda2.axes[,1] ~ ., data=mite.env))) # Normality test
mite.rda2.axis1.env <- lm(mite.rda2.axes[,1]~., data=mite.env)
summary(mite.rda2.axis1.env)
shapiro.test(resid(lm(mite.rda2.axes[,2] ~ ., data=mite.env))) # Normality test
mite.rda2.axis2.env <- lm(mite.rda2.axes[,2] ~ ., data=mite.env)
summary(mite.rda2.axis2.env)
shapiro.test(resid(lm(mite.rda2.axes[,3] ~ ., data=mite.env))) # Normality test
mite.rda2.axis3.env <- lm(mite.rda2.axes[,3] ~ ., data=mite.env)
summary(mite.rda2.axis3.env)

# Maps of the 8 significant dbMEM variables
# with homemade function sr.value()
dev.new(title="8 dbMEM variables - mites")
par(mfrow=c(2,4))
for(i in 1:ncol(dbmem.red))
{
	sr.value(mite.xy, dbmem.red[,i], sub=dbMEM.red[i], csub=2)
}


# dbMEM analysis of the mite data - broad scale
# *********************************************

(mite.dbmem.broad <- rda(mite.h.det ~ ., data=mite.dbmem[,c(1,3,4)]))
anova(mite.dbmem.broad)
(axes.broad <- anova(mite.dbmem.broad, by="axis"))
# Number of significant axes
(nb.ax.broad <- length(which(axes.broad[,ncol(axes.broad)] <= 0.05)))

# Plot of the two significant canonical axes
mite.dbmembroad.axes <- scores(mite.dbmem.broad, choices=c(1,2), 
	display="lc", scaling=1)
dev.new(title="dbMEM analysis of mite data - broad scale")
par(mfrow=c(1,2))
sr.value(mite.xy, mite.dbmembroad.axes[,1])
sr.value(mite.xy, mite.dbmembroad.axes[,2])

# Interpreting spatial variation: regression of the two 
# significant spatial canonical axes on the environmental variables
mite.dbmembroad.ax1.env <- lm(mite.dbmembroad.axes[,1] ~ ., data=mite.env)
summary(mite.dbmembroad.ax1.env)
mite.dbmembroad.ax2.env <- lm(mite.dbmembroad.axes[,2] ~ ., data=mite.env)
summary(mite.dbmembroad.ax2.env)


# dbMEM analysis of the mite data - medium scale
# **********************************************

(mite.dbmem.med <- rda(mite.h.det ~ ., data=mite.dbmem[,c(6,7,10,11)]))
anova(mite.dbmem.med)
(axes.med <- anova(mite.dbmem.med, by="axis"))
# Number of significant axes
(nb.ax.med <- length(which(axes.med[,ncol(axes.med)] <= 0.05)))

# Plot of the significant canonical axes
mite.dbmemmed.axes <- scores(mite.dbmem.med, choices=c(1,2), display="lc", 
	scaling=1)
dev.new(title="dbMEM analysis of mite data - medium scale")
par(mfrow=c(1,2))
sr.value(mite.xy, mite.dbmemmed.axes[,1])
sr.value(mite.xy, mite.dbmemmed.axes[,2])

# Interpreting spatial variation: regression of the significant 
# spatial canonical axes on the environmental variables
mite.dbmemmed.ax1.env <- lm(mite.dbmemmed.axes[,1] ~ ., data=mite.env)
summary(mite.dbmemmed.ax1.env)
mite.dbmemmed.ax2.env <- lm(mite.dbmemmed.axes[,2] ~ ., data=mite.env)
summary(mite.dbmemmed.ax2.env)


# dbMEM analysis of the mite data - fine scale
# ********************************************

(mite.dbmem.fine <- rda(mite.h.det ~ ., data=as.data.frame(mite.dbmem[,20])))
anova(mite.dbmem.fine)
# Analysis stops here, since the RDA is not significant.


# Single-step dbMEM analysis using function quickMEM()
# ****************************************************

# This function replaces the old quickPCNM()

dev.new(title="One-step dbMEM analysis of mite data")
mite.dbmem.quick <- quickMEM(mite.h, mite.xy)
summary(mite.dbmem.quick)
# Eigenvalues
mite.dbmem.quick[[2]]  # OR mite.dbmem.quick$eigenvalues
# Results of forward selection
mite.dbmem.quick[[3]]  # OR mite.dbmem.quick$fwd.sel

# Extract and plot RDA results from a quickMEM output (scaling 2)
dev.new(title="Biplot of dbMEM analysis result - scaling 2")
plot(mite.dbmem.quick$RDA, scaling=2)
sp.scores2 <- scores(mite.dbmem.quick$RDA, choices=1:2, scaling=2, display="sp")
arrows(0, 0, sp.scores2[,1]*0.9, sp.scores2[,2]*0.9, length=0, lty=1, col="red")


# Mite - trend - environment - dbMEM variation partitioning
# *********************************************************

## 1. Test trend. If significant, forward-select coordinates
mite.XY.rda <- rda(mite.h, mite.xy)
anova(mite.XY.rda)
(mite.XY.R2a <- RsquareAdj(mite.XY.rda)$adj.r.squared)
(mite.XY.fwd <- forward.sel(mite.h, as.matrix(mite.xy), 
	adjR2thresh=mite.XY.R2a))
XY.sign <- sort(mite.XY.fwd$order)
# Write the significant coordinates to a new object
XY.red <- as.data.frame(mite.xy[,c(XY.sign)])

## 2. Test and forward selection of environmental variables
# Recode environmental variables 3 to 5 into dummy binary variables
substrate <- model.matrix(~mite.env[,3])[,-1]
shrubs <- model.matrix(~mite.env[,4])[,-1]
topography <- model.matrix(~mite.env[,5])[,-1]
mite.env2 <- cbind(mite.env[,1:2], substrate, shrubs, topography)
colnames(mite.env2) <- c("SubsDens", "WatrCont", "Interface", "Litter", "Sphagn1"," Sphagn2", "Sphagn3", "Sphagn4", "Shrubs_Many", "Shrubs_None", "topography")
# Forward selection of the environmental variables
mite.env.rda <- rda(mite.h, mite.env2)
(mite.env.R2a <- RsquareAdj(mite.env.rda)$adj.r.squared)
mite.env.fwd <- forward.sel(mite.h, mite.env2, adjR2thresh=mite.env.R2a,
	nperm=9999)
env.sign <- sort(mite.env.fwd$order)
env.red <- mite.env2[,c(env.sign)]
colnames(env.red)

## 3. Test and forward selection of dbMEM variables
# Run the global dbMEM analysis on the *undetrended* mite data
mite.undet.dbmem.rda <- rda(mite.h, mite.dbmem)
anova(mite.undet.dbmem.rda)
# Since the analysis is significant, compute the adjusted R2
# and run a forward selection of the dbMEM variables
(mite.undet.dbmem.R2a <- RsquareAdj(mite.undet.dbmem.rda)$adj.r.squared)
(mite.undet.dbmem.fwd <- forward.sel(mite.h, as.matrix(mite.dbmem), 
	adjR2thresh=mite.undet.dbmem.R2a))
# Number of significant dbMEM
(nb.sig.dbmem <- nrow(mite.undet.dbmem.fwd))
# Identity of significant dbMEM in increasing order
(dbmem.sign <- sort(mite.undet.dbmem.fwd$order))
# Write the significant dbMEM to a new object
dbmem.red <- mite.dbmem[,c(dbmem.sign)]

## 4. Arbitrary split of the significant dbMEM into broad and fine scale
# Broad scale: dbMEM 1, 2, 3, 4, 6, 7, 8, 9
dbmem.broad <- dbmem.red[,1:8]
# Fine scale: dbMEM 16, 20
dbmem.fine <- dbmem.red[,9:10]

## 5. Mite - environment - trend - dbMEM variation partitioning
(mite.varpart <- varpart(mite.h, env.red, XY.red, dbmem.broad, dbmem.fine))
dev.new(title="Mite - environment - dbMEM variation partitioning", width=12, height=6)
par(mfrow=c(1,2))
showvarparts(4) # To show the symbols of the fractions
plot(mite.varpart, digits=2)
# Tests of the unique fractions [a], [b], [c] and [d]
# Fraction [a], pure environmental
anova(rda(mite.h, env.red, cbind(XY.red, dbmem.broad, dbmem.fine)))
# Fraction [b], pure trend
anova(rda(mite.h, XY.red, cbind(env.red, dbmem.broad, dbmem.fine)))
# Fraction [c], pure broad scale spatial
anova(rda(mite.h, dbmem.broad, cbind(env.red, XY.red, dbmem.fine)))
# Fraction [d], pure fine scale spatial
anova(rda(mite.h, dbmem.fine, cbind(env.red, XY.red, dbmem.broad)))



# MEM analysis of the detrended oribatid mite data
# ************************************************

# METHODOLOGICAL UPDATE: selection of MEM with positive autocorrelation only

## Selection of an optimal spatial weighting matrix

# 1. Search based on Delaunay triangulation.
#    We use mite.h.det as response data and mite.del as Delaunay
#    triangulation data.
#    No weighting matrix (binary weights only); function test.W selects among the
#    MEM variables constructed on the basis of the Delaunay triangulation.
# Delaunay triangulation and model selection
(mite.del <- tri2nb(mite.xy)) 
mite.del.res <- test.W(mite.h.det, mite.del, MEM.autocor="positive")
# Summary of the results for the best model
summary(mite.del.res$best)
# Unadjusted R^2 of best model
# This line returns the R^2 of the model with the smallest AICc value
(R2.del <- mite.del.res$best$AIC$R2[which.min(mite.del.res$best$AIC$AICc)])
# Adjusted R^2 of best model
RsquareAdj(R2.del, n=nrow(mite.h.det), m=which.min(mite.del.res$best$AIC$AICc))

# 2. Delaunay triangulation weighted by a function of distance.
#    Distances are ranged to maximum 1, and raised to power y
f2 <- function(D, dmax, y) { 1 - (D/dmax)^y }
# Largest Euclidean distance on links belonging to the Delaunay 
# triangulation
max.d1 <- max(unlist(nbdists(mite.del, as.matrix(mite.xy)))) 
# Power y is set from 2 to 10
mite.del.f2 <- test.W(mite.h.det, mite.del, MEM.autocor="positive", f=f2, 
    y=2:10, dmax=max.d1, xy=as.matrix(mite.xy))
# Unadjusted R^2 of best model
(R2.delW <- mite.del.f2$best$AIC$R2[which.min(mite.del.f2$best$AIC$AICc)])
# Adjusted R^2 of best model
RsquareAdj(R2.delW, n=nrow(mite.h.det), m=which.min(mite.del.f2$best$AIC$AICc))

# 3a. Connectivity matrix based on a distance (radius around points)
# Assessment of the relevant distances, based on a multivariate 
# variogram of the detrended mite data, with 20 distance classes.
(mite.vario <- variogmultiv(mite.h.det, mite.xy, nclass=20))
dev.new(title="Multivariate variogram, mites", width=12, height=6)
plot(mite.vario$d, mite.vario$var, ty='b', pch=20, xlab="Distance", 
	ylab="C(distance)")
# Construction of 10 neighbourhood matrices (class nb)
# Vector of 10 threshold distances
(thresh10 <- seq(give.thresh(dist(mite.xy)), 4, le=10))
# Create 10 neighbourhood matrices.
# Each matrix contains all connexions with lengths ² the threshold value
list10nb <- lapply(thresh10, dnearneigh, x=as.matrix(mite.xy), d1=0)
# Display an excerpt of the first neighbourhood matrix
print(listw2mat(nb2listw(list10nb[[1]], style="B"))[1:10,1:10], digits=1)
# Now we can apply the function test.W() to the 10 neighbourhood matrices.
# There are no weights on the links.
mite.thresh.res <- lapply(list10nb, 
    function(x) test.W(x, Y=mite.h.det, MEM.autocor="positive"))
# Lowest AICc, best model, threshold distance of best model
mite.thresh.minAIC <- sapply(mite.thresh.res, function(x) min(x$best$AIC$AICc, 
	na.rm=TRUE))
# Smallest AICc (best model among the 10)
min(mite.thresh.minAIC)
# Number of the model among the 10
which.min(mite.thresh.minAIC)
# Truncation threshold (distance)
thresh10[which.min(mite.thresh.minAIC)]

# 3b. Variant: same as above, but connections weighted by the complement
#     of the power of the distances, 1-(d/dmax)^y
mite.thresh.f2 <- lapply(list10nb, function(x) test.W(x, Y=mite.h.det, 
    MEM.autocor="positive", f=f2, y=2:10, dmax=max(unlist(nbdists(x, 
	as.matrix(mite.xy)))), xy=as.matrix(mite.xy)))
# Lowest AIC, best model
mite.f2.minAIC <- sapply(mite.thresh.f2, function(x) min(x$best$AIC$AICc, 
	na.rm=TRUE))
# Smallest AICc (best model among the 10)
min(mite.f2.minAIC)
# Number of the model among the 10
(nb.bestmod <- which.min(mite.f2.minAIC))
# Actual dmax of best model
(dmax.best <- mite.thresh.f2[nb.bestmod][[1]]$all[1,2])

# Extraction of the champion MEM model
mite.MEM.champ <- unlist(mite.thresh.f2[which.min(mite.f2.minAIC)], 
	recursive=FALSE)
summary(mite.MEM.champ)
# Number of MEM variables in best model
(nvars.best <- which.min(mite.MEM.champ$best$AIC$AICc))
# Eigenvalues
# mite.MEM.champ$best$AIC$values # No longer in output object
# MEM variables by order of added R2
mite.MEM.champ$best$AIC$ord
# MEM variables selected in the best model
MEMid <- mite.MEM.champ$best$AIC$ord[1:nvars.best]
sort(MEMid)
MEM.all <- mite.MEM.champ$best$MEM
MEM.select <- mite.MEM.champ$best$MEM[, sort(c(MEMid))]
colnames(MEM.select) <- sort(MEMid)
# Unadjusted R2 of best model
R2.MEMbest <- mite.MEM.champ$best$AIC$R2[nvars.best]
# Adjusted R2 of best model
RsquareAdj(R2.MEMbest, nrow(mite.h.det), length(MEMid))
# Plot the links using the function plot.links()
dev.new(title="Links, mite MEM champion model")
plot.links(mite.xy, thresh=dmax.best)

# RDA of the mite data constrained by the significant MEM, using vegan
(mite.MEM.rda <- rda(mite.h.det~., as.data.frame(MEM.select)))
(mite.MEM.R2a <- RsquareAdj(mite.MEM.rda)$adj.r.squared)
anova(mite.MEM.rda)
(axes.MEM.test <- anova(mite.MEM.rda, by="axis"))
# Number of significant axes
(nb.ax <- length(which(axes.MEM.test[,ncol(axes.MEM.test)] <= 0.05)))

# Plot maps of the significant canonical axes
mite.MEM.axes <- scores(mite.MEM.rda, choices=1:nb.ax, display="lc", 
	scaling=1)
dev.new(title="MEM analysis of mite data")
if(nb.ax <= 2) {
   par(mfrow=c(1,2))
   } else { par(mfrow=c(2,2)) }
for(i in 1:ncol(mite.MEM.axes)){
sr.value(mite.xy, mite.MEM.axes[,i])
}

# Maps of the significant MEM variables
dev.new(title="Maps of significant MEM eigenfunctions")
if(ncol(MEM.select) <= 6) {
   par(mfrow=c(2,3))
   } else { par(mfrow=c(3,3)) }
for(i in 1:ncol(MEM.select)){
	sr.value(mite.xy, MEM.select[,i], sub=sort(MEMid)[i], csub=2)
}

# Correlation of the retained MEM and dbMEM variables
cor(MEM.select, dbmem.red)


# Other connectivity matrices
# ***************************

# Examples of connectivity matrices in decreasing order of connectivity
# All these neighbourhood matrices are stored in objects of class nb
# Delaunay triangulation (as in the previous example)
mite.del <- tri2nb(mite.xy) 
# Gabriel graph
mite.gab <- graph2nb(gabrielneigh(as.matrix(mite.xy)), sym=TRUE)
# Relative neighbourhood
mite.rel <- graph2nb(relativeneigh(as.matrix(mite.xy)), sym=TRUE)
# Minimum spanning tree
mite.mst <- mst.nb(dist(mite.xy))

# Plots of the connectivity matrices
dev.new(title="Connectivity matrices", width=6, height=9)
par(mfrow=c(2,2))
plot(mite.del, mite.xy, col="red", pch=20, cex=1)
title(main="Delaunay triangulation ")
plot(mite.gab, mite.xy, col="purple", pch=20, cex=1)
title(main="Gabriel graph")
plot(mite.rel, mite.xy, col="dark green", pch=20, cex=1)
title(main="Relative neighbourhood")
plot(mite.mst, mite.xy, col="brown", pch=20, cex=1)
title(main="Minimum spanning tree")


# Link editing
# ************

# 1. Interactive:
dev.new(title="Delaunay triangulation")
plot(mite.del, mite.xy, col="red", pch=20, cex=2)
title(main="Delaunay triangulation")
mite.del2 <- edit.nb(mite.del, mite.xy)
	# To delete a link, click on its two nodes. Follow on-screen instructions.
	# Wait until you have finished editing before entering the next 
	# command line. Suggestion: edit a link to site 23 (see below why).

# 2. Alternatively, links can also be removed by command lines, 
# after having converted the nb object into an editable matrix:
mite.del.mat <- nb2mat(mite.del, style="B")
# Remove connection between objects 23 and 35:
mite.del.mat[23,35] <- 0
mite.del.mat[35,23] <- 0
# Back-conversion into nb object:
mite.del3 <- neig2nb(neig(mat01=mite.del.mat))
dev.new(title="Delaunay with edited links")
plot(mite.del3, mite.xy)

# Example: list of neighbours of site 23 for the Delaunay triangulation:
mite.del[[23]]		# Before editing
mite.del2[[23]]		# After interactive editing - depends on what you have edited above.
mite.del3[[23]]		# After command line editing


# Connectivity matrix based on a distance (radius around points)
# Using the same truncation distance dmin as in the dbMEM example (1.011187).
dmin = 1.011187
mite.thresh4 <- dnearneigh(as.matrix(mite.xy), 0, dmin*4)
# Display some values
nb2mat(mite.thresh4)[1:10,1:10]

# Using a shorter distance (1*dmin, 2*dmin)
mite.thresh1 <- dnearneigh(as.matrix(mite.xy), 0, dmin*1)
mite.thresh2 <- dnearneigh(as.matrix(mite.xy), 0, dmin*2)
# Using a longer distance
mite.thresh8 <- dnearneigh(as.matrix(mite.xy), 0, dmin*8)

# Plot of some connectivity matrices
dev.new(title="Connectivity matrices - threshold distances")
par(mfrow=c(1,2))
plot(mite.thresh1, mite.xy, col="red", pch=20, cex=0.8)
title(main="1 * dmin")
plot(mite.thresh4, mite.xy, col="red", pch=20, cex=0.8)
title(main="4 * dmin")

# Conversion of a "nb" object into a "listw" object
# Example: mite.thresh4 created above. "B" is for "binary"
mite.thresh4.lw <- nb2listw(mite.thresh4, style="B")
print(listw2mat(mite.thresh4.lw)[1:10,1:10], digits=1)

# Creation of a spatial weighting matrix W = Hadamard product of B and A
# Replace "1" by Euclidean distances in the connectivity matrix
mite.thresh4.d1 <- nbdists(mite.thresh4, as.matrix(mite.xy))
# Weights as function of inverse distance
mite.inv.dist <- lapply(mite.thresh4.d1, function(x) 1-x/max(dist(mite.xy)))
# Creation of spatial weighting matrix W. Argument "B" stands for 
# "binary" but concerns the links themselves, not their weights
mite.invdist.lw <- nb2listw(mite.thresh4, glist=mite.inv.dist, style="B")
print(listw2mat(mite.invdist.lw)[1:10,1:10], digits=2)

# Computation of MEM variables (from an object of class listw)
mite.invdist.MEM <- scores.listw(mite.invdist.lw)
summary(mite.invdist.MEM)
attributes(mite.invdist.MEM)$values
dev.new(title="Barplot of MEM eigenvalues")
barplot(attributes(mite.invdist.MEM)$values)

# Store all MEM vectors in new object
mite.invdist.MEM.vec <- as.matrix(mite.invdist.MEM)

# Test of Moran's I of each eigenvector
mite.MEM.Moran <- moran.randtest(mite.invdist.MEM.vec, mite.invdist.lw, 999)

# MEM with significant spatial correlation
which(mite.MEM.Moran$pvalue <= 0.05)
length(which(mite.MEM.Moran$pvalue <= 0.05))

# MEM with positive spatial correlation
MEM.Moran.pos <- which(mite.MEM.Moran$obs > -1/(nrow(mite.invdist.MEM.vec)-1))
mite.invdist.MEM.pos <- mite.invdist.MEM.vec[,MEM.Moran.pos]
# MEM with positive *and significant* spatial correlation
mite.invdist.MEM.pos.sig <- mite.invdist.MEM.pos[,which(mite.MEM.Moran$pvalue <= 0.05)]

# Plot of MEM eigenvalues vs Moran's I
dev.new(title="MEM eigenvalues vs Moran's I")
plot(attributes(mite.invdist.MEM)$values, mite.MEM.Moran$obs, ylab="Moran's I", 
	xlab="Eigenvalues")
text(0, 0.55, paste("Correlation=", cor(mite.MEM.Moran$obs, attributes(mite.invdist.MEM)$values)))



# AEM analysis
# ************

# Coding of a river arborescence. 
# See Legendre and Legendre (2012, p. 889).
node1 <- c(1,0,0,0,0,0,0,0)
n2lk6 <- c(0,1,0,0,0,0,0,0)
n3lk3 <- c(1,0,1,0,0,0,0,0)
n4lk2 <- c(1,0,0,1,0,0,0,0)
node5 <- c(0,1,0,0,1,0,0,0)
n6lk1 <- c(1,0,0,1,0,1,0,0)
ln7k4 <- c(0,1,0,0,1,0,1,0)
n8lk5 <- c(0,1,0,0,1,0,0,1)
arbor <- rbind(node1, n2lk6, n3lk3, n4lk2, node5, n6lk1, ln7k4, n8lk5)

# AEM construction
(arbor.aem <- aem(binary.mat=arbor))
arbor.aem.vec <- arbor.aem$vectors

# AEM eigenfunctions can also be obtained directly by singular value 
# decomposition (function svd()), which is what the function aem() does:
arbor.c <- scale(arbor, center=TRUE, scale=FALSE)
arbor.svd <- svd(arbor.c)
# Singular values of the construction bove
arbor.svd$d[1:7]
# AEM eigenfunctions of the construction above
arbor.svd$u[,1:7]

# Coding of sampling design: 10 cross-river transects, 4 traps 
# per transect. Edges weighted proportional to inverse squared distance.
# X-Y coordinates
xy <- cbind(1:40, expand.grid(1:4, 1:10))
# Object of class nb (spdep) containing links of chess type "queen"
nb <- cell2nb(4, 10, "queen")
# Site-by-edges matrix (produces a fictitious object "0")
edge.mat <- build.binary(nb, xy)
# Matrix of Euclidean distances
D1.mat <- as.matrix(dist(xy))
# Extract the edges, remove the ones directly linked to site 0
edges.b <- edge.mat$edges[-1:-4,]
# Construct a vector giving the length of each edge
length.edge <- vector(length=nrow(edges.b))
for(i in 1:nrow(edges.b))
{
	length.edge[i] <- D1.mat[edges.b[i,1], edges.b[i,2]]
}
# Weighting of edges based on inverse squared distance
weight.vec <- 1-(length.edge/max(length.edge))^2
# Construction of AEM eigenfunctions from edge.mat, of class build.binary
example.AEM <- aem(build.binary=edge.mat, weight=weight.vec, rm.link0=TRUE)
example.AEM$values
ex.AEM.vec <- example.AEM$vectors

# Construction of 5 fictitious species
# Two randomly distributed species
sp12 <- matrix(trunc(rnorm(80,5,2),0),40)
# One species restricted to the upper half of the stream
sp3 <- c(trunc(rnorm(20,8,2.5),0), rep(0,20))
# One species restricted to the left-hand half of the transect
sp4 <- t(matrix(c(trunc(rnorm(20,8,3),0), rep(0,20)),10))
sp4 <- as.vector(sp4)
# One species restricted to the 4 upper left-hand sites
sp5 <- c(4,7,0,0,3,8, rep(0,34))
# Build the species matrix
sp <- cbind(sp12, sp3, sp4, sp5)
colnames(sp) <- c("sp1", "sp2", "sp3", "sp4", "sp5")

# Global AEM analysis with 20 first AEM variables (for computation of R2a)
AEM.20 <- rda(sp ~ ., as.data.frame(ex.AEM.vec[,1:20]))
(R2a.AEM <- RsquareAdj(AEM.20)$adj.r.squared)

# Forward selection of the AEM variables
AEM.fwd <- forward.sel(sp, ex.AEM.vec, adjR2thresh=R2a.AEM)
(AEM.sign <- sort(AEM.fwd[,2]))
# Write significant AEM in a new object
AEM.sign.vec <- ex.AEM.vec[,c(AEM.sign)]
# RDA with signif. AEM
(sp.AEMsign.rda <- rda(sp ~ ., data=as.data.frame(AEM.sign.vec)))
anova(sp.AEMsign.rda)
(AEM.rda.axes.test <- anova(sp.AEMsign.rda, by="axis"))
# Number of significant axes
(nb.ax.AEM <- length(which(AEM.rda.axes.test[,4] <= 0.05)))
# Adjusted R-square
RsquareAdj(sp.AEMsign.rda)

# Plot of the significant canonical axes
AEM.rda.axes <- scores(sp.AEMsign.rda, choices=c(1,2), display="lc", scaling=1)
dev.new(title="AEM analysis of fictitious data, significant RDA axes")
par(mfrow=c(1,nb.ax.AEM))
for(i in 1:nb.ax.AEM) sr.value(xy[,c(2,3)], AEM.rda.axes[,i])


# Multiscale ordination (MSO)
# ***************************

# MSO of the undetrended mite data vs environment RDA
mite.undet.env.rda <- rda(mite.h~., mite.env2)
(mite.env.rda.mso <- mso(mite.undet.env.rda, mite.xy, grain=dmin, perm=999))
dev.new(title="MSO plot of the undetrended mite-environment RDA")
msoplot(mite.env.rda.mso, alpha=0.05/7)

# MSO of the undetrended mite data vs environment RDA, controlling for MEM
mite.undet.env.MEM <- rda(mite.h, mite.env2, as.data.frame(MEM.select))
(mite.env.MEM.mso <- mso(mite.undet.env.MEM, mite.xy, grain=dmin, perm=999))
dev.new(title="MSO plot of the undetrended mite-environment RDA controlling for MEM")
# msoplot(mite.env.MEM.mso, alpha=0.05/7)
msoplot(mite.env.MEM.mso, alpha=0.05/7, ylim=c(0, 0.0045)) # Expanded height to clear legend

# MSO on detrended mite and environmental data
# Detrend mite data on Y coordinate
mite.h.det2 <- resid(lm(as.matrix(mite.h) ~ mite.xy[,2]))
# Detrend environmental data on Y coordinate
env2.det <- resid(lm(as.matrix(mite.env2) ~ mite.xy[,2]))
# RDA and MSO
mitedet.envdet.rda <- rda(mite.h.det2, env2.det)
(miteenvdet.rda.mso <- mso(mitedet.envdet.rda, mite.xy, grain=dmin, perm=999))
dev.new(title="MSO plot of the detrended mite-environment RDA")
msoplot(miteenvdet.rda.mso, alpha=0.05/7, ylim=c(0, 0.006))

# MSO of the detrended mite data vs environment RDA, controlling for MEM
mite.det.env.MEM <- rda(mite.h.det2, env2.det, as.data.frame(MEM.select))
(mite.env.MEM.mso <- mso(mite.det.env.MEM, mite.xy, grain=dmin, perm=999))
dev.new(title="MSO plot of the detrended mite-environment RDA controlling for MEM")
msoplot(mite.env.MEM.mso, alpha=0.05/7, ylim=c(0, 0.005))
