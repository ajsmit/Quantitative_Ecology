################################################################################
### CHAPTER 4: CLUSTER ANALYSIS
###
### Online supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 30 January 2017
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(gclus)
library(cluster)
library(RColorBrewer)
library(labdsv)
library(mvpart)
library(MVPARTwrap)

# Load additionnal functions
# (files must be in the working directory)
source("Num_Ecol_R_book_ed1/hcoplot.R")
source("Num_Ecol_R_book_ed1/test.a.R")
source("Num_Ecol_R_book_ed1/coldiss.R")

# Import the data from CSV files
# (files must be in the working directory)
spe <- read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names=1)
env <- read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names=1)
spa <- read.csv("Num_Ecol_R_book_ed1/DoubsSpa.csv", row.names=1)
# Remove empty site 8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]



# Hierarchical agglomerative clustering of the species abundance data
# *******************************************************************

# Compute matrix of chord distance among sites followed by
# single linkage agglomerative clustering
spe.norm <- decostand(spe, "normalize")
spe.ch <- vegdist(spe.norm, "euc")
spe.ch.single <- hclust(spe.ch, method="single")
# Plot a dendrogram using the default options
dev.new(title="Fish - Chord - Single linkage", width=12, height=8)
plot(spe.ch.single)

# Compute complete-linkage agglomerative clustering
spe.ch.complete <- hclust(spe.ch, method="complete")
dev.new(title="Fish - Chord - Complete linkage", width=12, height=8)
plot(spe.ch.complete)

# Compute UPGMA agglomerative clustering
spe.ch.UPGMA <- hclust(spe.ch, method="average")
dev.new(title="Fish - Chord - UPGMA", width=12, height=8)
plot(spe.ch.UPGMA)

# Compute centroid clustering of the fish data
spe.ch.centroid <- hclust(spe.ch, method="centroid")
dev.new(title="Fish - Chord - Centroid", width=12, height=8)
plot(spe.ch.centroid)

# Compute Ward's minimum variance clustering
spe.ch.ward <- hclust(spe.ch, method="ward.D2")
# Note: in R 3.0.3 or older type:
#    spe.ch.ward <- hclust(spe.ch, method="ward")
#    spe.ch.ward$height <- sqrt(spe.ch.ward)
# The result will not be striclty the same, however.
# See the hclust help file and 'Numerical Ecology with R',
# 1st edition, p. 61-62.
dev.new(title="Fish - Chord - Ward", width=12, height=8)
plot(spe.ch.ward)


# Cophenetic correlations
# ***********************

# Single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)
# Complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)
# Average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)
# Ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)
cor(spe.ch, spe.ch.ward.coph, method="spearman")

# Shepard-like diagrams
dev.new(title="Cophenetic correlation", width=8, height=8)
par(mfrow=c(2,2))
plot(spe.ch, spe.ch.single.coph, xlab="Chord distance",
	ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
	main=c("Single linkage", paste("Cophenetic correlation =",
	round(cor(spe.ch, spe.ch.single.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.single.coph), col="red")
plot(spe.ch, spe.ch.comp.coph, xlab="Chord distance",
	ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
	main=c("Complete linkage", paste("Cophenetic correlation =",
	round(cor(spe.ch, spe.ch.comp.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.comp.coph), col="red")
plot(spe.ch, spe.ch.UPGMA.coph, xlab="Chord distance",
	ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
	main=c("UPGMA", paste("Cophenetic correlation =",
	round(cor(spe.ch, spe.ch.UPGMA.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.UPGMA.coph), col="red")
plot(spe.ch, spe.ch.ward.coph, xlab="Chord distance",
	ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)),
	ylim=c(0,max(spe.ch.ward$height)),
	main=c("Ward clustering", paste("Cophenetic correlation =",
	round(cor(spe.ch, spe.ch.ward.coph),3))))
abline(0,1)
lines(lowess(spe.ch, spe.ch.ward.coph), col="red")

# Gower (1983) distance
(gow.dist.single <- sum((spe.ch-spe.ch.single.coph)^2))
(gow.dist.comp <- sum((spe.ch-spe.ch.comp.coph)^2))
(gow.dist.UPGMA <- sum((spe.ch-spe.ch.UPGMA.coph)^2))
(gow.dist.ward <- sum((spe.ch-spe.ch.ward.coph)^2))


# Graphs of fusion level values
# *****************************

dev.new(title="Fusion levels", width=12, height=8)
par(mfrow=c(2,2))
# Plot the fusion level values of the single linkage clustering
plot(spe.ch.single$height, nrow(spe):2, type="S",
	main="Fusion levels - Chord - Single",
	ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(spe.ch.single$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
# Plot the fusion level values of the complete linkage clustering
plot(spe.ch.complete$height, nrow(spe):2, type="S",
	main="Fusion levels - Chord - Complete",
	ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(spe.ch.complete$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
# Plot the fusion level values of the UPGMA clustering
plot(spe.ch.UPGMA$height, nrow(spe):2, type="S",
	main="Fusion levels - Chord - UPGMA",
	ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(spe.ch.UPGMA$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)
# Plot the fusion level values of the Ward clustering
plot(spe.ch.ward$height, nrow(spe):2, type="S",
	main="Fusion levels - Chord - Ward",
	ylab="k (number of clusters)", xlab="h (node height)", col="grey")
text(spe.ch.ward$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)


# Cut the trees to obtain k groups and compare the group contents
# using contingency tables
# ***************************************************************

# Choose a common number of groups
k <- 4  # Number of groups where at least a small jump is present
        # in all four graphs of fusion levels
# Cut the dendrograms
spebc.single.g <- cutree(spe.ch.single, k)
spebc.complete.g <- cutree(spe.ch.complete, k)
spebc.UPGMA.g <- cutree(spe.ch.UPGMA, k)
spebc.ward.g <- cutree(spe.ch.ward, k)

# Compare classifications by constructing contingency tables
# Single vs complete linkage
table(spebc.single.g, spebc.complete.g)
# Single linkage vs UPGMA
table(spebc.single.g, spebc.UPGMA.g)
# Single linkage vs Ward
table(spebc.single.g, spebc.ward.g)
# Complete linkage vs UPGMA
table(spebc.complete.g, spebc.UPGMA.g)
# Complete linkage vs Ward
table(spebc.complete.g, spebc.ward.g)
# UPGMA vs Ward
table(spebc.UPGMA.g, spebc.ward.g)


# Optimal number of clusters according to silhouette widths
# (Rousseeuw quality index)
# *********************************************************

# Plot average silhouette widths (using Ward clustering) for all partitions
# except for the trivial partition in a single group (k=1)
# First, create an empty vector in which the average silhouette width (asw) values will be written
asw <- numeric(nrow(spe))
# Second, retrieve the asw values and write them into the vector
for (k in 2:(nrow(spe)-1)) {
	sil <- silhouette(cutree(spe.ch.ward, k=k), spe.ch)
	asw[k] <- summary(sil)$avg.width
	}
k.best <- which.max(asw)
dev.new(title="Silhouettes - Ward - k = 2 to n-1")
plot(1:nrow(spe), asw, type="h",
	main="Silhouette-optimal number of clusters, Ward",
	xlab="k (number of groups)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
	col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n",
	"with an average silhouette width of", max(asw), "\n")


# Optimal number of clusters according to Mantel statistic (Pearson)
# ******************************************************************

# Function to compute a binary dissimilarity matrix from groups
grpdist <- function(X)
{
	require(cluster)
	gr <- as.data.frame(as.factor(X))
	distgr <- daisy(gr, "gower")
	distgr
}

# Run based on the Ward clustering
kt <- data.frame(k=1:nrow(spe), r=0)
for (i in 2:(nrow(spe)-1)) {
	gr <- cutree(spe.ch.ward, i)
	distgr <- grpdist(gr)
	mt <- cor(spe.ch, distgr, method="pearson")
	kt[i,2] <- mt
}
kt
k.best <- which.max(kt$r)
dev.new(title="Optimal number of clusters - Mantel")
plot(kt$k, kt$r, type="h", main="Mantel-optimal number of clusters - Ward",
	xlab="k (number of groups)", ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
	col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n",
	"with a matrix linear correlation of", max(kt$r), "\n")


# Silhouette plot of the final partition
# **************************************

# Choose the number of clusters
k <- 4
# Silhouette plot
cutg <- cutree(spe.ch.ward, k=k)
sil <- silhouette(cutg, spe.ch)
rownames(sil) <- row.names(spe)
dev.new(title="Silhouette plot - Ward - k=4")
plot(sil, main="Silhouette plot - Chord - Ward",
	cex.names=0.8, col=2:(k+1), nmax=100)


# Final dendrogram with the selected groups
# *****************************************

spe.chwo <- reorder.hclust(spe.ch.ward, spe.ch)

# Plot reordered dendrogram with group labels
dev.new(title="Final dendrogram", width=8, height=6)
plot(spe.chwo, hang=-1, xlab="4 groups", sub="",
	ylab="Height", main="Chord - Ward (reordered)",
	labels=cutree(spe.chwo, k=k))
rect.hclust(spe.chwo, k=k)

# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters,
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
hcoplot(spe.ch.ward, spe.ch, k=4)


# Plot of the 4 Ward clusters on a map of the Doubs river
# *******************************************************

# Map of the Doubs river (see Chapter 2)
dev.new(title="Four Ward clusters on river")
plot(spa, asp=1, type="n", main="Four Ward clusters along the Doubs river",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
text(65, 20, "Upstream", cex=1.2)
text(15, 32, "Downstream", cex=1.2)
# Add the four groups
grw <- spebc.ward.g
k <- length(levels(factor(grw)))
for (i in 1:k)
{
	points(spa[grw==i,1], spa[grw==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
}
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("Cluster", 1:k), pch=(1:k)+20, col=2:(k+1),
	pt.bg=2:(k+1), pt.cex=2, bty="n")


# Heat map
# ********

# Heat map of the dissimilarity matrix ordered with the dendrogram
dend <- as.dendrogram(spe.chwo)
dev.new(title="Heatmap - sites")
heatmap(as.matrix(spe.ch), Rowv=dend, symm=TRUE, margin=c(3,3))

# Ordered community table
# Species are ordered by their weighted averages on site scores
or <- vegemite(spe, spe.chwo)

# Heat map of the doubly ordered community table, with dendrogram
dev.new(title="Heatmap - species")
heatmap(t(spe[rev(or$species)]), Rowv=NA, Colv=dend,
	col=c("white", brewer.pal(5,"Greens")), scale="none", margin=c(4,4),
	ylab="Species (weighted averages of sites)", xlab="Sites")



# k-means partitioning of the pre-transformed species data
# ********************************************************

# With 4 groups
spe.kmeans <- kmeans(spe.norm, centers=4, nstart=100)

# Comparison with the 4-group classification derived from Ward clustering:
table(spe.kmeans$cluster, spebc.ward.g)

# k-means partitioning, 2 to 10 groups
spe.KM.cascade <- cascadeKM(spe.norm, inf.gr=2, sup.gr=10, iter=100,
	criterion="ssi")
summary(spe.KM.cascade)
spe.KM.cascade$results
spe.KM.cascade$partition
dev.new(title="CascadeKM", width=10, height=6)
plot(spe.KM.cascade, sortg=TRUE)

# Reorder the sites according to the k-means result
spe[order(spe.kmeans$cluster),]

# Reorder sites and species using function vegemite()
ord.KM <- vegemite(spe, spe.kmeans$cluster)
spe[ord.KM$sites, ord.KM$species]



# Partitioning around medoids (PAM)
# Computed on the chord distance matrix
# *************************************

# Choice of the number of clusters
# Loop to compute average silhouette width for 2 to 28 clusters.
asw <- numeric(nrow(spe))
for (k in 2:(nrow(spe)-1))
	asw[k] <- pam(spe.ch, k, diss=TRUE)$silinfo$avg.width
k.best <- which.max(asw)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n",
	"with an average silhouette width of", max(asw), "\n")
dev.new(title="PAM")
plot(1:nrow(spe), asw, type="h", main="Choice of the number of clusters",
	xlab="k (number of clusters)", ylab="Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2,
	col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)

# PAM for k = 4 clusters
spe.ch.pam <- pam(spe.ch, k=4, diss=TRUE)
summary(spe.ch.pam)
spe.ch.pam.g <- spe.ch.pam$clustering
spe.ch.pam$silinfo$widths

# Compare with classification from Ward clustering and from k-means
table(spe.ch.pam.g, spebc.ward.g)
table(spe.ch.pam.g, spe.kmeans$cluster)

# Silhouette profile for k = 4 groups, k-means and PAM
dev.new(title="Silhouettes - k-means and PAM", width=12, height=8)
par(mfrow=c(1,2))
k <- 4
sil <- silhouette(spe.kmeans$cluster, spe.ch)
rownames(sil) <- row.names(spe)
plot(sil, main="Silhouette plot - k-means",
	cex.names=0.8, col=2:(k+1))
plot(silhouette(spe.ch.pam), main="Silhouette plot - PAM", cex.names=0.8,
	col=2:(k+1))


# Plot of the 4 k-means clusters on a map of the Doubs river
# **********************************************************

dev.new(title="Four k-means groups on river")
plot(spa, asp=1, type="n", main="Four k-means groups",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
text(80, 10, "Upstream", cex=1.2)
text(15, 70, "Downstream", cex=1.2)
grKM <- spe.kmeans$cluster
k <- length(levels(factor(grKM)))
for (i in 1:k)
{
	points(spa[grKM==i,1], spa[grKM==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
}
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("Cluster",1:k), pch=(1:k)+20, col=2:(k+1),
	pt.bg=2:(k+1), pt.cex=2, bty="n")


# Relationships between fish clusters and 4 environmental variables
# based on the k-means clustering results (four groups)
# *****************************************************************

attach(env)

# Boxplots of quantitative environmental variables:
# Altitude, Slope, Oxygen, and Ammonium
dev.new(title="Boxplots of quantitative environmental variables")
par(mfrow=c(2,2))
boxplot(sqrt(alt) ~ spe.kmeans$cluster, main="Altitude", las=1,
	ylab="sqrt(alt)", col=2:5, varwidth=TRUE)
boxplot(log(slo) ~ spe.kmeans$cluster, main="Slope", las=1,
	ylab="log(slo)", col=2:5, varwidth=TRUE)
boxplot(oxy ~ spe.kmeans$cluster, main="Oxygen", las=1,
	ylab="oxy", col=2:5, varwidth=TRUE)
boxplot(sqrt(amm) ~ spe.kmeans$cluster, main="Ammonium", las=1,
	ylab="sqrt(amm)", col=2:5, varwidth=TRUE)

# Test of ANOVA assumptions
# Normality of residuals
shapiro.test(resid(aov(sqrt(alt) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(aov(log(slo) ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(aov(oxy ~ as.factor(spe.kmeans$cluster))))
shapiro.test(resid(aov(sqrt(amm) ~ as.factor(spe.kmeans$cluster))))

# Homogeneity of variances
bartlett.test(sqrt(alt), as.factor(spe.kmeans$cluster))
bartlett.test(log(slo), as.factor(spe.kmeans$cluster))
bartlett.test(oxy, as.factor(spe.kmeans$cluster))
bartlett.test(sqrt(amm), as.factor(spe.kmeans$cluster))

# ANOVA of the testable variables
summary(aov(log(slo) ~ as.factor(spe.kmeans$cluster)))
summary(aov(oxy ~ as.factor(spe.kmeans$cluster)))
summary(aov(sqrt(amm) ~ as.factor(spe.kmeans$cluster)))

# Kruskal-Wallis test of variable alt
kruskal.test(alt ~ as.factor(spe.kmeans$cluster))

detach(env)


# Contingency table of two typologies
# ***********************************

# Environment-based typology (see Chapter 2)
env2 <- env[,-1]
env.de <- vegdist(scale(env2), "euc")
env.kmeans <- kmeans(env.de, centers=4, nstart=100)
env.KM.4 <- env.kmeans$cluster
# Table crossing the k-means and environment 4-group typologies
table(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))

# Test the relationship using a chi-square test
chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster))
# Change the testing procedure to a permutation test
chisq.test(as.factor(spe.kmeans$cluster), as.factor(env.kmeans$cluster),
	simulate.p.value=TRUE)


# Mean abundances on k-means site clusters
# ****************************************

groups <- as.factor(spe.kmeans$cluster)
spe.means <- matrix(0, ncol(spe), length(levels(groups)))
row.names(spe.means) <- colnames(spe)
for(i in 1:ncol(spe))
{
	spe.means[i,] <- tapply(spe[,i], spe.kmeans$cluster, mean)
}
# Mean species abundances of the four groups
group1 <- round(sort(spe.means[,1], decreasing=TRUE), 2)
group2 <- round(sort(spe.means[,2], decreasing=TRUE), 2)
group3 <- round(sort(spe.means[,3], decreasing=TRUE), 2)
group4 <- round(sort(spe.means[,4], decreasing=TRUE), 2)
# Species with abundances greater than group mean species abundance
group1.domin <- which(group1 > mean(group1))
group1
group1.domin
#... same for other groups


# Kendall's W coefficient of concordance
# **************************************

# Extraction of the most abundant species
sp.sum <- apply(spe, 2, sum)
spe.sorted <- spe[,order(sp.sum, decreasing=TRUE)]
spe.small <- spe.sorted[,1:20]

# Transformation of species data and transposition
spe.small.hel <- decostand(spe.small, "hellinger")
spe.small.std <- decostand(spe.small.hel, "standardize")
spe.small.t <- t(spe.small.std)

# k-means partitioning of species
spe.t.kmeans.casc <- cascadeKM(spe.small.t, inf.gr=2, sup.gr=8,
	iter=100, criterion="calinski")
dev.new(title="Kmeans on species", width=10, height=6)
plot(spe.t.kmeans.casc, sortg=TRUE)

# The partition into 2 groups is found in column 1 of the object $partition
clusters <- spe.t.kmeans.casc$partition[,1]
clusters

# Concordance analysis
(spe.kendall.global <- kendall.global(spe.small.hel, clusters))

# A posteriori tests
(spe.kendall.post <- kendall.post(spe.small.hel, clusters, nperm=9999))


# Species assemblages on presence-absence data
# ********************************************

# Transform the data to presence-absence
spe.pa <- decostand(spe, "pa")
# Test the co-occurrence of species
# Many permutations to have enough power, but this may take some time!
res <- test.a(spe.pa, nperm=99999)
summary(res)

# Compute a Holm correction on the matrix of p-values
(res.p.vec <- as.vector(res$p.a.dist))
(adjust.res <- p.adjust(res.p.vec, method="holm"))
# Check that the number of permutations allows significant p-values
# after Holm's correction. Look for values <=0.05:
range(adjust.res)
# Significance threshold: among the corrected Holm p-values, find 0.05 or
# the closest value smaller than 0.05.
(adj.sigth <- max(adjust.res[adjust.res<=0.05]))
# Now find the uncorrected p-value corresponding to adj.sigth:
(sigth <- max(res.p.vec[adjust.res<=0.05]))
# Assign 1 to nonsignificant p-values after Holm's correction
res.pa.dist <- res$p.a.dist
res.pa.dist[res.pa.dist > sigth] <- 1
# How many (unadjusted) values are equal to or smaller than sigth?
length(which(res.p.vec<=sigth))

# Heat map of significant p-values
dev.new(title="0-1 species associations, heat map", width=14, height=7)
coldiss(res.pa.dist, nc=16, byrank=TRUE, diag=TRUE)

# Fuzzy clustering
res.pa.fuz <- fanny(res.pa.dist, k=5, memb.exp=1.5)
summary(res.pa.fuz)
dev.new(title="0-1 species associations, fuzzy clustering", width=8, height=8)
plot(silhouette(res.pa.fuz), main="Silhouette plot - Fuzzy clustering",
	cex.names=0.8, col=res.pa.fuz$silinfo$widths+1)


# Species indicator values (Dufrene and Legendre)
# ***********************************************

# Divide the sites into 4 groups depending on the distance from
# the source of the river
dfs.D1 <- dist(data.frame(dfs=env[,1], row.names=rownames(env)))
dfsD1.kmeans <- kmeans(dfs.D1, centers=4, nstart=100)
# Cluster delimitation and numbering
dfsD1.kmeans$cluster
# The numbering is machine-dependent. To avoid confusions let us construct a
# vector of groups with consecutive numbers on the basis of dfsD1.kmeans$cluster:
grps <- rep(1:4, c(8, 10, 6, 5))
# Indicator species for this typology of the sites
(iva <- indval(spe, grps))

# Table of the significant indicator species
gr <- iva$maxcls[iva$pval <= 0.05]
iv <- iva$indcls[iva$pval <= 0.05]
pv <- iva$pval[iva$pval <= 0.05]
fr <- apply(spe > 0, 2, sum)[iva$pval <= 0.05]
fidg <- data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg <- fidg[order(fidg$group, -fidg$indval),]
fidg
# Export the result to a CSV file (to be opened in a spreadsheet)
write.csv(fidg, "IndVal-dfs.csv")


# Multivariate regression trees
# *****************************

dev.new(title="Multivariate regression tree - all explanatory variables", width=14, height=7)
par(mfrow=c(1,2))
spe.ch.mvpart <- mvpart(data.matrix(spe.norm) ~ ., env, margin=0.08, cp=0,
	xv="pick", xval=nrow(spe), xvmult=100, which=4)
# Here, click on the desired number of groups (for example 4)

summary(spe.ch.mvpart)
printcp(spe.ch.mvpart)

# Residuals of MRT
dev.new(title="Residuals of MRT", width=10, height=6)
par(mfrow=c(1,2))
hist(residuals(spe.ch.mvpart), col="grey")
plot(predict(spe.ch.mvpart, type="matrix"), residuals(spe.ch.mvpart), # changed!
	main="Residuals vs Predicted")
abline(h=0, lty=3, col="grey")

# Group composition
spe.ch.mvpart$where
# Group identity
(groups.mrt <- levels(as.factor(spe.ch.mvpart$where)))
# Fish composition of first leaf
spe.norm[which(spe.ch.mvpart$where==groups.mrt[1]),]
# Environmental variables of first leaf
env[which(spe.ch.mvpart$where==groups.mrt[1]),]

# Table and pie charts of fish composition of leaves
leaf.sum <- matrix(0, length(groups.mrt), ncol(spe))
colnames(leaf.sum) <- colnames(spe)
for(i in 1:length(groups.mrt))
{
	leaf.sum[i,] <- apply(spe.norm[which(spe.ch.mvpart$where==groups.mrt[i]),],
		2, sum)
}
leaf.sum
dev.new(title="Fish composition of 4 leaves")
par(mfrow=c(2,2))
for(i in 1:length(groups.mrt))
{
	pie(which(leaf.sum[i,] > 0), radius=1, main=paste("leaf #", groups.mrt[i]))
}

# Extracting MRT results from an mvpart object
# Packages MVPARTwrap and rdaTest must have been loaded
spe.ch.mvpart.wrap <- MRT(spe.ch.mvpart, percent=10, species=colnames(spe))
summary(spe.ch.mvpart.wrap)

# Indicator species search on the MRT result
spe.ch.MRT.indval <- indval(spe.norm, spe.ch.mvpart$where)
spe.ch.MRT.indval$pval		# Probability

# For each significant species, find the leaf with the highest IndVal
spe.ch.MRT.indval$maxcls[which(spe.ch.MRT.indval$pval <= 0.05)]

# IndVal value in the best leaf for each significant species
spe.ch.MRT.indval$indcls[which(spe.ch.MRT.indval$pval <= 0.05)]



# Clustering with sequential constraint
# *************************************

dev.new(title="MRT with sequential constraint", width=14, height=7)
par(mfrow=c(1,2))
spe.ch.seq <- mvpart(as.matrix(spe.norm) ~ dfs, env, cp=0, xv="pick", margin=0.08,
	xval=nrow(spe.norm), xvmult=100, which=4)
# Here, click on the desired number of groups

summary(spe.ch.seq)

# Group composition (labels of terminal nodes)
(gr <- spe.ch.seq$where)

# Renumber clusters sequentially
aa <- 1
gr2 <- rep(1, length(gr))
for (i in 2:length(gr))
{
	if (gr[i]!=gr[i-1]) aa <- aa+1
	gr2[i] <- aa
}

# Plot the clusters on a map of the Doubs river
dev.new(title="MRT groups on river")
plot(spa, asp=1, type="n", main="MRT groups",
	xlab="x coordinate (km)", ylab="y coordinate (km)")
lines(spa, col="light blue")
text(80, 10, "Upstream", cex=1.2)
text(15, 70, "Downstream", cex=1.2)
k <- length(levels(factor(gr2)))
for (i in 1:k)
{
	points(spa[gr2==i,1], spa[gr2==i,2], pch=i+20, cex=3, col=i+1, bg=i+1)
}
text(spa, row.names(spa), cex=0.8, col="white", font=2)
legend("bottomright", paste("Cluster", 1:k), pch=(1:k)+20, col=2:(k+1),
	pt.bg=2:(k+1), pt.cex=2, bty="n")



# Fuzzy c-means clustering of the fish species data
# *************************************************

k <- 4		# Choose the number of clusters
spe.fuz <- fanny(spe.ch, k=k, memb.exp=1.5)
summary(spe.fuz)

# Site fuzzy membership
spe.fuz$membership
# Nearest crisp clustering
spe.fuz$clustering
spefuz.g <- spe.fuz$clustering

# Silhouette plot
dev.new(title="Fuzzy clustering of fish data - Silhouette plot")
plot(silhouette(spe.fuz), main="Silhouette plot - Fuzzy clustering",
	cex.names=0.8, col=spe.fuz$silinfo$widths+1)

# Ordination of fuzzy clusters (PCoA)
# -----------------------------------
# Step 1: ordination (PCoA) of the fish chord distance matrix
dc.pcoa <- cmdscale(spe.ch)
dc.scores <- scores(dc.pcoa, choices=c(1,2))

# Step 2: ordination plot of fuzzy clustering result
dev.new(title="Fuzzy clustering of fish data - Ordination plot")
plot(scores(dc.pcoa), asp=1, type="n",
	main="Ordination of fuzzy clusters (PCoA)")
abline(h=0, lty="dotted")
abline(v=0, lty="dotted")

# Step 3: representation of fuzzy clusters
for (i in 1:k)
{
	gg <- dc.scores[spefuz.g==i,]
	hpts <- chull(gg)
	hpts <- c(hpts, hpts[1])
	lines(gg[hpts,], col=i+1)
}
stars(spe.fuz$membership, location=scores(dc.pcoa), draw.segments=TRUE,
	add=TRUE, scale=FALSE, len=0.1, col.segments=2:(k+1))
legend(locator(1), paste("Cluster", 1:k, sep=" "),
	pch=15, pt.cex=2, col=2:(k+1), bty="n")
# Click on the graph to position legend

