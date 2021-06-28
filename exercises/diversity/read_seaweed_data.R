library(vegan)

spp <- read.csv('exercises/diversity/SeaweedsSpp.csv')
spp <- dplyr::select(spp, -1)
sor <- vegdist(spp, binary = TRUE, diag = TRUE) # this is jaccard
sor.df <- as.data.frame(as.matrix(sor))
write.csv(sor.df, file = "exercises/diversity/SeaweedSpp_dis_matrix.csv")

env <- read.csv('exercises/diversity/SeaweedsEnv.csv', header = TRUE)
# Calculate z-scores:
E1 <- decostand(env, method = "standardize")
E1.euc <- vegdist(E1, method = "euclidian", upper = TRUE)
E1.df <- as.data.frame(as.matrix(E1.euc))
write.csv(E1.df, file = "exercises/diversity/SeaweedEnv_dist_matrix.csv")

# I select only some of the thermal vars; the rest
# are collinear with some of the ones I import:
E1 <- dplyr::select(env, febMean, febRange, febSD, augMean,
                    augRange, augSD, annMean, annRange, annSD)

library(geodist)
geo_dist <- read.csv("exercises/diversity/sites.csv")
# calc pairwise distances in m
dists <- geodist(geo_dist, paired = TRUE, measure = "geodesic")
dists_df <- as.data.frame(as.matrix(dists))
colnames(dists_df) <- seq(1:58)
write.csv(dists_df, file = "exercises/diversity/Seaweed_geodist.csv")

dists/1000

# do a PCoA on environmental data
env.pcoa <- cmdscale(E1.euc, k = (nrow(E1) - 1), eig = TRUE)
env.pcoa

# importance of first eigenvalue?
env.pcoa$eig[1]/sum(env.pcoa$eig) * 100

# do a PCoA on seaweed data
sor.pcoa <- cmdscale(sor, eig = TRUE)

sor.pcoa$eig[1]/sum(sor.pcoa$eig) * 100


ordiplot(scores(env.pcoa, choices = c(1,2), type = "t"))
abline(v = 0, h = 0, lty = 3)
env.wa <- wascores(env.pcoa$points[,1:2], env)
text(env.wa, rownames(env.wa), cex = 0.7, col = "red")

# do an nMDS
env.nmds <- metaMDS(env[,c(2:7)], distance = "euclidian")
env.nmds
plot(env.nmds, type = "t")
env.nmds$species
par(mfrow = c(1,2))
stressplot(env.nmds, main = "Shepard plot")
gof <- goodness(env.nmds)
plot(env.nmds, type = "t", main = "Goodness of fit")
points(env.nmds, display = "sites", cex = gof * 300)

