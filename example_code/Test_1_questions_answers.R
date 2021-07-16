# BCB743 Quantitative Ecology practice test 				Total 80 marks

# Included, please find three files:
# 1.	‘seaweeds.csv’ – presence/absence of macroalgal species along a coastline
# 2.	‘env.csv’ – environmental conditions at the macroalgal species collection sites
# 3.	‘sites.csv’ – the geographic coordinates of the sites
#
# The environmental conditions include the mean and minimum August temperatures, the mean and maximum February temperatures, and the mean August and February water column chlorophyll-a biomass in the water column at the sites.

library(tidyverse)
library(ggpubr)
library(vegan)


# Question 1 (2 marks) ----------------------------------------------------
# Specify the dimensions of the three data sets.

env <- as.tibble(read.csv("data/seaweed/env.csv")[, -1])
seaweeds <- as.tibble(read.csv("data/seaweed/seaweeds.csv")[, -1])
sites <- as.tibble(read.csv("data/seaweed/sites.csv"))

# dimensions given as rows x columns
dim(env)
dim(seaweeds)
dim(sites)


# Question 2 (12 marks) ---------------------------------------------------
# Provide a i) table and ii) figure(s) of the descriptive statistics of the environmental variables.

# i)
desc_stats <- function(x) {
  tibble(min = min(x),
         q2 = quantile(x, 0.25),
         med = quantile(x, 0.5),
         mean = mean(x),
         q4 = quantile(x, 0.75),
         max = max(x))
}

plyr::ldply(apply(env, 2, desc_stats))

# or simply

summary(env)

# ii)
# the temperature vars
plt1 <- env %>%
  select(-augChl, -febChl) %>%
  gather(key = var, value = measurement) %>%
  ggplot(aes(x = var, y = measurement)) +
    geom_boxplot() +
    labs(x = "Temperature variable",
         y = "Temperature (°C)",
         title = "Thermal variables") +
    theme_bw()

plt2 <- env %>%
  select(augChl, febChl) %>%
  gather(key = var, value = measurement) %>%
  ggplot(aes(x = var, y = measurement)) +
    geom_boxplot() +
    labs(x = "Chlorophyll-a variable",
         y = "Chlorophyll-a (mg/m3)",
         title = "Chlorophyll-a variables") +
    theme_bw()

ggarrange(plt1, plt2, ncol = 2)

# or
boxplot(env[,1:6])


# Question 3 (6 marks) ----------------------------------------------------
# Provide a profile plot (variable vs. site number) of the variables recorded for the month of February.

plt3 <- env %>%
  select(febMean, febMax) %>%
  gather(key = var, value = measurement) %>%
  group_by(var) %>%
  mutate(site = seq(1:nrow(env))) %>%
  ungroup() %>%
  ggplot(aes(x = site, y = measurement)) +
    geom_line(aes(col = var, group = var)) +
    labs(x = "Site along coast",
         y = "Temperature (°C)",
         title = "Thermal variables") +
    theme_bw() +
    theme(legend.position = c(0.8, 0.2),
          legend.title = element_blank())

plt4 <- env %>%
  select(febChl) %>%
  gather(key = var, value = measurement) %>%
  group_by(var) %>%
  mutate(site = seq(1:nrow(env))) %>%
  ungroup() %>%
  ggplot(aes(x = site, y = measurement)) +
  geom_line(aes(col = var, group = var), show.legend = FALSE) +
  labs(x = "Site along coast",
       y = "Chlorophyll-a (mg/m3)",
       title = "Chlorophyll-a variable") +
  theme_bw()

ggarrange(plt3, plt4, nrow = 2)


# Question 4 (5 marks) ----------------------------------------------------
# i.	Calculate an association matrix for the species data. Show the code used to produce the calculation. (1)
# ii.	Briefly describe the meaning of the data represented in the association matrix, make use of specific references to cells within the matrix in your explanation. (4)

# i)
seaweeds.t <- t(seaweeds)
seaweeds.t[1:20, 1:20]
seaweeds.t.S7 <- vegdist(seaweeds.t, binary = TRUE)
round(as.matrix(seaweeds.t.S7)[1:15, 1:15], 2)

# ii)
# I'll explain again in words in the class. Can't be bothered typing it out now...

# Question 5 (20) ---------------------------------------------------------
# i.	Undertake a PCA on the environmental data, and provide a written description of each line of code used (in addition to showing the code). (2)
# ii.	How is the total inertia calculated? (2)
# iii.	What do the Eigenvalues represent? (2)
# iv.	What do the Eigenvectors represent? (2)
# v.	Provide the R code for how to derive the “Proportion Explained”. (2)
# vi.	What proportion of the variation is explained by the first three PC axes individually? (3)
# vii.	What proportion of the variation is cumulatively explained by the first three axes? (1)
# viii.	How would you interpret the “Species Scores”? (3)
# ix.	How would you interpret the “Site Scores”? (3)

# i)
# use the environmental data and apply a principal components analysis; the `scale = TRUE` argument standardises the data prior to the calculation in the same way that `decostand(env, method = "standardize")` would; this creates a correlation matrix, which is necessary because the variables each have a different measurement scale
cor(env)
env.pca <- rda(env, scale = TRUE)

# ii)
# the total iniertia is the sum of the diagonal of the correlation matrix that feeds into the PCA

# iii)
# The eigenvalues indicate the importance of the PC axes, with the first PC axis always having the highest eigenvalue, diminishing until the smallest one is found for the last PC axis; the sum of all the eigenvalues equals the total inertia

# iv)
# The eigenvectors represent the 'loadings' of the original, untransformed environmental variables on the corresponding eigenvalue (PC axis); they indicate the importance of the original variable along the new reduced dimension captured by the PC axes

# v)
# For the first PC axis, the proportion of variation explained is given by
(p1 <- env.pca$CA$eig[1] / sum(env.pca$CA$eig))

# vi)
# For the first, see v); for the rest, the code is
(p2 <- env.pca$CA$eig[2] / sum(env.pca$CA$eig))
(p3 <- env.pca$CA$eig[3] / sum(env.pca$CA$eig))

# vii)
# The cumulative inertia (%) of the first three axes is
((p1 + p2 + p3) * 100)

# viii)
# In this case, it being an environmental data set, the species scores refer to the position of the arrow heads on the PCA biplot; they indicate the direction and relative magnitude of the environmental  variable it is associated with as a linear gradient across the 'landscape', i.e. the reduced space represented by the corresponding PC axes

# ix)
# The coordinates of the sites (rows) on the ordination (Euclidian) space, i.e.
# the ordination diagram


# Question 6 (10) ---------------------------------------------------------
# i.	Provide biplots using both scaling 1 and scaling 2 of the environmental data PCA, and show the necessary code. (4)
# ii.	Interpret the figures (6)

# i)
source("Num_Ecol_R_book_ed1/cleanplot.pca.R")
par(mfrow = c(2, 2))
biplot(env.pca, scaling = 1, main = "PCA - scaling 1 (sites)")
biplot(env.pca, main = "PCA - scaling 2 ('species')") # default scaling 2
cleanplot.pca(env.pca, scaling = 1, mar.percent = 0.08) # sites scaling
cleanplot.pca(env.pca, scaling = 2, mar.percent = 0.04) # species scaling

# ii)
# In the left figure we use sites scaling (scaling 1) so the sites are scaled by eigenvalues. Here the distances in multidimensional space are more accurately reflected on the graph plane, so it shows the relationships between sites better. Sites that share similar environmental characteristcs plots closer together and the ones that are further appart on the graph are also more dissimilar in their actual multidimensional space. The circle of equilibrium on the bottom left panel shows that the February mean temperature and the two chlorophyll-a variables can be interpreted as certainly affecting the spread of sites. The higher chl-a concentrations are generally associated with the spread of sites towards the left of the vertical zero line, although their influence is felt rather obliquely (not perpectly alligned with the direction of spread in sites); this might indicate the another important variable is lacking that can account for this direction of spread, but certainly the chl-a can be used to explain some of the spread of sites in this direction. Towards the right of the graph (again talking about the bottom left figure), the sites here are primarily influenced by the thermal variables. There sites have higher temperatures during both the coldest and warmest months of the year. Overall, the sites are alligned along a gradient (low site numbers on the left, high number to the right). Plot a map of the sites and see where they are, and you will understand why the gradient exists... The species scaling plots on the left show how the influential variables are tightly correlated with each other. On the left of the bottom right-hand graph, the chl-a values are strongly correlated, and on the right the thermal ones are. This could (and probably does) indicate some collinearity and it shoulld be examined and thought about.


# Question 7 (25) ---------------------------------------------------------
# Analyse the species data using the multivariate technique of your choice.
# i.	Justify the use of your chosen multivariate method, and describe the reason for all the arguments specified in the function call, if necessary. (5)
# ii.	Provide biplots for the species ordination, and show the code that you used to produce your figures. Show figures for both scalings. (4)
# iii.	Explain your findings. (6)
# iv.	In the light of the environmental PCA performed earlier, give possible explanations for the patterns that you find. Demonstrate that you are able to use the combined understanding from all aforementioned analyses to influence your arguments. (10)

# i) and ii)
# A dissimilarity matrix can be accepted by both a NMDS and a PCoA. I choose the Jaccard dissimilarity because it is suitable and widely used for presence/absence data.

# Both PCoA and NMDS chosen because they can accept dissimilarity matrices. NMDS is, however, very well suitable for ordination as this rank-based approach is more 'flexible' regarding how the distances in the observed multivariate space are mapped to the ordination distance. The relationship is non-linear, as can be seen in the Shepard diagram. NMDS more faithfully retains the pairwise dissimilarity between objects in a low-dimensional space. Here I start with a Jaccard dissmilarity matrix, which will be supplied to both the NMDS and the PCoA.
seaweeds.jac <- vegdist(seaweeds, method = "jaccard", binary = TRUE)

# a NMDS...
seaweeds.nmds <- metaMDS(seaweeds.jac)

# ...or a PCoA
seaweeds.pcoa <- cmdscale(seaweeds.jac, k = (nrow(seaweeds) - 1), eig = TRUE)

par(mfrow = c(2, 2))
stressplot(seaweeds.nmds, seaweeds.jac)
ordiplot(seaweeds.nmds, type = "t", display = c("sites"),
         main = "NMDS with site scores")
ordiplot(scores(seaweeds.pcoa, choices = c(1, 2)), type = "t",
         main = "PCoA with site scores")
abline(v = 0, h = 0, lty = 3)

# iii)
# The NMDS plot has a low stress, which means that the pairwise dissimilarities are well represented in ordination space (there is little scatter around the red line in the Shepard diagram). We can interpret the graph with confidence. The shape of the arrangement of points in the both the NMDS and the PCoA is similar, but the 'horsehoe effect' is more visible in the PCoA. The curvature in the string of sites reflect real environmental gradients that act strongly from left to right, and somewhat less strongly top to bottom. If you had done a plot of the sites on a map you can see how this arrangement of sites is linked to the actual environment: the sites here, which were plotted in the complete absence of geographic information, nevertheless are arranged according to environmental constraints that are tighly geographically constrained!

# iv)
# For a comprehensive assessment of the situation, see the full findings and explanation in the paper Smit, A. J., Bolton, J. J., & Anderson, R. J. (2017). Seaweeds in Two Oceans: Beta-Diversity. Frontiers in Marine Science, 4(December), 1–15. http://doi.org/10.3389/fmars.2017.00404
