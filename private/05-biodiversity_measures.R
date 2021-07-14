# Quantitative Ecology (BCB743)
# Measures of Biodiversity
# Author: AJ Smit
# Date: 29 June 2021


# Load the packages -------------------------------------------------------

library(vegan)
library(ggplot2)


# LOAD THE SEAWEED DATA ---------------------------------------------------

spp <- read.csv('exercises/diversity/SeaweedsSpp.csv')


# Understand the data -----------------------------------------------------

# how many rows and columns?
dim(spp)

# view the first 10 rows and columns
spp[1:10, 1:10]

# remove the first column
spp <- dplyr::select(spp, -1)


# Look at species richness (alpha-diversity) ------------------------------

library(BiodiversityR)
spp_richness <- diversityresult(spp, index = 'richness', method = 'each site')
# read the help file: it can also be used to calculate Shannon, Simpson, etc.
# setting method = pooled calculate gamma diversity

# or of you cannot install/run BiodiversityR, use vegan...

specnumber(spp, MARGIN = 1)

ggplot(data = spp_richness, (aes(x = 1:58, y = richness))) +
  geom_line(col = "red4") + xlab("Coastal section, west to east") + ylab("Species richness")
# we have presence/absence data, so we


# Univariate diversity indices --------------------------------------------

# only wirks with abundance data, so
# load a different dataset
light <- read.csv("exercises/diversity/light_levels.csv")
light

# calculate the jaccard dissimilarity
light_div <- data.frame(
  site = c("low_light", "mid_light", "high_light"),
  richness = specnumber(light[, 2:7], MARGIN = 1),
  shannon = round(diversity(light[, 2:7], MARGIN = 1, index = "shannon"), 2),
  simpson = round(diversity(light[, 2:7], MARGIN = 1, index = "simpson"), 2)
)
light_div


# Dissmilarity indices ----------------------------------------------------

# back to seaweed data...
# presence absence data, therefore
# calculate the sørensen dissmilarity
# by setting binary = TRUE
sor <- round(vegdist(spp, binary = TRUE, diag = TRUE), 4)

# what are the class and dimensions?
class(sor)
str(sor)

# create a dataframe for easy use and check it out
sor_df <- as.data.frame(as.matrix(sor))
class(sor_df)
dim(sor_df)
sor_df[1:10, 1:10]

# can write to csv if necessary
# write.csv(sor.df, file = "exercises/diversity/SeaweedSpp_dis_matrix.csv")


# Gamma diversity ---------------------------------------------------------

# the number of columns gives the total
# number of species in this example
ncol(spp)

# we may also use
diversityresult(spp, index = 'richness', method = 'pooled')
# why this difference?


# Beta-diversity ----------------------------------------------------------

# Whittaker's concept of β-diversity
# true beta
true_beta <- data.frame(
  beta = specnumber(spp, MARGIN = 1) / ncol(spp),
  section_no = c(1:58)
)

ggplot(data = true_beta, (aes(x = section_no, y = beta))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("True beta-diversity")

# absolute species turnover
abs_beta <- data.frame(
  beta = ncol(spp) - specnumber(spp, MARGIN = 1),
  section_no = c(1:58)
)

ggplot(data = abs_beta, (aes(x = section_no, y = beta))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("Absolute beta-diversity")


# Contemporary definitions β-diversity ------------------------------------

library(betapart)

# decompose total Sørensen dissimilarity into
# turnover and nestedness-resultant components

# compute the basic quantities needed for computing
# the multiple-site beta diversity measures and pairwise
# dissimilarity matrices
Y.core <- betapart.core(spp)
str(Y.core)

# compute 3 distance matrices accounting for the
# (i) turnover (replacement),
# (ii) nestedness-resultant component, and
# (iii) total dissimilarity (i.e. the sum of both components)
Y.pair <- beta.pair(Y.core, index.family = "sor")
str(Y.pair)

# Let Y1 be the turnover component (beta-sim):
Y1 <- as.matrix(Y.pair$beta.sim)

# Let Y2 be the nestedness-resultant component (beta-sne):
Y2 <- as.matrix(Y.pair$beta.sne)
