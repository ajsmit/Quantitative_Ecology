# Quantitative Ecology (BCB743)
# Environmental Distance
# Author: AJ Smit
# Date: 29 June 2021


# Load the packages -------------------------------------------------------

library(vegan)
library(ggplot2)
library(geodist) # for calculating geographical distances between lats/lons
library(ggpubr) # to arrange the multipanel graphs


# LOAD THE DATA -----------------------------------------------------------

xyz <- read.csv("exercises/diversity/Euclidian_distance_demo_data_xyz.csv")


# Understand the data -----------------------------------------------------

# the dimensions
dim(xyz)

# what's inside?
xyz


# Calculate Euclidian distances -------------------------------------------

xyz_euc <- round(vegdist(xyz[, 2:4], method = "euclidian", upper = FALSE, diag = TRUE), 4) # select only cols 2, 3 and 4
xyz_euc
class(xyz_euc)

# convert to dataframe
xyz_df <- as.data.frame(as.matrix(xyz_euc))
class(xyz)
xyz_df

# now imagine these are actual environmental data
env_fict <- read.csv("exercises/diversity/Euclidian_distance_demo_data_env.csv")
env_fict

env_fict_euc <- round(vegdist(env_fict[, 2:4], method = "euclidian", upper = FALSE, diag = TRUE), 4) # select only cols 2, 3 and 4
env_fict_euc


# A LOOK AT THE SEAWEED ENVIRONMENTAL DATA --------------------------------

load("exercises/diversity/SeaweedEnv.RData")


# Understand the data -----------------------------------------------------

# lets look at the data's dimensions
dim(env)

# the top five lines and first five columns
round(env[1:5, 1:5], 4)

# the bottom five lines and last five columns
round(env[(nrow(env) - 5):nrow(env), (ncol(env) - 5):ncol(env)], 4)

# the names of the variables
colnames(env)

# select only some of the thermal variables
env1 <- dplyr::select(env, febMean, febRange, febSD, augMean,
                      augRange, augSD, annMean, annRange, annSD)
dim(env1)


# Calculate z-scores ------------------------------------------------------

E1 <- round(decostand(env1, method = "standardize"), 4)
E1[1:5, 1:5]


# Calculate Euclidian distance --------------------------------------------

E1_euc <- round(vegdist(E1, method = "euclidian", upper = TRUE), 4)
E1_df <- as.data.frame(as.matrix(E1_euc))
E1_df[1:10, 1:10]


# Make a plot -------------------------------------------------------------

ggplot(data = E1_df, (aes(x = 1:58, y = `1`))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("Environmental distance")


# EUCLIDIAN DISTANCE OF GEOGRAPHICAL COORDINATES --------------------------

geo <- read.csv("exercises/diversity/sites.csv")
dim(geo)
head(geo)


# Calculate geographic distances between sections -------------------------

dists <- geodist(geo, paired = TRUE, measure = "geodesic")
dists_df <- as.data.frame(as.matrix(dists))
colnames(dists_df) <- seq(1:58)
dists_df[1:5, 1:5]


# Make a plot -------------------------------------------------------------

plt1 <- ggplot(data = dists_df, (aes(x = 1:58, y = `1`/1000))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("Distance (km)") + ggtitle("Actual geographic distance")


# Calculate Euclidian distances btw sections ------------------------------

dists_euc <- vegdist(geo, method = "euclidian")
dists_euc_df <- round(as.data.frame(as.matrix(dists_euc)), 4)
dists_euc_df[1:5, 1:5]


# Make a plot -------------------------------------------------------------

plt2 <- ggplot(data = dists_euc_df, (aes(x = 1:58, y = `1`))) +
  geom_line() + xlab("Coastal section, west to east") + ylab("Euclidian distance") + ggtitle("Euclidian distance")


# Arrange plots -----------------------------------------------------------

ggarrange(plt1, plt2, nrow = 2)



