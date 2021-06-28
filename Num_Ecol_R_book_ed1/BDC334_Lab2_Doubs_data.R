################################################################################
### CHAPTER 2: EXPLORATORY DATA ANALYSIS
###
### Online supporting material for:
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 30 January 2017
################################################################################

# Load required packages
library(vegan)
library(tidyverse)

# Import the data from CSV files
# (files must be in the working directory)
# Species (community) data frame (fish abundances)
spe <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsSpe.csv", row.names = 1))
spe <- spe[-8,] # this row has all zeroes; remove it
spe

# calculate the species dissimilarities based on the Jaccard Index
spe_dist <- round(vegdist(spe, method = "jaccard", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist))

# Environmental data frame
env <- as.tibble(read.csv("Num_Ecol_R_book_ed1/DoubsEnv.csv", row.names = 1))
env <- env[-8,] # this row has all zeroes in the species table; remove it, but why?
env

# this will calculate the dissimilarities based on ALL of the environmental
# data in this table
# but first we standardize; why?
(env_std <- as.tibble(round(decostand(env, method = "standardize"), 2)))
env_dist_all <- round(vegdist(env_std, method = "euclidian", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist_all))

# we can also calculate distances based on only one (or a few) variables, e.g.
env_dist_alt <- round(vegdist(env_std[, "alt"], method = "euclidian", diag = TRUE, upper = TRUE), 2)
as.tibble(as.matrix(env_dist_alt))

# the matrices can be exported to CSV files for plotting in Excel; for added
# levels of amazingness you'll of course want to create the plots in R;
# although a nice thought, I won't expect that you soar to the lofty
# heights where my high levels of expectations reside... :-(

# At this point AJ will recap the matrices
