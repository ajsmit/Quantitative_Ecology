# load the libraries required for the analysis:
library(vegan)
library(tidyverse)

# read in the plantations data
plantations <- read.csv("plantations.csv")

# replace 'NA' with zero, because NAs don't work with our analysis
plantations[is.na(plantations)] <- 0

# ... base rest of script on the 'BDC334_Lab2_Doubs_data.R' file...
