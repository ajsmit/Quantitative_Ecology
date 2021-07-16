# Quantitative Ecology (BCB743)
# Principal Coordinate Analysis - mtcars example using prcomp() and ggbiplot()
# Author: AJ Smit
# Date: 16 July 2021

# install.packages("devtools")
# library(devtools)
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# install_github("vqv/ggbiplot")

library(ggbiplot)
library(vegan)

data(mtcars)
head(mtcars)

# using the built-in prcomp() function, not vegan's rda()
cars_pca <- prcomp(mtcars[,c(1:7,10,11)], center = TRUE,scale. = TRUE)
summary(cars_pca)

# using ggbiplot to create ordination plots
ggbiplot(cars_pca)
ggbiplot(cars_pca, labels = rownames(mtcars))

# assign the cars to the countries within which they are used
mtcars_country <- c(rep("Japan", 3), rep("US",4), rep("Europe", 7),
                    rep("US",3), "Europe", rep("Japan", 3), rep("US",4),
                    rep("Europe", 3), "US", rep("Europe", 3))

# more plots...
ggbiplot(cars_pca, ellipse = TRUE, labels = rownames(mtcars),
         groups = mtcars_country)

ggbiplot(cars_pca, ellipse = TRUE, choices = c(3,4), labels = rownames(mtcars),
         groups = mtcars_country)

ggbiplot(cars_pca, ellipse = TRUE, obs.scale = 1, var.scale = 1, labels = rownames(mtcars),
         groups = mtcars_country)

ggbiplot(cars_pca, ellipse = TRUE, obs.scale = 1, var.scale = 1, labels = rownames(mtcars),
         groups = mtcars_country) +
  scale_colour_manual(name = "Origin", values= c("forest green", "red3", "dark blue")) +
  ggtitle("PCA of mtcars dataset") +
  theme_minimal() +
  theme(legend.position = "bottom")
