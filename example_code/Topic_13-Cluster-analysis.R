library(tidyverse)
# library(GGally)
library(cluster)
# library(dendextend)
library(ggcorrplot)
library(factoextra)
# library(gridExtra)
library(vegan)
# library(Rtsne) # for t-SNE plot

SDGs <- read_csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/SDG_complete.csv")
head(SDGs)

unique(SDGs$ParentLocation)
length(unique(SDGs$ParentLocation))
length(SDGs$Location)

colnames(SDGs[, 3:ncol(SDGs)])

corr <- round(cor(SDGs[3:ncol(SDGs)]), 1)
corr[lower.tri(corr)] <- NA
ggcorrplot(corr, type = 'upper', outline.col = "white",
           colors = c("#1679a1", "white", "#f8766d"),
           lab = TRUE)

SDGs_std <- decostand(SDGs[3:ncol(SDGs)], method = "standardize")
rownames(SDGs_std) <- SDGs$Location


# pam clustering ----------------------------------------------------------

# number of clusters?

# using silhouette analysis
fviz_nbclust(SDGs_std, cluster::pam, method = "silhouette") + theme_grey()

# total within sum of square
fviz_nbclust(SDGs_std, cluster::pam, method = "wss") + theme_grey()

# gap statistics
fviz_nbclust(SDGs_std, cluster::pam, method = "gap_stat") + theme_grey()

SDGs_pam <- pam(SDGs_std, metric = "euclidean", k = 3)
SDGs_pam

# Create plots ------------------------------------------------------------

fviz_cluster(SDGs_pam, geom = "point", ellipse.type = "convex") +
  geom_text(aes(label = SDGs$Location), size = 2.5)


# scale SA bigger for plotting

SDGs <- SDGs |>
  mutate(col_vec = ifelse(Location == "South Africa", "black", "grey50"),
         scale_vec = ifelse(Location == "South Africa", 3.5, 2.5))

fviz_cluster(SDGs_pam, geom = "point", ellipse.type = "convex") +
  geom_text(aes(label = SDGs$Location), size = SDGs$scale_vec, col = SDGs$col_vec)


SDGs <- SDGs |>
  mutate(cluster = SDGs_pam$clustering)

SDGs_clustered2 <- SDGs |>
  group_by(cluster) |>
  summarise_at(vars(other_1:SDG3.b_5), mean, na.rm = TRUE)

pairs(SDGs[, 3:15], col = c("#FC4E07", "navy")[SDGs_pam$clustering])

SDGs_pam

