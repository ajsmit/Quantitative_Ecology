# Quantitative Ecology (BCB743)
# non-Metric Multidimensional Scaling
# Author: AJ Smit
# Date: 20 July 2021

library(betapart)
library(vegan)
library(gridExtra)
library(grid)
library(gridBase)
library(tidyr)

spp <- read.csv('/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/diversity/SeaweedsSpp.csv')
spp <- dplyr::select(spp, -1)
dim(spp)

Y.core <- betapart.core(spp)
Y.pair <- beta.pair(Y.core, index.family = "sor")

# Let Y1 be the turnover component (beta-sim):
Y1 <- as.matrix(Y.pair$beta.sim)

load("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/diversity/SeaweedEnv.RData")
dim(env)

# I select only some of the thermal vars; the rest are collinear with some of the ones I import:
E1 <- dplyr::select(env, febMean, febRange, febSD, augMean,
                    augRange, augSD, annMean, annRange, annSD)

# Calculate z-scores:
E1 <- decostand(E1, method = "standardize")

bioreg <- read.csv('/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/diversity/bioregions.csv', header = TRUE)
head(bioreg)

sites <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/diversity/sites.csv")
sites <- sites[, c(2, 1)]
head(sites)
dim(sites)

# fit the full model:
rda_full <- capscale(Y1 ~., E1)
rda_full
# summary(rda_full)

anova(rda_full, parallel = 4) # ... yes!

rda_full_R2 <- RsquareAdj(rda_full)$adj.r.squared
round(rda_full_R2, 2)

round(sum(rda_full$CCA$eig) / rda_full$tot.chi * 100, 2) # in %

vif.cca(rda_full)

E2 <- dplyr::select(E1, -annMean)
rda_sel1 <- capscale(Y1 ~., E2)
vif.cca(rda_sel1)

E3 <- dplyr::select(E2, -febMean)
rda_sel2 <- capscale(Y1 ~., E3)
vif.cca(rda_sel2)

rda_final <- rda_sel2

# is the fit significant?
anova(rda_final, parallel = 4) # ... yes!

anova(rda_final, by = "axis", parallel = 4) # ... yes!

(rda_final_axis_test <- anova(rda_final, by = "terms", parallel = 4))

rda_final_ax <- which(rda_final_axis_test[, 4] < 0.05)
rda_final_sign_ax <- colnames(E3[,rda_final_ax])
rda_final_sign_ax

round(rda_final_R2 <- RsquareAdj(rda_final)$adj.r.squared, 2) # %

round(sum(rda_final$CCA$eig) / rda_final$tot.chi * 100, 2)

scores(rda_final, display = "bp", choices = c(1:2))

# The ordiplots in Fig. 2:
# use scaling = 1 or scaling = 2 for site and species scaling, respectively
rda_final_scrs <- scores(rda_final, display = c("sp","wa","lc","bp","cn"))
# see ?plot.cca for insight into the use of lc vs wa scores
# below I splot the wa (site) scores
site_scores <- data.frame(rda_final_scrs$constraints)
site_scores$bioreg <- bioreg$bolton
site_scores$section <- seq(1:58)
colnames(site_scores) <- c("x", "y", "Bioregion", "Section")

biplot_scores <- data.frame(rda_final_scrs$biplot)
biplot_scores$labels <- rownames(biplot_scores)
biplot_scores_sign <- biplot_scores[biplot_scores$labels %in% rda_final_sign_ax,]

ggplot(data = site_scores, aes(x, y, colour = Bioregion)) +
  geom_point(size = 5.0, shape = 24, fill = "white") +
  geom_text(aes(label = Section), size = 3.0, col = "black") +
  geom_label(data = biplot_scores_sign,
            aes(CAP1, CAP2, label = rownames(biplot_scores_sign)),
            color = "black") +
  geom_segment(data = biplot_scores_sign,
               aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "lightseagreen", alpha = 1, size = 0.7) +
  xlab("CAP1") + ylab("CAP2") +
  labs(colour = "Bioregions") +
  ggtitle(expression(paste("Significant thermal variables and ", beta[sim]))) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 0.8)


