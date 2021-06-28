# multivariate_analysis.R

# This script reproduces the analysis in Mayombo et al (2019).

library(vegan)
#library(BiodiversityR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(tibble)

# https://www.rdocumentation.org/packages/BiodiversityR/versions/2.10-1/topics/nested.anova.dbrda

# Reviewer 1
# The design of the observational study includes 2 treatments - age (young versus old) and host species (Laminaria versus Ecklonia), 4 replicates (4 primary blades from each combination of host algae and age), and 3 subsamples from each blade (pseudoreplicates, if treated incorrectly as replicates). The experimental design is analogous to a 2-way ANOVA, but with community data instead of a single individual response variable. This design can evaluate interactive effects between the two treatments (age and species). The authors’ experimental design is most suited to analyses using PERMANOVA, which is the community statistics version of the ANOVA.

# Please indicate for the readers why the data were transformed and standardised using the stated procedures. Definitely a good idea to transform data, but the readers need to understand why particular procedures were employed. Please describe the Wisconsin double standardisation (row/column standardised by row/column total – to produce relative abundance to total and column/row standardised by column/row max – to produce abundance relative to species max abundance). Why a double standardisation + square-root transformation, as opposed to a single row/column standardization by row/column total + square-root transformation?

# Please indicate for the readers why the data were transformed and standardised using the stated procedures. Definitely a good idea to transform data, but the readers need to understand why particular procedures were employed. Please describe the Wisconsin double standardisation:
# * row/column standardised by row/column total–to produce relative abundance to total and column/row standardised; vs.
# * column/row max–to produce abundance relative to species max abundance.
# Why a double standardisation + square-root transformation, as opposed to a single row/column standardisation by row/column total + square-root transformation?

# AJS: About ANOSIM and PERMANOVA
# "Overall, ANOSIM and the Mantel test were very sensitive to heterogeneity in dispersions, with ANOSIM generally being more sensitive than the Mantel test. In contrast, PERMANOVA and Pillai’s trace were largely unaffected by heterogeneity for balanced designs. [...]. PERMANOVA was also unaffected by differences in correlation structure. [...] PERMANOVA was generally, but not always, more powerful than the others to detect changes in community structure"
# Anderson, M. J., & Walsh, D. C. I. (2013). PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological Monographs, 83(4), 557–574. http://doi.org/10.1890/12-2010.1

# AJS: About data transformation... Useful when the range of data values is very large. Data are square root transformed, and then submitted to Wisconsin double standardization, or species divided by their maxima, and stands standardized to equal totals. These two standardizations often improve the quality of ordinations, but we forgot to think about them in the initial analysis.


# Get the data in ---------------------------------------------------------

# The diatom species data include the following:
# columns: diatom genera
# rows: samples (samples taken from two species of kelp; equivalent to sites in other species x sites tables)
# row names correspond to combinations of the factors in the columns inside 'PB_diat_env.csv'
# where "host_size" is 'A' for adult kelp plant (host), 'J' for juvenile kelp plant (host), "host_spp" is 'Lp' for
# kelp species Laminaria pallida (host), 'Em' for kelp plant Ecklonia maxima (host), "plant" is the unique number
# identifying a specific kelp plant, and "rep" is the replicate tissue sample from each kelp host plant from which
# the diatoms were extracted.

# with shortened name to fix nMDS overplotting
spp <- read.csv(file = "exercises/diatoms/PB_data_matrix_abrev.csv",
                row.names = "Replicate", sep = ",", header = TRUE)
# with full names
spp2 <- read.csv(file = "exercises/diatoms/PB_data_matrix.csv",
                 row.names = "Replicate", sep = ",", header = TRUE)
# remove ".spp" from column header name
colnames(spp) <- str_replace(colnames(spp), "\\.spp", "")
colnames(spp2) <- str_replace(colnames(spp2), "\\.spp", "")
# Logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of the logarithm; zeros are left as zeros. Higher bases give less weight to quantities and more to presences.
spp.log <- decostand(spp, method = "log")
spp.log.dis <- vegdist(spp.log, method = "bray")

# load the 'environmental' data:
# the content is described above; these variables are categorical vars -- they are not actually
# 'environmental' data, but their purpose in the analysis is analogous to true environmental data;
# it's simply data that describe where the samples were taken from.
env <- as.tibble(read.csv(file = "exercises/diatoms/PB_diat_env.csv",
                          sep = ",", header = TRUE))
env$plant <- as.factor(env$plant)
env$rep <- as.factor(env$rep)

# With the environmental data (factors), the following analyses can be done:
# Discriminant Analysis (DA)
# ✘ Analysis of Similarities (ANOSIM)
# ✔︎ Permutational Analysis of Variance (PERMANOVA)
# Mantel test

# Below we will do an nMDS and PERMANOVA.


# betadisper --------------------------------------------------------------

# Before doing the PERMANOVA, first check to see if the dispersion is the same
# Homogeneity of groups
# betadisper studies the differences in group homogeneities
# analogous to Levene’s test of the equality of variances;
# the null hypothesis that the population variances are equal;
# can only use one factor as an independent variable

# we test the H0 that the dispersion (variance) in diatom community structure does not differ between
# the two host species (the two species of kelp):
(mod.spp <- with(env, betadisper(spp.log.dis, host_spp)))
anova(mod.spp) # there is no difference in dispersion between the diatom communities on the two host species

# make some plots to see...
par(mfrow = c(2, 2))
plot(mod.spp, sub = NULL)
boxplot(mod.spp)

# the PERMANOVA tests the H0 that there is no difference in the diatom community structure taken from
# Ecklonia maxima and Laminaria pallida
permutest(mod.spp) # there is in fact no difference

# now apply the same procedure to see if host size has an effect:
(mod.size <- with(env, betadisper(spp.log.dis, host_size)))
anova(mod.size)
plot(mod.size)
boxplot(mod.size)
permutest(mod.size) # nope...


# PERMANOVA ---------------------------------------------------------------

# Permutational multivariate analysis of variance using distance matrices
# (Bray-Curtis similarities by default). ANOSIM uses only ranks of Bray-Curtis,
# so the former preserves more information. PERMANOVA allows for variation
# partitioning and allows for more complex designs (multiple factors, nested
# factors, interactions, covariates, etc.)
# adonis2 studies the differences in the group means
# analogous to multivariate analysis of variance

# note that nestedness should be stated in the model formulation
# as well as in the strata
# "If you have a nested error structure, so that you do not want your data be
# shuffled over classes (strata), you should define strata in your permutation"
# -- Jari Oksannen
(perm.1 <- adonis2(spp.log.dis ~ (host_spp * host_size) / plant,
                   strata = plant,
                   method = p, data = env))

# within subjects effects
# (perm.2 <- adonis2(spp.log.dis ~ host_spp * host_size + plant,
#                    strata = plant,
#                    method = p, data = env))


# nMDS --------------------------------------------------------------------

# START MAKING FIGURE -----------------------------------------------------

spp.nmds <- metaMDS(spp.log, k = 2,trymax = 100,
                    distance = "bray", wascores = TRUE)

scores(spp.nmds, display = "species")
scores(spp.nmds, display = "sites")

(ef <- envfit(spp.nmds, env, permu = 999))
plot(spp.nmds, display = "sites", tck = .02, mgp = c(1.8, 0.5, 0))
plot(ef, p.max = 0.1)

# set things up for the panel of plots (Figure 3)
# output will appear in a PDF in the path specified below (change on own computer)
pdf(file = "exercises/diatoms/Figure__3.pdf", width = 8, height = 7)
col <- c("indianred", "darkturquoise")
pch <- c(17, 19)
opar <- par()
plt1 <- layout(rbind(c(1, 1, 2, 2, 3, 3),
                     c(4, 4, 4, 5, 5, 5)),
               heights = c(2, 3),
               respect = TRUE)
layout.show(plt1)
par(mar = c(3,3,1,1))

# plot 1
plot(mod.spp, main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0), col = col, pch = pch,
     sub = NULL)
# plot 2
plot(mod.size, main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0), col = col, pch = pch,
     sub = NULL)

# plot 3
stressplot(spp.nmds,
           tck = .05, mgp = c(1.8, 0.5, 0))

# plot 4
par(mar = c(3,3,2,1))
plot(spp.nmds, display = "sites", type = "n",
     main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0),
     xlim = c(-2, 2), ylim = c(-1, 2))
with(env,
     points(spp.nmds, display = "sites", col = col[host_spp],
            pch = pch[host_spp]))
with(env,
     ordispider(spp.nmds, groups = host_spp,
                label = TRUE,
                col = col))
with(env, ordiellipse(spp.nmds, groups = host_spp,
                      col = col, label = FALSE))
points(spp.nmds, display = "species", pch = 1, col = "deeppink")
orditorp(spp.nmds, display = "species", cex = 0.8,
         col = "black", air = 0.01)

# plot 5
par(mar = c(3, 3, 2, 1))
plot(spp.nmds, display = "sites", type = "n",
     main = NULL,
     tck = .05, mgp = c(1.8, 0.5, 0),
     xlim = c(-2, 2), ylim = c(-1, 2))
with(env,
     points(spp.nmds, display = "sites", col = col[host_size],
            pch = pch[host_size]))
with(env,
     ordispider(spp.nmds, groups = host_size,
                label = TRUE,
                col = col))
with(env, ordiellipse(spp.nmds, groups = host_size,
                      col = col, label = FALSE))
points(spp.nmds, display = "species", pch = 1, col = "deeppink")
orditorp(spp.nmds, display = "species", cex = 0.8,
         col = "black", air = 0.01)
dev.off()
par(opar)

# END MAKING FIGURE -------------------------------------------------------


# mvabund -----------------------------------------------------------------

library(mvabund)
diat_spp <- mvabund(spp2)

# look at the spread of the data using the boxplot function
# not used in paper
par(mar = c(2, 10, 2, 2)) # adjusts the margins
boxplot(spp, horizontal = TRUE, las = 2, main = "Abundance")
# dev.off()
# par(opar)

# check the mean-variance relationship
# it shows that spp with a high mean also have a high var
meanvar.plot(diat_spp)

# Are there differences in the species composition of the diatom spp. sampled?
# Do some of them specialise on particular spp of kelp, while others are more
# generalised?
# Do some occur more on juveniles, while some are on adults, and which ones
# indiscriminately live across age classes?
# Which species?

# plot(diat_spp ~ env$host_spp, cex.axis = 0.8, cex = 0.8, n.vars = 18)
# plot(diat_spp ~ env$host_size, cex.axis = 0.8, cex = 0.8,
#      type = "bx", n.vars = 18)
# replicated in ggplot2 below...

# scale manually for ggplot2 custom plot
log_fun <- function(x) {
  min_x <- min(x[x != 0], na.rm = TRUE)
  a <- log(x) / min_x
  a[which(!is.finite(a))] <- 0
  return(a)
}

plt1 <- spp2 %>%
  mutate(host_size = env$host_size) %>%
  gather(key = species, value = abund, -host_size) %>%
  as_tibble() %>%
  group_by(species) %>%
  mutate(log.abund = log_fun(abund)) %>%
  ungroup() %>%
  ggplot(aes(x = fct_reorder(species, abund, .fun = mean), y = log.abund)) +
  geom_boxplot(aes(colour = host_size), size = 0.4, outlier.size = 0,
               fill = "grey90") +
  geom_point(aes(colour = host_size, shape = host_size),
             position = position_dodge2(width = 0.8),
             alpha = 0.6, size = 2.5) +
  scale_colour_manual(name = "Age", values = c("red3", "blue3")) +
  scale_shape_manual(name = "Age", values = c(17, 19)) +
  annotate("text", x = 15, y = 3, size = 4.5,
           label = expression(paste(italic("p"), "=0.017"))) +
  annotate("text", x = 14, y = 3, size = 4.5,
           label = expression(paste(italic("p"), "=0.004"))) +
  scale_y_continuous(name = "Log abundance") +
  coord_flip() + theme_bw() +
  theme(panel.grid.major = element_line(linetype = "dashed",
                                        colour = "turquoise", size = 0.2),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, color = "black",
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        axis.text.y = element_text(size = 13, color = "black", face = "italic",
                                   margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")),
        axis.title.x = element_text(size = 14, vjust = 5.75, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.ticks = element_line(color = "black", size = 0.5))
ggsave(filename = "exercises/diatoms/Figure__4.pdf", plot = plt1,
       width = 6, height = 10, scale = 0.9)

# try with poisson distribution, but residuals no good
# size_mod1 <- manyglm(diat_spp ~ env$host_spp * env$host_size * env$plant, family = "poisson")
# plot(size_mod1)

# use negative binomial
size_mod2 <- manyglm(diat_spp ~ (env$host_spp * env$host_size) / env$plant,
                     family = "negative binomial")
plot(size_mod2) # better residuals...
anova(size_mod2, test = "wald")
out <- anova(size_mod2, p.uni = "adjusted", test = "wald")
out$table
prop.contrib <- data.frame(spp = colnames(out$uni.test),
                           prop = out$uni.test[3, ],
                           row.names = NULL)
prop.contrib <- prop.contrib %>%
  mutate(perc = round((prop / sum(prop)) * 100, 1)) %>%
  arrange(desc(perc)) %>%
  mutate(cum = cumsum(perc))
