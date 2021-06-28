library(vegan)

# Dataset 1
species <- data.frame(site = c("site_A", "site_B", "site_C", "site_D", "site_E", "site_F", "site_G", "site_H"),
           sp_A = c(1, 1, 4, 10, 0, 0, 1, 10),
           sp_B = c(1, 2, 4, 11, 0, 0, 1, 10),
           sp_C = c(1, 1, 5, 10, 0, 0, 1, 10),
           sp_D = c(2, 1, 4, 10, 0, 0, 1, 10),
           sp_E = c(1, 2, 5, 10, 1, 1, 1, 10),
           sp_F = c(10, 1, 4, 11, 1, 10, 1, 10))


species_pa <- decostand(species[, c(2:ncol(species))], method = "pa", MARGIN = 1)
species_pa <- cbind(species[, 1], species_pa)
colnames(species_pa) <- c("site", c(colnames(species[, 2:ncol(species)])))

specnumber(species[, 2:9], MARGIN = 2)
round(diversity(species[, 2:9], MARGIN = 2, index = "shannon"), 2)
round(diversity(species[, 2:9], MARGIN = 2, index = "simpson"), 2)

# Dataset 2 (plants and light levels)
light <- read.csv("exercises/diversity/light_levels.csv")
light_pa <- decostand(light[, c(2:7)], method = "pa", MARGIN = 1)
light_pa <- cbind(light[, 1], light_pa)
colnames(light_pa) <- c("site", c(colnames(light[, 2:7])))

light_div <- data.frame(
  site = c("low_light", "mid_light", "high_light"),
  richness = specnumber(light[, 2:6], MARGIN = 1),
  shannon = round(diversity(light[, 2:6], MARGIN = 1, index = "shannon"), 2),
  simpson = round(diversity(light[, 2:6], MARGIN = 1, index = "simpson"), 2)
)


# Dataset 3
plantations <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/data/plantations.csv")
plantations[is.na(plantations)] <- 0

specnumber(plantations[, 2:5], MARGIN = 2)
round(diversity(plantations[, 2:5], MARGIN = 2, index = "shannon"), 2)
round(diversity(plantations[, 2:5], MARGIN = 2, index = "simpson"), 2)

plantations[1:5, 1:5]
