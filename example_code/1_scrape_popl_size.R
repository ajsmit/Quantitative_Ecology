# Scrape population data from https://www.worldometers.info/world-population/population-by-country/

library(tidyverse)
library(rvest)

webpage <- read_html("https://www.worldometers.info/world-population/population-by-country/")
tables <- webpage %>%
  html_table(fill = TRUE)
world_pop <- tables[[1]] %>%
  select(-1) %>%
  clean_names() %>%
  rename(Location = country_or_dependency) %>%
  mutate(population_2020 = as.numeric(gsub(",", "", population_2020)),
         yearly_change = as.numeric(gsub(" %", "", yearly_change)),
         net_change = as.numeric(gsub(",", "", net_change)),
         density_p_km2 = as.numeric(gsub(",", "", density_p_km2)),
         land_area_km2 = as.numeric(gsub(",", "", land_area_km2)),
         migrants_net = as.numeric(gsub(",", "", migrants_net)),
         fert_rate = as.numeric(fert_rate),
         med_age = as.numeric(med_age),
         urban_pop_percent = as.numeric(gsub(" %", "", urban_pop_percent)),
         world_share = as.numeric(gsub(" %", "", world_share)))

# Edit names of countries/regions to match that in the JHU data
world_pop[world_pop == "South Korea"] <- "Democratic People's Republic of Korea" # correct
world_pop[world_pop == "Taiwan"] <- "Taiwan*" # ???
world_pop[world_pop == "United States"] <- "United States of America" # correct
world_pop[world_pop == "Czech Republic (Czechia)"] <- "Czechia" # correct
# West Bank and Gaza ???
# 'Congo (Kinshasa)' in JSU data changed to 'DR Congo'
world_pop[world_pop == "Côte d'Ivoire"] <- "Côte d’Ivoire" # correct
world_pop[world_pop == "St. Vincent & Grenadines"] <- "Saint Vincent and the Grenadines" # correct
# 'Congo (Brazzaville)' in JSU data changed to 'DR Congo'
world_pop[world_pop == "Saint Kitts & Nevis"] <- "Saint Kitts and Nevis" # correct
# 'Kosovo' and 'Serbia' combined in JSU data
world_pop[world_pop == "Myanmar"] <- "Burma" # correct
world_pop[world_pop == "Sao Tome & Principe"] <- "Sao Tome and Principe" # correct
world_pop[world_pop == "Venezuela"] <- "Venezuela (Bolivarian Republic of)" # correct

# Serbia in mortalities data


# Update % urban population from https://data.worldbank.org/indicator/SP.URB.TOTL.in.zs
# except where noted otherwise
world_pop %>%
  select(country_or_dependency, urban_pop_percent) %>%
  filter(is.na(urban_pop_percent))

world_pop$urban_pop_percent[world_pop$country_or_dependency == "Venezuela"] <- 88
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Hong Kong"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Singapore"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Kuwait"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Puerto Rico"] <- 94
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Macao"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Guadeloupe"] <- 98.4 # updated from http://data.un.org/
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Monaco"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Gibraltar"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Palau"] <- 80
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Anguilla"] <- 100 # updated from https://www.cia.gov/
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Nauru"] <- 100
world_pop$urban_pop_percent[world_pop$country_or_dependency == "Holy See"] <- 100 # updated from https://www.cia.gov/
