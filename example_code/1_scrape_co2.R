# Scrape CO2 emissions data from https://www.worldometers.info/co2-emissions/

library(rvest)
library(tidyverse)
library(janitor)

webpage <- read_html("https://www.worldometers.info/co2-emissions/")
tables <- webpage %>%
  html_table(fill = TRUE)
co2_table <- tables[[1]] %>%
  select(-1) %>%
  clean_names()

# Edit names of countries/regions to match that in the JHU data
co2_table[co2_table == "South Korea"] <- "Korea, South"
co2_table[co2_table == "Taiwan"] <- "Taiwan*"
co2_table[co2_table == "United States"] <- "US"
co2_table[co2_table == "Czech Republic (Czechia)"] <- "Czechia"
# West Bank and Gaza ???
# 'Congo (Kinshasa)' in JSU data changed to 'DR Congo'
co2_table[co2_table == "Côte d'Ivoire"] <- "Cote d'Ivoire"
co2_table[co2_table == "St. Vincent & Grenadines"] <- "Saint Vincent and the Grenadines"
# 'Congo (Brazzaville)' in JSU data changed to 'DR Congo'
co2_table[co2_table == "Saint Kitts & Nevis"] <- "Saint Kitts and Nevis"
# 'Kosovo' and 'Serbia' combined in JSU data
co2_table[co2_table == "Myanmar"] <- "Burma"
co2_table[co2_table == "Sao Tome & Principe"] <- "Sao Tome and Principe"
