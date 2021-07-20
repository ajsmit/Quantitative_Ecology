# Quantitative Ecology (BCB743)
# Per-country life-expectancy and mortality due to various causes (2016)
# Data from the World Health Organization
# Author: AJ Smit
# Date: 18 July 2021


# Load packages -----------------------------------------------------------

library(tidyverse)


# Read data, filter, and subset -------------------------------------------

# SDG 1.a
# Domestic general government health expenditure (GGHE-D) as percentage of general government expenditure (GGE) (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/domestic-general-government-health-expenditure-(gghe-d)-as-percentage-of-general-government-expenditure-(gge)
SDG1.a <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG1.a_domestic_health_expenditure.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG1.a")

# SDG 2.2
#

# SDG 3.1
# Maternal mortality ratio (per 100 000 live births)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/maternal-mortality-ratio-(per-100-000-live-births)
SDG3.1_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.1_maternal_mort.csv") %>%
  filter(Period == 2016,
         Indicator == "Maternal mortality ratio (per 100 000 live births)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.1_1")

# SDG 3.1
# Births attended by skilled health personnel (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/births-attended-by-skilled-health-personnel-(-)
SDG3.1_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.1_skilled_births.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.1_2")

# SDG 3.2
# Number of neonatal deaths (Child mortality)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/number-of-neonatal-deaths
SDG3.2_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.2_neonatal_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_1")

# SDG 3.2
# Number of under-five deaths (Child mortality)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/number-of-under-five-deaths
SDG3.2_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.2_under_5_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_2")

# SDG 3.2
# Number of infant deaths (Child mortality)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/number-of-infant-deaths
SDG3.2_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.2_infant_deaths.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.2_3")

# SDG 3.3
# New HIV infections (per 1000 uninfected population)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/new-hiv-infections-(per-1000-uninfected-population)
SDG3.3_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.3_new_HIV_infections.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_1")

# SDG 3.3
# Incidence of tuberculosis (per 100 000 population per year)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/incidence-of-tuberculosis-(per-100-000-population-per-year)
SDG3.3_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.3_TB.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_2")

# SDG 3.3
# Malaria incidence (per 1 000 population at risk)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/malaria-incidence-(per-1-000-population-at-risk)
SDG3.3_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.3_malaria.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_3")

# SDG 3.3
# Hepatitis B surface antigen (HBsAg) prevalence among children under 5 years
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/hepatitis-b-surface-antigen-(hbsag)-prevalence-among-children-under-5-years
SDG3.3_4 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.3_hepatitis_B.csv") %>%
  filter(Period == 2015) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_4")

# SDG 3.3
# Reported number of people requiring interventions against NTDs
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/reported-number-of-people-requiring-interventions-against-ntds
SDG3.3_5 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.3_NCD_interventions.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.3_5")

# SDG 3.4
# Adult mortality rate (probability of dying between 15 and 60 years per 1000 population)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/adult-mortality-rate-(probability-of-dying-between-15-and-60-years-per-1000-population)
SDG3.4_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_adult_death_prob.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_1")

# SDG 3.4
# Number of deaths attributed to non-communicable diseases, by type of disease and sex
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/number-of-deaths-attributed-to-non-communicable-diseases-by-type-of-disease-and-sex
SDG3.4_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Diabetes mellitus") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_2")

SDG3.4_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Cardiovascular diseases") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_3")

SDG3.4_4 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_NCD_by_cause.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes",
         Dim2 == "Respiratory diseases") %>%
  mutate(Indicator = Dim2) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_4")

# SDG 3.4
# Crude suicide rates (per 100 000 population) (SDG 3.4.2)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/crude-suicide-rates-(per-100-000-population)
SDG3.4_5 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_suicides.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_5")

# SDG3.4
# Total NCD Deaths (in thousands)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/gho-ghe-ncd-deaths-in-thousands
SDG3.4_6 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.4_NCD_data_total.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.4_6")

# SDG 3.5
# Alcohol, total per capita (15+) consumption (in litres of pure alcohol) (SDG Indicator 3.5.2)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/total-(recorded-unrecorded)-alcohol-per-capita-(15-)-consumption
SDG3.5 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.5_alcohol_consumption.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.5")

# SDG 3.6
# Estimated road traffic death rate (per 100 000 population)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/estimated-road-traffic-death-rate-(per-100-000-population)
SDG3.6 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.6_traffic_deaths_prop.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.6")

# SDG 3.7
# Adolescent birth rate (per 1000 women aged 15-19 years)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/adolescent-birth-rate-(per-1000-women-aged-15-19-years)
SDG3.7 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.7_adolescent_births.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.7")

# SDG 3.8
# UHC Index of service coverage (SCI)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/uhc-index-of-service-coverage
SDG3.8_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.8_UHC_data_availability.csv") %>%
  filter(Period == "2013-2017") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.8_1")

# SDG 3.8
# Data availability for UHC index of essential service coverage (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/data-availability-for-uhc-index-of-essential-service-coverage-(-)
SDG3.8_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.8_UHC_index_of_service_coverage.csv") %>%
  filter(Period == 2017) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.8_2")

# SDG 3.8
# Population with household expenditures on health greater than 10% of total household expenditure or income (SDG 3.8.2) (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/population-with-household-expenditures-on-health-greater-than-10-of-total-household-expenditure-or-income-(sdg-3-8-2)-(-)
# SDG3.8_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.8_UHC_percent_of_expenditure_1.csv") %>%
#   filter(Period == 2016) %>%
#   select(Indicator, ParentLocation, Location, FactValueNumeric)

# SDG 3.8
# Population with household expenditures on health greater than 25% of total household expenditure or income ( SDG indicator 3.8.2) (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/population-with-household-expenditures-on-health-greater-than-25-of-total-household-expenditure-or-income-(-sdg-indicator-3-8-2)-(-)
# SDG3.8_4 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.8_UHC_percent_of_expenditure_2.csv") %>%
#   filter(Period == 2016) %>%
#   select(Indicator, ParentLocation, Location, FactValueNumeric)

# SDG 3.9
# Poison control and unintentional poisoning
# https://www.who.int/data/gho/data/themes/topics/indicator-groups/poison-control-and-unintentional-poisoning
SDG3.9_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.9_unintentional_poisoning_prop.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.9_1")

# SDG 3.9
# Indicator 3.9.1: Mortality rate attributed to household and ambient air pollution (per 100 000 population)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/ambient-and-household-air-pollution-attributable-death-rate-(per-100-000-population)
# SDG3.9_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.9_ambient_air_pollution.csv") %>%
#   filter(Period == 2016,
#          Dim1 == "Both sexes") %>%
#   select(Indicator, ParentLocation, Location, FactValueNumeric)

# SDG 3.9
# Mortality rate attributed to exposure to unsafe WASH services (per 100 000 population) (SDG 3.9.2)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/mortality-rate-attributed-to-exposure-to-unsafe-wash-services-(per-100-000-population)-(sdg-3-9-2)
SDG3.9_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.9_WASH_mortalities.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.9_3")

# SDG 16.1
# Estimates of rate of homicides (per 100 000 population)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/estimates-of-rates-of-homicides-per-100-000-population
SDG16.1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG16.1_homicides.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG16.1")

# SDG 3.a
# Prevalence of current tobacco use among persons aged 15 years and older (age-standardized rate)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/age-standardized-prevalence-of-current-tobacco-smoking-among-persons-aged-15-years-and-older
SDG3.a <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.a_tobacco_control.csv") %>%
  filter(Period == 2016,
         Dim1 == "Both sexes") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.a")

# SDG 3.b
# Total net official development assistance to medical research and basic health sectors per capita (US$), by recipient country
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/total-net-official-development-assistance-to-medical-research-and-basic-health-sectors-per-capita
SDG3.b_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.b_dev_assistence_for_med_research.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_1")

# SDG 3.b
# Measles-containing-vaccine second-dose (MCV2) immunization coverage by the nationally recommended age (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/measles-containing-vaccine-second-dose-(mcv2)-immunization-coverage-by-the-nationally-recommended-age-(-)
SDG3.b_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.b_measles_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_2")

# SDG 3.b
# Diphtheria tetanus toxoid and pertussis (DTP3) immunization coverage among 1-year-olds (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/diphtheria-tetanus-toxoid-and-pertussis-(dtp3)-immunization-coverage-among-1-year-olds-(-)
SDG3.b_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.b_diphtheria_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_3")

# SDG 3.b
# Pneumococcal conjugate vaccines (PCV3) immunization coverage among 1-year-olds (%)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/pneumoccocal-conjugate-vaccines-(pcv3)-immunization-coverage-among-1-year-olds-(-)
SDG3.b_4 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.b_pneumococcal_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_4")

# SDG 3.b
# Girls aged 15 years old that received the recommended doses of HPV vaccine
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/girls-aged-15-years-old-that-received-the-recommended-doses-of-hpv-vaccine
SDG3.b_5 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.b_HPV_vaccine.csv") %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.b_5")

# SDG 3.b
# Proportion of health facilities with a core set of relevant essential medicines available and affordable on a sustainable basis
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/proportion-of-health-facilities-with-a-core-set-of-relevant-essential-medicines-available-and-affordable-on-a-sustainable-basis
# Full data not available

# SDG 3.c
# SDG Target 3.c | Health workforce: Substantially increase health financing and the recruitment, development, training and retention of the health workforce in developing countries, especially in least developed countries and small island developing States
# https://www.who.int/data/gho/data/themes/topics/indicator-groups/indicator-group-details/GHO/sdg-target-3.c-health-workforce
SDG3.c_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Medical doctors (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_1")

SDG3.c_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Nursing and midwifery personnel (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_2")

SDG3.c_3 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Dentists (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_3")

SDG3.c_4 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.c_health_workforce.csv")  %>%
  filter(Period == 2016,
         Indicator == "Pharmacists  (per 10,000)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.c_4")

# SDG 3.d
# SDG Target 3.d | National and global health risks: Strengthen the capacity of all countries, in particular developing countries, for early warning, risk reduction and management of national and global health risks
# https://www.who.int/data/gho/data/themes/topics/indicator-groups/indicator-group-details/GHO/sdg-target-3.d-national-and-global-health-risks
# Data not available

# SDG 3.d
# Average of 13 International Health Regulations core capacity scores, SPAR version
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/-average-of-13-international-health-regulations-core-capacity-scores-spar-version
SDG3.d_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_SDG3.d_health_risks.csv")  %>%
  filter(Period == 2016) %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "SDG3.d_1")

# Other
# Life expectancy at birth (years)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/life-expectancy-at-birth-(years)
other_1 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_Other_life_expectancy.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes",
         Indicator == "Life expectancy at birth (years)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "other_1")

# Other
# Life expectancy at age 60 (years)
# https://www.who.int/data/gho/data/indicators/indicator-details/GHO/life-expectancy-at-age-60-(years)
other_2 <- read.csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_Other_life_expectancy.csv") %>%
  filter(Period == 2015,
         Dim1 == "Both sexes",
         Indicator == "Life expectancy at age 60 (years)") %>%
  select(Indicator, ParentLocation, Location, FactValueNumeric) %>%
  mutate(SDG = "other_2")


# rbind the data ----------------------------------------------------------

health <- do.call("rbind", lapply(ls(), get))


# Create list of SDGs used ------------------------------------------------

unique(health[, c(5, 1)])
write_csv(unique(health[, c(5, 1)]), file = "/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/SDG_description.csv")


# Listing of mortality indicators -----------------------------------------

unique(health$Indicator)

# [1] "Mortality rate attributed to unintentional poisoning (per 100 000 population)"
# [2] "Estimates of rates of homicides per 100 000 population"
# [3] "Data availability for UHC index of essential service coverage (%)"
# [4] "Total net official development assistance to medical research and basic health sectors per capita (US$), by recipient country"
# [6] "Measles-containing-vaccine second-dose (MCV2) immunization coverage by the nationally recommended age (%)"
# [7] "Ambient and household air pollution attributable death rate (per 100 000 population)"
# [8] "Ambient and household air pollution attributable death rate (per 100 000 population, age-standardized)"
# [9] "UHC index of essential service coverage"
# [11] "Diphtheria tetanus toxoid and pertussis (DTP3) immunization coverage among 1-year-olds (%)"
# [12] "Mortality rate attributed to exposure to unsafe WASH services (per 100 000 population) (SDG 3.9.2)"
# [13] "Population with household expenditures on health greater than 10% of total household expenditure or income (SDG 3.8.2) (%)"
# [15] "Pneumoccocal conjugate vaccines (PCV3) immunization coverage among 1-year-olds (%)"
# [16] "Girls aged 15 years old that received the recommended doses of HPV vaccine (%)"
# [17] "Life expectancy at birth (years)"
# [18] "Average of 13 International Health Regulations core capacity scores"
# [19] "Maternal mortality ratio (per 100 000 live births)"
# [22] "Life expectancy at age 60 (years)"
# [23] "Births attended by skilled health personnel (%)"
# [25] "Crude suicide rates (per 100 000 population)"
# [26] "Alcohol, total per capita (15+) consumption (in litres of pure alcohol) (SDG Indicator 3.5.2)"
# [27] "Estimated road traffic death rate (per 100 000 population)"
# [28] "Adolescent birth rate (per 1000 women aged 15-19 years)"
# [29] "Medical doctors (per 10,000)"
# [30] "New HIV infections (per 1000 uninfected population)"
# [31] "Domestic general government health expenditure (GGHE-D) as percentage of general government expenditure (GGE) (%)"
# [32] "Nursing and midwifery personnel (per 10,000)"
# [33] "Incidence of tuberculosis (per 100 000 population per year)"
# [34] "Dentists (per 10,000)"
# [35] "Malaria incidence (per 1 000 population at risk)"
# [36] "Pharmacists  (per 10,000)"
# [37] "Hepatitis B surface antigen (HBsAg) prevalence among children under 5 years (%)"

# To standardise to unit population size
# [5] "Number of neonatal deaths"
# [10] "Number of under-five deaths"
# [14] "Number of infant deaths"
# [20] "Diabetes mellitus"
# [21] "Age-standardized prevalence of current tobacco smoking among persons aged 15 years and older"
# [24] "Total NCD Deaths (in thousands)"
# [38] "Reported number of people requiring interventions against NTDs"


# Pivot wider -------------------------------------------------------------

health_wide <- health %>%
  arrange(Location) %>%
  select(-Indicator) %>%
  pivot_wider(names_from = SDG, values_from = FactValueNumeric) %>%
  as_tibble()


# Add world population data -----------------------------------------------

popl <- read_csv("/Users/ajsmit/Dropbox/R/workshops/Quantitative_Ecology/exercises/WHO/WHO_population.csv") %>%
  filter(Year == 2016) %>%
  rename(popl_size = `Population (in thousands) total`,
         Location = Country) %>%
  select(Location, popl_size) %>%
  mutate(popl_size = as.numeric(gsub("[[:space:]]", "", popl_size)) * 1000)

health_wide <- health_wide %>%
  left_join(popl)


# Standardise to popl size ------------------------------------------------

# standardise the remaining unstandardised variables
# *here normalised to unit population size
health_wide <- health_wide %>%
  mutate(SDG3.4_4 = SDG3.4_4 / popl_size * 100000,
         SDG3.4_3 = SDG3.4_3 / popl_size * 100000,
         SDG3.4_2 = SDG3.4_2 / popl_size * 100000,
         SDG3.4_6 = SDG3.4_6 / 100,
         SDG3.2_2 = SDG3.2_2 / popl_size * 100000,
         SDG3.2_3 = SDG3.2_3 / popl_size * 100000,
         SDG3.2_1 = SDG3.2_1 / popl_size * 100000)


# Correlations ------------------------------------------------------------

# a histogram of missing value occurrences
health_wide$na_count <- apply(health_wide[, 3:(ncol(health_wide) - 1)], 1, function(x) sum(is.na(x)))
hist(health_wide$na_count, breaks = 14, plot = TRUE)

# remove rows where there are more than 10 NAs
health_wide <- health_wide %>%
  filter(na_count <= 10) %>%
  select(-na_count)

# compute a correlation matrix
corr <- round(cor(health_wide[, 3:(ncol(health_wide) - 1)]), 1)

# visualization of the correlation matrix
ggcorrplot(corr, type = 'upper', outline.col = "grey60",
           colors = c("#1679a1", "white", "#f8766d"),
           lab = TRUE)


# Impute remaining NAs ----------------------------------------------------

# there are still many remaining NAs
# I impute them with the imputePCA() method in the msssMDA package
# see Dray and Josse (2015) Principal component analysis with missing values:
# a comparative survey of methods. Plant Ecol 216: 657-667
library(missMDA)

health_wide_complete <- imputePCA(health_wide[, 3:(ncol(health_wide) - 1)])$completeObs


# Scale and center the data -----------------------------------------------

health_wide_complete_std <- decostand(health_wide_complete, method = "standardize")


# Do the PCA --------------------------------------------------------------

health_pca <- rda(health_wide_complete_std)
summary(health_pca)

par(mfrow = c(1, 2))
biplot(health_pca, scaling = 1, main = "PCA scaling 1", choices = c(1, 2))
biplot(health_pca, scaling = 2, main = "PCA scaling 2", choices = c(1, 2))

pl1 <- ordiplot(health_pca, type = "none", scaling = 1, main = "PCA WHO/SDG")
points(pl1, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl1, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl1, "species", col = "blue4", cex = 0.9)
text(pl1, "sites", col = "red4", cex = 0.9)

pl2 <- ordiplot(health_pca, type = "none", scaling = 2, main = "PCA WHO/SDG")
points(pl2, "sites", pch = 21, cex = 1.75, col = "grey80", bg = "grey80")
points(pl2, "species", pch = 21, col = "turquoise", arrows = TRUE)
text(pl2, "species", col = "blue4", cex = 0.9)
text(pl2, "sites", col = "red4", cex = 0.9)

site_scores <- tibble(ParentLocation = health_wide$ParentLocation,
                      Location = health_wide$Location)
site_scores <- tibble(cbind(site_scores, scores(health_pca, display = "sites", choices = c(1:7))))
species_scores <- data.frame(scores(health_pca, display = "species", choices = c(1:7)))
species_scores$species <- rownames(species_scores)
species_scores <- tibble(species_scores)

ggplot(data = site_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = ParentLocation)) +
  geom_segment(data = species_scores,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.4, "cm"), type = "closed"),
               color = "lightseagreen", alpha = 1, size = 0.3) +
  geom_text(data = species_scores,
            aes(x = PC1, y = PC2, label = species),
            color = "black") +
  xlab("PC1") + ylab("PC2") +
  ggtitle("WHO SDGs, Scaling 2")

