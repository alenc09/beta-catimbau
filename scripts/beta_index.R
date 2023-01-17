# Script to calculate descriptive beta-diversity metrics

#libraries----
library(readxl)
library(here)
library(dplyr)
library(vegan)
library(betapart)

#data----
read_xlsx(here("data/herbáceas_pp.xlsx"))-> herbs
read_xlsx(here("data/lenhosas_pp.xlsx"))-> len

#analysis----
##beta diversity metrics----
herbs %>% 
  mutate(parcela = as.factor(parcela)) %>% 
  rename(total = `Total Geral`) %>% 
  glimpse -> herbs

beta.pair
