# Script to calculate descriptive beta-diversity metrics

#libraries----
library(readxl)
library(here)
library(dplyr)
library(vegan)
library(betapart)

#data----
read_xlsx(here("data/herbáceas_pp.xlsx"))-> herbs
read_xlsx(here("data/lenhosas_pp.xlsx"))-> wood

#analysis----
##beta diversity metrics----
###herbs----
#### abundance-based dissimilarity indexes----
herbs %>% 
  mutate(parcela = as.factor(parcela)) %>% 
  rename(total = `Total Geral`) %>% 
  glimpse -> herbs

decostand(herbs[,-c(1,59)], 'hellinger') -> herbs_hell #First line is the column names and 59 is the total abundance of each site

beta.pair.abund(herbs_hell, index.family = "bray")-> herbs.bray #bray-curtis dissimilarity index
herbs.bray$beta.bray -> herb.beta.tot #total beta-diversity
herbs.bray$beta.bray.bal -> herb.beta.tu #turnover share
herbs.bray$beta.bray.gra -> herb.beta.ne #nestedness share

#### true-diversity indexes----
