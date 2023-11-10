# Script to calculate descriptive beta-diversity metrics

#libraries----
library(readxl)
library(here)
library(dplyr)
library(vegan)
library(betapart)
library(entropart)
library(data.table)

#data----
read_xlsx(here("data/herbáceas_pp.xlsx"))-> herb
read_xlsx(here("data/lenhosas_pp.xlsx"))-> wood

#analysis----
##beta diversity metrics----
###herb----
#### abundance-based dissimilarity indexes----
herb %>% 
  mutate(parcela = as.factor(parcela)) %>% 
  rename(total = `Total Geral`) %>% 
  glimpse -> herb

decostand(herb[,-c(1,59)], 'hellinger') -> herb_hell #First line is the column names and 59 is the total abundance of each site

beta.pair.abund(herb_hell, index.family = "bray")-> herb.bray #bray-curtis dissimilarity index
herb.bray$beta.bray -> herb.beta.tot #total beta-diversity
herb.bray$beta.bray.bal -> herb.beta.tu #turnover share
herb.bray$beta.bray.gra -> herb.beta.ne #nestedness share

#### true-diversity indexes----
herb %>% 
  select(-total) %>% 
  transpose(make.names = T) %>% 
  glimpse %>% 
  MetaCommunity()-> mc.herb 

summary(mc.herb)  
plot(mc.herb)

BetaDiversity(mc.herb, q = 0, Correction = "None") %>% 
summary()

BetaDiversity(mc.herb, q = 1, Correction = "None") %>% 
summary()

BetaDiversity(mc.herb,q = 2, Correction = "None") %>% 
summary()


###wood----
#### abundance-based dissimilarity indexes----
wood %>% 
  mutate(parcela = as.factor(parcela)) %>% 
  rename(total = `Total Geral`) %>% 
  glimpse -> wood

decostand(wood[,-c(1,59)], 'hellinger') -> wood_hell #First column is the plot names and 59 is the total abundance of each site

beta.pair.abund(wood_hell, index.family = "bray")-> wood.bray #bray-curtis dissimilarity index
wood.bray$beta.bray -> wood.beta.tot #total beta-diversity
wood.bray$beta.bray.bal -> wood.beta.tu #turnover share
wood.bray$beta.bray.gra -> wood.beta.ne #nestedness share

#figures----
##herb----
### dissimilarity ----

###true diversity