# Fri Oct 27 12:23:10 2023 ------------------------------
# Script para análises de beta-diversidade de herbaceas com dados de presença/ausencia

#Libraries----
library(readxl)
library(dplyr)
library(geosphere)
library(vegan)
library(betapart)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/lenhosas_pp.xlsx") -> wood

##manipulation----
wood -> wood_pa
wood_pa[,-1][wood_pa[,-1] >= 1] <- 1
glimpse(wood_pa)

#analysis----
##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##standardize----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf
decostand(wood_pa[-3,-1], 'hellinger') -> wood_pa_hell

##RDA model----
set.seed(123)
rda(wood_pa_hell ~ ., data = plot_pcnm_transf) -> RDAmod_wood_pa
vif.cca(RDAmod_wood_pa)

rda(
  wood_pa_hell ~ 
    PPI +
    LPI +
    WEI +
    #alt + #highest VIF value
    prec +
    fert_sol +
    #PCNM1 + #highest VIF value among PCNMs
    PCNM2 + 
    PCNM3 +
    PCNM4 +
    PCNM5 +
    PCNM6 +
    PCNM7 +
    PCNM8 + 
    PCNM9 + 
    PCNM10 + 
    PCNM11,
  data = plot_pcnm_transf
) -> RDAmod2_wood_pa
vif.cca(RDAmod2_wood_pa) #all VIF values are below 10. 

##dbRDA----
beta.pair(wood_pa[-3,-1]) -> wood.pair.pa
wood.pair.pa$beta.sor -> wood.pa.tot
wood.pair.pa$beta.sim -> wood.pa.tu
wood.pair.pa$beta.sne -> wood.pa.ne

## Total beta-diversity----
set.seed(123)
ordistep(
  capscale(wood.pa.tot ~ 1, plot_pcnm_transf),
  scope = formula(RDAmod2_wood_pa),
  direction = 'both',
  pstep = 1000)

capscale(wood.pa.tot ~ prec + WEI + PCNM2 + PCNM10 + PCNM5, plot_pcnm_transf) -> mod.wood.pa.tot
anova(mod.wood.pa.tot)
RsquareAdj(mod.wood.pa.tot)

## Turnover component-----
set.seed(123)
ordistep(
  capscale(wood.pa.tu ~ 1, plot_pcnm_transf),
  scope = formula(RDAmod2_wood_pa),
  direction = 'both',
  pstep = 1000)

capscale(wood.pa.tu ~ WEI + prec + PCNM2, plot_pcnm_transf) -> mod.wood.pa.tu
anova(mod.wood.pa.tu)
RsquareAdj(mod.wood.pa.tu)

## Nestedness component ----
ordistep(
  capscale(wood.pa.ne ~ 1, plot_pcnm_transf),
  scope = formula(RDAmod2_wood_pa),
  direction = 'both',
  pstep = 1000)

capscale(wood.pa.ne ~ PCNM11, plot_pcnm_transf) -> mod.wood.pa.ne
anova(mod.wood.pa.ne)
RsquareAdj(mod.wood.pa.ne)

