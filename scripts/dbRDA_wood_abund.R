# Tue Nov  1 10:56:42 2022 ------------------------------
#Script for dbRDA analysis on woody plant species

#Libraries----
library(readxl)
library(dplyr)
library(geosphere)
library(vegan)
library(betapart)

#data----
read_xlsx("data/lenhosas_pp.xlsx") -> wood
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot

#analysis----
##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

#data transformation####
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(wood[-3,-1], 'hellinger') -> wood_abund_hell

##dbRDA data----
beta.pair.abund(wood_abund_hell) ->wood.pair.abund
wood.pair.abund$beta.bray ->wood.abund.tot
wood.pair.abund$beta.bray.bal -> wood.abund.tu
wood.pair.abund$beta.bray.gra ->wood.abund.ne

#Analysis----
##dbRDA base model----
rda(wood_abund_hell ~ ., plot_pcnm_transf) -> modT.abund.wood
vif.cca(modT.abund.wood)

rda(
  wood_abund_hell ~ 
    PPI +
    LPI + 
    WEI + 
    #alt + #Highest vif value
    prec + 
    fert_sol + 
    #PCNM1 + #PCNM with the highest VIF value. We cut this var. and recheck if the VIF values of other variables
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
) ->modT2.abund.wood
vif.cca(modT2.abund.wood)

## Total beta-diversity----
set.seed(123)
ordistep(
  capscale(wood.abund.tot ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.wood),
  direction = 'both',
  pstep = 1000
)
capscale(wood.abund.tot ~ WEI + prec + PCNM2 + PCNM5, plot_pcnm_transf) -> mod.wood.abund.tot
anova(mod.wood.abund.tot)
RsquareAdj(mod.wood.abund.tot)

## Turnover component----
set.seed(123)
ordistep(capscale(wood.abund.tu~1, plot_pcnm_transf), scope = formula(modT2.abund.wood),
         direction = 'both', pstep=1000)
capscale(wood.abund.tu ~ WEI + prec + PCNM2 + PCNM5, plot_pcnm_transf) -> mod.wood.abund.tu
anova(mod.wood.abund.tu)
RsquareAdj(mod.wood.abund.tu)

## Nestedness component ----
ordistep(capscale(wood.abund.ne~1, plot_pcnm_transf), scope = formula(modT2.abund.wood),
         direction = 'both', pstep=1000)
capscale(wood.abund.ne ~ PCNM3, plot_pcnm_transf) -> mod.wood.abund.ne
anova(mod.wood.abund.ne)
RsquareAdj(mod.wood.abund.ne)

