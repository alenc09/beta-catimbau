# Fri Oct 27 10:33:37 2023 ------------------------------
# Script para análises de beta-diversidade de herbaceas com dados de presença/ausencia

#Libraries----
library(readxl)
library(dplyr)
library(geosphere)
library(vegan)
library(betapart)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/herbáceas_pp.xlsx") -> herb

##manipulation----
herb -> herb_pa
herb_pa[,-1][herb_pa[,-1] >= 1] <- 1
glimpse(herb_pa)


#analysis----
##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##standardize----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf
decostand(herb_pa[,-c(1,59)], 'hellinger') -> herb_pa_hell

##RDA model----
set.seed(123)
rda(herb_pa_hell ~ ., data = plot_pcnm_transf) -> RDAmod_herb_pa
vif.cca(RDAmod_herb_pa)

rda(
  herb_pa_hell ~ #removed variables without theoretical relation with herbaceous beta-div (e.g. wood extraction)
    PPI +
    LPI +
    prec +
    fert_sol +
    PCNM1 +
    #PCNM2 + 
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
) -> RDAmod2_herb_pa
vif.cca(RDAmod2_herb_pa) #all VIF values are below 10. 

## dbRDA
herb.pair.pa<- beta.pair(herb_pa[,-1])
herb.pa.tot<- herb.pair.pa$beta.sor
herb.pa.tu<- herb.pair.pa$beta.sim
herb.pa.ne<- herb.pair.pa$beta.sne

###dbRDA - total----
set.seed(123)
ordistep(capscale(herb.pa.tot~1, plot_pcnm_transf), scope = formula(RDAmod2_herb_pa),
         direction = 'both', pstep=1000)
capscale(herb.pa.tot ~ PCNM1 + PCNM3 + PCNM5 + prec, plot_pcnm_transf) -> mod.herb.pa.tot
anova(mod.herb.pa.tot)
RsquareAdj(mod.herb.pa.tot)
plot(mod.herb.pa.tot)

###dbRDA - turnover
ordistep(capscale(herb.pa.tu~1, plot_pcnm_transf), scope = formula(RDAmod2_herb_pa),
         direction = 'both', pstep=1000)
mod.herb.pa.tu<- capscale(herb.pa.tu ~ PCNM5, plot_pcnm_transf)
anova(mod.herb.pa.tu)
RsquareAdj(mod.herb.pa.tu)
plot(mod.herb.pa.tu)

###dbRDA - nestedness----
ordistep(capscale(herb.pa.ne~1, plot_pcnm_transf), scope = formula(RDAmod2_herb_pa),
         direction = 'both', pstep=1000)
mod.herb.pa.ne<- capscale(herb.pa.ne ~ PCNM1 + PCNM3 + PCNM7 + PCNM10 + PCNM8, plot_pcnm_transf)
anova(mod.herb.pa.ne)
RsquareAdj(mod.herb.pa.ne)
plot(mod.herb.pa.ne)

