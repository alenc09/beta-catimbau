# Fri Nov 10 10:41:23 2023 ------------------------------

#libraries----
library(readxl)
library(geosphere)
library(vegan)
library(betapart)
library(dplyr)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/herbáceas_pp.xlsx") -> herb

##manipulation----
herb -> herb_pa
herb_pa[,-1][herb_pa[,-1] >= 1] <- 1
glimpse(herb_pa)

##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##data transformation----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(herb_pa[,-c(1,59)], 'hellinger')-> herb_pa_hell

##Herb - dbRDA data----
herb.pair.pa<- beta.pair(herb_pa[,-1])
herb.pa.tot<- herb.pair.pa$beta.sor
herb.pa.tu<- herb.pair.pa$beta.sim
herb.pa.ne<- herb.pair.pa$beta.sne

### total beta-diversity ---- #Models were defined in previous analysis (dbRDA_herb_pa.R)
anova(capscale(
  herb.pa.tot ~ 
    PCNM1 +
    PCNM3 +
    PCNM5 +
    Condition(prec), #only spatial effects
  data = plot_pcnm_transf)) #Analysis to check the significance of each variable group in explaining the beta-diversity components

anova(capscale(
  herb.pa.tot ~ 
    prec + 
    Condition(PCNM1 + PCNM3 + PCNM5), #only precipitation effect
  data = plot_pcnm_transf))

varpart(herb_pa_hell,                    #Proper partitioning analysis
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(8, 10, 12)]) -> varpart.pa.herb.tot #herb.abund.tot ~ prec + PCNM1 + PCNM3 + PCNM5
summary(varpart.pa.herb.tot)
plot(varpart.pa.herb.tot)

### Turnover component ----
#not necessary because it was only explained by spatial variables (PCNM 5)

### Nestedness component ----
#not necessary because it was only explained by spatial variables (PCNM1 + PCNM3 + PCNM7 + PCNM10 + PCNM8)