# Fri Nov 10 11:34:59 2023 ------------------------------

#libraries----
library(readxl)
library(geosphere)
library(vegan)
library(betapart)

#data----
read_xlsx("data/lenhosas_pp.xlsx") -> len
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot

##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

#data transformation####
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(len[-3,-1], method = "pa")-> len_pa
decostand(len_pa, 'hellinger') -> len_pa_hell

beta.pair(len_pa[-3,-1]) -> wood.pair.pa
wood.pair.pa$beta.sor -> wood.pa.tot
wood.pair.pa$beta.sim -> wood.pa.tu
wood.pair.pa$beta.sne -> wood.pa.ne

### Total beta-diversity ----
anova(capscale(wood.pa.tot ~ 
                 PCNM2 +
                 PCNM10 +
                 PCNM5 +
                 Condition(WEI + prec), #only spatial effects
               data = plot_pcnm_transf[-3,])) #remove the extra plot wood vegetation data has

anova(capscale(wood.pa.tot ~ 
                 WEI + 
                 Condition(prec + PCNM2 + PCNM10 + PCNM5), #only CAD
               data = plot_pcnm_transf[-3,]))

anova(capscale(wood.pa.tot ~ 
                 prec + 
                 Condition(WEI + PCNM2 + PCNM10 + PCNM5), #only precipitation
               data = plot_pcnm_transf[-3,]))

varpart(len_pa_hell, 
        plot_pcnm_transf[, 3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(9,17,12)]) -> varpart.pa.len.tot #len.abund.tot ~ WEI + prec + PCNM2 + PCNM10 + PCNM5
summary(varpart.pa.len.tot)
plot(varpart.abund.len.tot)

### Turnover component ----
anova(capscale(wood.pa.tot ~ 
                 PCNM2 +
                 Condition(WEI + prec),
               data = plot_pcnm_transf[-3,]))

anova(capscale(wood.pa.tot ~ 
                 WEI + 
                 Condition(prec + PCNM2),
               data = plot_pcnm_transf[-3,]))

anova(capscale(wood.pa.tot ~ 
                 prec + 
                 Condition(WEI + PCNM2),
               data = plot_pcnm_transf[-3,]))

varpart(len_pa_hell, 
        plot_pcnm_transf[,3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, 9]) -> varpart.pa.len.tu #len.abund.tot ~ WEI + prec + PCNM2 + PCNM5
summary(varpart.pa.len.tu)
plot(varpart.pa.len.tu)

### Nestedness component ----
#not necessary because it was only explained by spatial variables