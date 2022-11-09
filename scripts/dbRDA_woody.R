# Tue Nov  1 10:56:42 2022 ------------------------------
#Script for dbRDA analysis on woody plant species

#Libraries----

#data----
read_xlsx("data/lenhosas_pp.xlsx") -> len
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot

#data transformation####
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(len[-3,-1], 'hellinger') -> len_abund_hell

##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##dbRDA data----
beta.pair.abund(len_abund_hell) ->len.pair.abund
len.pair.abund$beta.bray ->len.abund.tot
len.pair.abund$beta.bray.bal -> len.abund.tu
len.pair.abund$beta.bray.gra ->len.abund.ne

#Analysis----
##dbRDA base model----
rda(len_abund_hell ~ ., plot_pcnm_transf) -> modT.abund.len
vif.cca(modT.abund.len)

rda(
  len_abund_hell ~ 
    LPI + 
    WEI + 
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
) ->modT2.abund.len
vif.cca(modT2.abund.len)

## Total beta-diversity----
set.seed(123)
ordistep(
  capscale(len.abund.tot ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.len),
  direction = 'both',
  pstep = 1000
)
capscale(len.abund.tot ~ WEI + prec + PCNM2 + PCNM5, plot_pcnm_transf) -> mod.len.abund.tot
anova(mod.len.abund.tot)
RsquareAdj(mod.len.abund.tot)

## Turnover component-----
set.seed(123)
ordistep(capscale(len.abund.tu~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
capscale(len.abund.tu ~ WEI + prec + PCNM2 + PCNM5, plot_pcnm_transf) -> mod.len.abund.tu
anova(mod.len.abund.tu)
RsquareAdj(mod.len.abund.tu)

## Nestedness component ----
ordistep(capscale(len.abund.ne~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
capscale(len.abund.ne ~ PCNM3, plot_pcnm_transf) -> mod.len.abund.ne
anova(mod.len.abund.ne)
RsquareAdj(mod.len.abund.ne)