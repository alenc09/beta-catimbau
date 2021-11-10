# Wed Nov 10 17:58:45 2021 ------------------------------
#script para análises usando dados de abundância

#data----
library(betapart)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

#Analysis----
## Herbs----
set.seed(123)
rda(herb_abund_hell ~ ., data = plot_pcnm_transf) -> modT.abund.herb #base model with all variabels
vif.cca(modT.abund.herb)
rda(
  herb_abund_hell ~ 
    LPI +
    prec +
    fert_sol +
    PCNM1 +
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
) -> modT2.abund.herb #Base model with VIF<10 variables, removed PCNM2 because, empirically, it was the PCNM with highest VIF
vif.cca(modT2.abund.herb) 

### Total beta-diversity ----
ordistep(
  capscale(herb.abund.tot ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.herb),
  direction = 'both',
  pstep = 1000
) 

capscale(herb.abund.tot ~ 
           PCNM3 + 
           LPI + 
           prec + 
           PCNM1 +
           PCNM10, 
         plot_pcnm_transf) -> mod.herb.abund.tot #model with variables selected by stepwise selection
anova(mod.herb.abund.tot)
RsquareAdj(mod.herb.abund.tot)
plot(mod.herb.abund.tot)
summary(mod.herb.abund.tot)

### Turnover component----
set.seed(123)
ordistep(
  capscale(herb.abund.tu ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.herb),
  direction = 'both',
  pstep = 1000
)

capscale(herb.abund.tu ~ 
           PCNM3 + 
           prec +
           PCNM5 +
           PCNM1, 
         data = plot_pcnm_transf) -> mod.herb.abund.tu #model with variables selected by stepwise selection

anova(mod.herb.abund.tu)
RsquareAdj(mod.herb.abund.tu)
plot(mod.herb.abund.tu)
summary(mod.herb.abund.tu)

###Nestedness component----
set.seed(123)
ordistep(
  capscale(herb.abund.ne ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.herb),
  direction = 'both',
  pstep = 1000
)
capscale(herb.abund.ne ~ 
           PCNM3 + 
           PCNM7 +
           PCNM8,
         data = plot_pcnm_transf) -> mod.herb.abund.ne
anova(mod.herb.abund.ne)
RsquareAdj(mod.herb.abund.ne)
plot(mod.herb.abund.ne)
summary(mod.herb.abund.ne)

##Wood----
rda(len_abund_hell ~ ., plot_pcnm_transf) -> modT.abund.len
vif.cca(modT.abund.len)

rda(
    len_abund_hell ~ 
      LPI + 
      WEI + 
      prec + 
      fert_sol + 
      PCNM1 + 
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

### Total beta-diversity----
set.seed(123)
ordistep(
  capscale(len.abund.tot ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.len),
  direction = 'both',
  pstep = 1000
)
capscale(len.abund.tot ~ PCNM1 + prec, plot_pcnm_transf) -> mod.len.abund.tot
anova(mod.len.abund.tot)
RsquareAdj(mod.len.abund.tot)

### Turnover component-----
#PAREI AQUI----
ordistep(capscale(len.abund.tu~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
mod.len.abund.tu<- capscale(len.abund.tu ~ WEI + prec + PCNM1 + PCNM5, plot_pcnm_transf)
anova(mod.len.abund.tu)
RsquareAdj(mod.len.abund.tu)