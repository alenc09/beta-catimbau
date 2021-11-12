# Wed Nov 10 17:58:45 2021 ------------------------------
#script para análises usando dados de abundância

#data----
library(betapart)
library(vegan)

# dbRDA Analysis ----
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
    #PCNM2 + #PCNM with the highest VIF value. We cut this var. and recheck if the VIF values of other variables
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
set.seed(123)
ordistep(
  capscale(herb.abund.tot ~ 1, plot_pcnm_transf),
  scope = formula(modT2.abund.herb),
  direction = 'both',
  pstep = 1000
) 

capscale(herb.abund.tot ~ 
           PCNM3 + 
           LPI + 
           prec, 
         plot_pcnm_transf) -> mod.herb.abund.tot #model with variables selected by stepwise selection
anova(mod.herb.abund.tot)
RsquareAdj(mod.herb.abund.tot)


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

### Total beta-diversity----
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

### Turnover component-----
set.seed(123)
ordistep(capscale(len.abund.tu~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
capscale(len.abund.tu ~ WEI + prec + PCNM2 + PCNM5, plot_pcnm_transf) -> mod.len.abund.tu
anova(mod.len.abund.tu)
RsquareAdj(mod.len.abund.tu)

### Nestedness component ----
ordistep(capscale(len.abund.ne~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
capscale(len.abund.ne ~ PCNM3, plot_pcnm_transf) -> mod.len.abund.ne
anova(mod.len.abund.ne)
RsquareAdj(mod.len.abund.ne)

# Variation partitioning analysis----
## herbs ----
### total bet-diversity ----
anova(capscale(
  herb.abund.tot ~ 
    PCNM3 +
    PCNM1 + 
    PCNM10 + 
    Condition(LPI + prec), #only spatial effects
  data = plot_pcnm_transf
)) #Analysis to check the significance of each variable group in explaining the beta-diversity components
anova(capscale(
  herb.abund.tot ~ 
    LPI  + 
    Condition(PCNM3 + PCNM1 + PCNM10 + prec), #only disturbance effect
  data = plot_pcnm_transf
))
anova(capscale(
  herb.abund.tot ~ 
    prec + 
    Condition(PCNM3 + PCNM1 + PCNM10 + LPI), #only precipitation effect
  data = plot_pcnm_transf
))
varpart(herb_abund_hell,                    #Proper partitioning analysis
        plot_pcnm_transf[, 2],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(8, 10, 17)]) -> varpart.abund.herb.tot #herb.abund.tot ~ LPI + prec + PCNM1 + PCNM3 + PCNM10
plot(varpart.abund.herb.tot)

### Turnover component ----
anova(capscale(herb.abund.tu ~ 
                 PCNM3 + 
                 PCNM1 + 
                 PCNM5 + 
                 Condition(prec),
               data = plot_pcnm_transf))
anova(capscale(herb.abund.tu ~ 
                 prec + 
                 Condition(PCNM3 + PCNM1 + PCNM5), data = plot_pcnm_transf))

varpart(herb_abund_hell,
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(8, 10, 12)]) -> varpart.abund.herb.tu #herb.abund.tu ~ prec + PCNM3 + PCNM1 + PCNM5
plot(varpart.abund.herb.tu)

### Nestedness component ----
#not necessary because it was only explained by spatial variables

## Wood ----
### Total beta-diversity ----
anova(capscale(len.abund.tot ~ 
                 PCNM2 +
                 PCNM5 +
                 Condition(WEI + prec),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tot ~ 
                 WEI + 
                 Condition(prec + PCNM2 + PCNM5),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tot ~ 
                 prec + 
                 Condition(WEI + PCNM2 + PCNM5),
               data = plot_pcnm_transf))

varpart(len_abund_hell, 
        plot_pcnm_transf[,3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(9,12)]) -> varpart.abund.len.tot #len.abund.tot ~ WEI + prec + PCNM2 + PCNM5
plot(varpart.abund.len.tot)

### Turnover component ----
anova(capscale(len.abund.tu ~ 
                 PCNM2 +
                 PCNM5 +
                 Condition(WEI + prec),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~ 
                 WEI + 
                 Condition(prec + PCNM2 + PCNM5),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~ 
                 prec + 
                 Condition(WEI + PCNM2 + PCNM5),
               data = plot_pcnm_transf))

varpart(len_abund_hell, 
        plot_pcnm_transf[,3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(9,12)]) -> varpart.abund.len.tu #len.abund.tot ~ WEI + prec + PCNM2 + PCNM5
plot(varpart.abund.len.tu)

### Nestedness component ----
#not necessary because it was only explained by spatial variables