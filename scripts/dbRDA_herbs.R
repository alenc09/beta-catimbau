# Wed Nov 10 17:58:45 2021 ------------------------------
#script para análises usando dados de abundância

#libraries----
library(betapart)
library(vegan)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/herbáceas_pp.xlsx") -> herb

##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##data transformation----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(herb[,-c(1,59)], 'hellinger')-> herb_abund_hell

##Herb - dbRDA data----
beta.pair.abund(herb_abund_hell) -> herb.pair.abund
herb.pair.abund$beta.bray -> herb.abund.tot
herb.pair.abund$beta.bray.bal -> herb.abund.tu
herb.pair.abund$beta.bray.gra -> herb.abund.ne

# dbRDA Analysis ----
## bray-curtir dissimilarity index----
vegdist(herb[,-c(1,59)], method = "bray") -> brayDiss_herb
range(brayDiss_herb)

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
