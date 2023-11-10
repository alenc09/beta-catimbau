# Fri Nov 10 13:07:14 2023 ------------------------------

#libraries----
library(readxl)
library(dplyr)
library(geosphere)
library(vegan)
library(betapart)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/lenhosas_pp.xlsx") -> wood

#Analysis----
##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##data transformation----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf 
decostand(wood[-3,-1], method = "pa") -> wood_pa
decostand(wood_pa, 'hellinger') -> wood_pa_hell

#Exponential decay models ----
##distance matrices ----
vegdist(plot_pcnm$LPI, method = "euclidean") -> dist.LPI
vegdist(plot_pcnm$prec, method = "euclidean") -> dist.prec
vegdist(plot_pcnm$WEI, method = "euclidean") -> dist.wei

vegdist(plot_pcnm_transf$LPI, method = "euclidean") -> dist.LPI.transf
vegdist(plot_pcnm_transf$prec, method = "euclidean") -> dist.prec.transf
vegdist(plot_pcnm_transf$WEI, method = "euclidean") -> dist.wei.transf

##dbRDA data----
beta.pair(wood_pa) -> wood.pair.pa
wood.pair.pa$beta.sor -> wood.pa.tot
wood.pair.pa$beta.sim -> wood.pa.tu
wood.pair.pa$beta.sne -> wood.pa.ne

## Wood ----
### turnover ----
decay.model(
  y = wood.pa.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.wood.prec


decay.model(
  y = wood.pa.tu,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.wood.lpi


decay.model(
  y = wood.pa.tu,
  x = dist.wei,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.wood.wei

decay.model(
  y = wood.pa.tu,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.wood.precL

decay.model(
  y = wood.pa.tu,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.wood.precW
