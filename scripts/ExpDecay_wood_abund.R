# Tue Nov  1 11:08:57 2022 ------------------------------
#Script for exponential decay analysis on woody communities

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
decostand(wood[-3,-1], 'hellinger') -> wood_abund_hell

#Exponential decay models ----
##distance matrices ----
vegdist(plot_pcnm$LPI, method = "euclidean") -> dist.LPI
vegdist(plot_pcnm$prec, method = "euclidean") -> dist.prec
vegdist(plot_pcnm$WEI, method = "euclidean") -> dist.wei

vegdist(plot_pcnm_transf$LPI, method = "euclidean") -> dist.LPI.transf
vegdist(plot_pcnm_transf$prec, method = "euclidean") -> dist.prec.transf
vegdist(plot_pcnm_transf$WEI, method = "euclidean") -> dist.wei.transf

##dbRDA data----
beta.pair.abund(wood_abund_hell) ->wood.pair.abund
wood.pair.abund$beta.bray ->wood.abund.tot
wood.pair.abund$beta.bray.bal -> wood.abund.tu
wood.pair.abund$beta.bray.gra ->wood.abund.ne

## Wood ----
### turnover ----
decay.model(
  y = wood.abund.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.prec


decay.model(
  y = wood.abund.tu,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.lpi


decay.model(
  y = wood.abund.tu,
  x = dist.wei,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.wei

decay.model(
  y = wood.abund.tu,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.precL

decay.model(
  y = wood.abund.tu,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.precW

### nestedness ----
decay.model(
  y = wood.abund.ne,
  x = dist.prec,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.prec.ne

decay.model(
  y = wood.abund.ne,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.lpi.ne

decay.model(
  y = wood.abund.ne,
  x = dist.wei,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.wei.ne

decay.model(
  y = wood.abund.ne,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.precL.ne

decay.model(
  y = wood.abund.ne,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.abund.wood.precW.ne
