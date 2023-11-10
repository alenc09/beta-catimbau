# Fri Nov 10 12:21:32 2023 ------------------------------

#libraries----
library(readxl)
library(dplyr)
library(tidyr)
library(geosphere)
library(vegan)
library(betapart)
library(reshape2)

#data----
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot
read_xlsx("data/herbáceas_pp.xlsx") -> herb

#Analysis----
##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm

##data transformation----
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf 
decostand(herb[,-c(1,59)], method = "pa") -> herb_pa
decostand(herb_pa, 'hellinger')-> herb_pa_hell

#Exponential decay models ----
##distance matrices ----
vegdist(plot_pcnm$LPI, method = "euclidean") -> dist.LPI
vegdist(plot_pcnm$prec, method = "euclidean") -> dist.prec
vegdist(plot_pcnm$WEI, method = "euclidean") -> dist.wei

vegdist(plot_pcnm_transf$LPI, method = "euclidean") -> dist.LPI.transf
vegdist(plot_pcnm_transf$prec, method = "euclidean") -> dist.prec.transf
vegdist(plot_pcnm_transf$WEI, method = "euclidean") -> dist.wei.transf

##beta-diversity----
herb.pair.pa<- beta.pair(herb_pa)
herb.pa.tot<- herb.pair.pa$beta.sor
herb.pa.tu<- herb.pair.pa$beta.sim
herb.pa.ne<- herb.pair.pa$beta.sne

## turnover ----
decay.model(
  y = herb.pa.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.herb.prec

decay.model(
  y = herb.pa.tu,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.herb.lpi

decay.model(
  y = herb.pa.tu,
  x = dist.wei,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.herb.wei

decay.model(
  y = herb.pa.tu,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.herb.precL

decay.model(
  y = herb.pa.tu,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "similarities"
) -> decay.pa.herb.precW
