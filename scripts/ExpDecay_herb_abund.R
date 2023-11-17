# Tue Nov  1 11:07:45 2022 ------------------------------
#Script for exponential decay analysis on herbaceous communities

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
decostand(herb[,-c(1,59)], 'hellinger')-> herb_abund_hell

#Exponential decay models ----
##distance matrices ----
vegdist(plot_pcnm$LPI, method = "euclidean") -> dist.LPI
vegdist(plot_pcnm$prec, method = "euclidean") -> dist.prec
vegdist(plot_pcnm$WEI, method = "euclidean") -> dist.wei

vegdist(plot_pcnm_transf$LPI, method = "euclidean") -> dist.LPI.transf
vegdist(plot_pcnm_transf$prec, method = "euclidean") -> dist.prec.transf
vegdist(plot_pcnm_transf$WEI, method = "euclidean") -> dist.wei.transf
for(i in plot_pcnm_transf[,8:18]){
  vegdist(i, method = "euclidean")
} -> list

## beta-diversity----
beta.pair.abund(herb_abund_hell) -> herb.pair.abund
herb.pair.abund$beta.bray -> herb.abund.tot
herb.pair.abund$beta.bray.bal -> herb.abund.tu
herb.pair.abund$beta.bray.gra -> herb.abund.ne

## turnover ----
decay.model(
  y = herb.abund.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.prec

plot(decay.abund.herb.prec)
plot(decay.abund.wood.prec, add = T, pty = 2)

decay.model(
  y = herb.abund.tu,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.lpi

decay.model(
  y = herb.abund.tu,
  x = dist.wei,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.wei

decay.model(
  y = herb.abund.tu,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.precL

decay.model(
  y = herb.abund.tu,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.precW

### nestedness ----
decay.model(
  y = herb.abund.ne,
  x = dist.prec,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.prec.ne

decay.model(
  y = herb.abund.ne,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.lpi.ne

decay.model(
  y = herb.abund.ne,
  x = dist.wei,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.wei.ne

decay.model(
  y = herb.abund.ne,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.precL.ne

decay.model(
  y = herb.abund.ne,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.precW.ne

#Figures----

### Variables matrices----
prec.dist.mat<- as.matrix(dist.prec) 
prec.dist.mat[upper.tri(prec.dist.mat, diag=T)]<-NA
prec.dist.melt<- drop_na(as_tibble(melt(prec.dist.mat)))

lpi.dist.mat<- as.matrix(dist.LPI)
lpi.dist.mat[upper.tri(lpi.dist.mat, diag=T)]<-NA
lpi.dist.melt<- drop_na(as_tibble(melt(lpi.dist.mat)))

wei.dist.mat<- as.matrix(dist.wei)
wei.dist.mat[upper.tri(wei.dist.mat, diag=T)]<-NA
wei.dist.melt<- drop_na(as_tibble(melt(wei.dist.mat)))

precL.dist.mat<- dist.prec.transf*dist.LPI.transf
precL.dist.mat<- as.matrix(precL.dist.mat)
precL.dist.mat[upper.tri(precL.dist.mat, diag=T)]<-NA
precL.dist.melt<- drop_na(as_tibble(melt(precL.dist.mat)))

precW.dist.mat<- dist.prec.transf*dist.wei.transf
precW.dist.mat<- as.matrix(precW.dist.mat)
precW.dist.mat[upper.tri(precW.dist.mat, diag=T)]<-NA
precW.dist.melt<- drop_na(as_tibble(melt(precW.dist.mat)))
