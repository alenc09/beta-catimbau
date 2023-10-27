# Tue Nov  1 11:07:45 2022 ------------------------------
#Script for exponential decay analysis on herbaceous communities

#libraries----
#data----

##PCNMs----
distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid) -> mat_dist
pcnm(mat_dist) -> pcnms
cbind(plot, pcnms$vectors) -> plot_pcnm
#Analysis----

#Exponential decay models ----
## data for distance matrixes ----
dist.LPI<- vegdist(plot_pcnm$LPI, method = "euclidean")
dist.prec<- vegdist(plot_pcnm$prec, method = "euclidean")
dist.wei<- vegdist(plot_pcnm$WEI, method = "euclidean")
dist.gmdi<- vegdist(plot_pcnm$GMDI, method = "euclidean")

dist.LPI.transf<- vegdist(plot_pcnm_transf$LPI, method = "euclidean")
dist.prec.transf<- vegdist(plot_pcnm_transf$prec, method = "euclidean")
dist.wei.transf<- vegdist(plot_pcnm_transf$WEI, method = "euclidean")
dist.gmdi.transf<- vegdist(plot_pcnm_transf$GMDI, method = "euclidean")

## Herbs ----
### turnover ----
decay.model(
  y = herb.abund.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.herb.prec

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
##pivot long the matrices----
### species dissimilarity matrices----
herb.mat.tu<-as.matrix(herb.abund.tu)
herb.mat.tu[upper.tri(herb.mat.tu,diag=T)]<-NA
herb.melt.tu<- drop_na(as_tibble(melt(herb.mat.tu)))



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