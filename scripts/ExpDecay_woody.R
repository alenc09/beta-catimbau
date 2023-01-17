# Tue Nov  1 11:08:57 2022 ------------------------------
#Script for exponential decay analysis on woody communities

#libraries----
#data----
#Analysis----

## Wood ----
### turnover ----
decay.model(
  y = len.abund.tu,
  x = dist.prec,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.prec


decay.model(
  y = len.abund.tu,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.lpi


decay.model(
  y = len.abund.tu,
  x = dist.wei,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.wei

decay.model(
  y = len.abund.tu,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.precL

decay.model(
  y = len.abund.tu,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.precW

### nestedness ----

decay.model(
  y = len.abund.ne,
  x = dist.prec,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.prec.ne

decay.model(
  y = len.abund.ne,
  x = dist.LPI,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.lpi.ne

decay.model(
  y = len.abund.ne,
  x = dist.wei,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.wei.ne

decay.model(
  y = len.abund.ne,
  x = dist.prec.transf * dist.LPI.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.precL.ne

decay.model(
  y = len.abund.ne,
  x = dist.prec.transf * dist.wei.transf,
  model.type = "exponential",
  y.type = "dissimilarities"
) -> decay.abund.len.precW.ne

#Figure----
##pivot long the matrices----
### species dissimilarity matrices----
len.mat.tu<-as.matrix(len.abund.tu)
len.mat.tu[upper.tri(len.mat.tu,diag=T)]<-NA
len.melt.tu<- drop_na(as_tibble(melt(len.mat.tu)))

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