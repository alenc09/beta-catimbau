# Tue Nov  1 11:07:45 2022 ------------------------------
#Script for exponential decay analysis on herbaceous communities

#libraries----
#data----
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