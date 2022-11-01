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

