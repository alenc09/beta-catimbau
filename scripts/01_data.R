# Tue Oct 08 14:01:57 2019 ------------------------------
# Script para importar e organizar os dados necessários para as análises


#libraries####
library(readxl)
library(geosphere)
library(vegan)


#data####
read_xlsx("exp_pp_sem_P7.xlsx")-> plot
read_xlsx("herbáceas_pp.xlsx") -> herb
read_xlsx("lenhosas_pp.xlsx") -> len


#PCNMs####
mat_dist<- distm(plot[,c('lon','lat')], plot[,c('lon','lat')], fun=distVincentyEllipsoid)
pcnms<- pcnm(mat_dist)
plot_pcnm<- cbind(plot, pcnms$vectors)

#data transformation####
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(herb[,-c(1,59)], 'hellinger')-> herb_abund_hell
decostand(len[-3,-1], 'hellinger') -> len_abund_hell

##Herb - dbRDA data----
beta.pair.abund(herb_abund_hell) -> herb.pair.abund
herb.pair.abund$beta.bray -> herb.abund.tot
herb.pair.abund$beta.bray.bal -> herb.abund.tu
herb.pair.abund$beta.bray.gra -> herb.abund.ne

##Lenhosas - dbRDA data----
beta.pair.abund(len_abund_hell) ->len.pair.abund
len.pair.abund$beta.bray ->len.abund.tot
len.pair.abund$beta.bray.bal -> len.abund.tu
len.pair.abund$beta.bray.gra ->len.abund.ne


#Lenhosas



#varpart.abund.len.ne<- not possible to do a partioning for it is explained only by geographic variation (PCNM3)

#Analysis with presence/absence data############################################################################
#species data transformation####
herb.pa<- read.csv(file = "herb_pa.csv")
len.pa<- read.csv(file="len_pa.csv")
herb.pah<- decostand(herb.pa[,-1], method = "hellinger")
len.pah<- decostand(len.pa[,-1], method = "hellinger")



# Tue Nov 05 13:40:00 2019 ------------------------------


#P/A
decay.pa.herb.prec<- decay.model(y = herb.pa.tu, x = dist.prec, model.type = "exponential", y.type = "dissimilarities")
decay.pa.herb.lpi <- decay.model(y = herb.pa.tu, x = dist.LPI,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.herb.wei <- decay.model(y = herb.pa.tu, x = dist.wei,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.herb.gmdi<- decay.model(y = herb.pa.tu, x = dist.gmdi, model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.prec <- decay.model(y = len.pa.tu,  x = dist.prec, model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.lpi  <- decay.model(y = len.pa.tu,  x = dist.LPI,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.wei  <- decay.model(y = len.pa.tu,  x = dist.wei,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.gmdi <- decay.model(y = len.pa.tu,  x = dist.gmdi, model.type = "exponential", y.type = "dissimilarities")

#Decay models - interaction#####
#abundance
decay.abund.herb.precL<-decay.model(y = herb.abund.tu, x = dist.prec.transf*dist.LPI.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.precL<- decay.model(y = len.abund.tu,  x = dist.prec.transf*dist.LPI.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.herb.precW<-decay.model(y = herb.abund.tu, x = dist.prec.transf*dist.wei.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.precW<- decay.model(y = len.abund.tu,  x = dist.prec.transf*dist.wei.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.herb.precG< decay.model(y = herb.abund.tu, x = dist.prec.transf*dist.gmdi.transf, model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.precG<- decay.model(y = len.abund.tu,  x = dist.prec.transf*dist.gmdi.transf, model.type = "exponential", y.type = "dissimilarities")

#P/A
decay.pa.herb.precL<- decay.model(y = herb.pa.tu, x = dist.prec.transf*dist.LPI.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.precL<-  decay.model(y = len.pa.tu,  x = dist.prec.transf*dist.LPI.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.herb.precW<- decay.model(y = herb.pa.tu, x = dist.prec.transf*dist.wei.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.precW<-  decay.model(y = len.pa.tu,  x = dist.prec.transf*dist.wei.transf,  model.type = "exponential", y.type = "dissimilarities")
decay.pa.herb.precG<- decay.model(y = herb.pa.tu, x = dist.prec.transf*dist.gmdi.transf, model.type = "exponential", y.type = "dissimilarities")
decay.pa.len.precG<- decay.model(y = len.pa.tu, x = dist.prec.transf*dist.gmdi.transf, model.type = "exponential", y.type = "dissimilarities")

# Tue Nov 12 13:45:43 2019 ------------------------------
#bootstrap to compare slope and intercept####
set.seed(123)
#abundance
#precipitation
a<- boot.coefs.decay(decay.abund.herb.prec, 1000)
b<- boot.coefs.decay(decay.abund.len.prec, 1000)
sum(a$boot.coefs[,2] > b$boot.coefs[,2])/1000
sum(a$boot.coefs[,1] < b$boot.coefs[,1])/1000

#LPI
c<- boot.coefs.decay(decay.abund.herb.wei, 1000)
d<- boot.coefs.decay(decay.abund.len.wei, 1000)
sum(c$boot.coefs[,1] < d$boot.coefs[,1])/1000
sum(c$boot.coefs[,2] > d$boot.coefs[,2])/1000

#GMID
e<- boot.coefs.decay(decay.abund.herb.gmdi, 1000)
f<- boot.coefs.decay(decay.abund.len.gmdi, 1000)
sum(e$boot.coefs[,1] > f$boot.coefs[,1])/1000
sum(e$boot.coefs[,2] < f$boot.coefs[,2])/1000

#Prec*LPI
g<- boot.coefs.decay(decay.abund.herb.precL, 1000)
h<- boot.coefs.decay(decay.abund.len.precL, 1000)
sum(g$boot.coefs[,1] > h$boot.coefs[,1])/1000
sum(g$boot.coefs[,2] < h$boot.coefs[,2])/1000

#prec*WEI
i<- boot.coefs.decay(decay.abund.herb.precW, 1000)
j<- boot.coefs.decay(decay.abund.len.precW, 1000)
sum(i$boot.coefs[,1] < j$boot.coefs[,1])/1000
sum(i$boot.coefs[,2] > j$boot.coefs[,2])/1000

#prec*GMDI
k<- boot.coefs.decay(decay.abund.herb.precG, 1000)
l<- boot.coefs.decay(decay.abund.len.precG, 1000)
sum(k$boot.coefs[,1] < l$boot.coefs[,1])/1000
sum(k$boot.coefs[,2] > l$boot.coefs[,2])/1000

#P/A
##precipitation
m<- boot.coefs.decay(decay.pa.herb.prec, 1000)
n<- boot.coefs.decay(decay.pa.len.prec, 1000)
sum(m$boot.coefs[,1] > n$boot.coefs[,1])/1000
sum(m$boot.coefs[,2] > b$boot.coefs[,2])/1000

#wei
o<- boot.coefs.decay(decay.pa.herb.wei, 1000)
p<- boot.coefs.decay(decay.pa.len.wei, 1000)
sum(o$boot.coefs[,1] < p$boot.coefs[,1])/1000
sum(o$boot.coefs[,2] > p$boot.coefs[,2])/1000

#gmdi
q<- boot.coefs.decay(decay.pa.herb.gmdi, 1000)
r<- boot.coefs.decay(decay.pa.len.gmdi, 1000)
sum(q$boot.coefs[,1] > r$boot.coefs[,1])/1000
sum(q$boot.coefs[,2] > r$boot.coefs[,2])/1000

#prec*lpi
s<- boot.coefs.decay(decay.pa.herb.precL, 1000)
t<- boot.coefs.decay(decay.pa.len.precL, 1000)
sum(s$boot.coefs[,1] > t$boot.coefs[,1])/1000
sum(s$boot.coefs[,2] > t$boot.coefs[,2])/1000

#prec*wei
u<- boot.coefs.decay(decay.pa.herb.precW, 1000)
v<- boot.coefs.decay(decay.pa.len.precW, 1000)
sum(u$boot.coefs[,1] < v$boot.coefs[,1])/1000
sum(u$boot.coefs[,2] > v$boot.coefs[,2])/1000

#prec*gmdi
w<- boot.coefs.decay(decay.pa.herb.precG, 1000)
x<- boot.coefs.decay(decay.pa.len.precG, 1000)
sum(w$boot.coefs[,1] > x$boot.coefs[,1])/1000
sum(w$boot.coefs[,2] > x$boot.coefs[,2])/1000

# Mon Dec 02 16:14:13 2019 ------------------------------
        

# Thu Dec 05 16:16:07 2019 ------------------------------
# Decay model in ggplot2####
#melt the matrices
herb.mat.tu<-as.matrix(herb.abund.tu)
herb.mat.tu[upper.tri(herb.mat.tu,diag=T)]<-NA
herb.melt.tu<- drop_na(as_tibble(melt(herb.mat.tu)))

len.mat.tu<-as.matrix(len.abund.tu)
len.mat.tu[upper.tri(len.mat.tu,diag=T)]<-NA
len.melt.tu<- drop_na(as_tibble(melt(len.mat.tu)))

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


herb.melt.tu <- dplyr::rename(herb.melt.tu, herb.tu = "value")
len.melt.tu <- dplyr::rename(len.melt.tu, len.tu = "value")
prec.dist.melt <- dplyr::rename(prec.dist.melt, prec = "value")
lpi.dist.melt <- dplyr::rename(lpi.dist.melt, lpi = "value")
wei.dist.melt <- dplyr::rename(wei.dist.melt, wei = "value")
precL.dist.melt <- dplyr::rename(precL.dist.melt, precL = "value")
precW.dist.melt<- dplyr::rename(precW.dist.melt, precW = "value")

tab_melt<- cbind(herb.melt.tu, len.melt.tu[,c(-1,-2)], prec.dist.melt[,c(-1,-2)],
                 lpi.dist.melt[,c(-1,-2)], wei.dist.melt[,c(-1,-2)], 
                 precL.dist.melt[,c(-1,-2)], precW.dist.melt[,c(-1,-2)])
head(tab_melt)

#decay graphics and panel
ggplot(tab_melt, mapping = aes(x = prec, y= 1- herb.tu))+
  geom_point(shape = 1, color = "black")+
stat_function(fun=function(x){coef(decay.abund.herb.prec)["a.intercept"][["(Intercept)"]] * coef(decay.abund.herb.prec)["b.slope"][["x"]] ** x})
  #stat_smooth(mapping = aes(x = prec, y= 1 - herb.tu), formula = y ~ a*b^x,
  #            fun = "nls",
   #           method = "lm")
#  geom_point(mapping = aes(x=prec, y= 1 -len.tu), shape = 1, color = "grey60")
  