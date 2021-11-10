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





ordistep(capscale(len.abund.ne~1, plot_pcnm_transf), scope = formula(modT2.abund.len),
         direction = 'both', pstep=1000)
mod.len.abund.ne<- capscale(len.abund.ne ~ PCNM3, plot_pcnm_transf)
anova(mod.len.abund.ne)
RsquareAdj(mod.len.abund.ne)


#variation partitioning - Abundance####
#herb?ceas

varpart.abund.herb.rda<- varpart(herb_abund_hell, plot_pcnm_transf[,6], plot_pcnm_transf[, 10]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.abund.herb.rda)

varpart.pa.herb.tot<- varpart(herb_abund_hell, plot_pcnm_transf[,6], plot_pcnm_transf[,c(8, 10, 12)]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.herb.tot)

anova(capscale(herb.abund.tot ~ PCNM3 + PCNM1 + PCNM10 + Condition(LPI + prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.tot ~ LPI  + Condition(PCNM3 + PCNM1 + PCNM10+ prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.tot ~ prec + Condition(PCNM3 + PCNM1 + PCNM10+ LPI), data = plot_pcnm_transf))
varpart.abund.herb.tot<- varpart(herb_abund_hell, plot_pcnm_transf[,2], plot_pcnm_transf[,6], plot_pcnm_transf[,c(8, 10, 17)]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.abund.herb.tot)

anova(capscale(herb.abund.tu ~ PCNM3 + PCNM1 + Condition(LPI + prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.tu ~ LPI  + Condition(PCNM3 + PCNM1+ prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.tu ~  prec + Condition(PCNM3 + PCNM1+ LPI), data = plot_pcnm_transf))
varpart.abund.herb.tu<- varpart(herb_abund_hell, plot_pcnm_transf[,2], plot_pcnm_transf[,6], plot_pcnm_transf[,c(8, 10)]) #herb.abund.tu ~ PCNM3 + LPI + prec + PCNM1
plot(varpart.abund.herb.tu)

anova(capscale(herb.abund.ne ~ PCNM7 + PCNM3 + Condition(LPI + prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.ne ~ LPI + Condition(PCNM7 + PCNM3+ prec), data = plot_pcnm_transf))
anova(capscale(herb.abund.ne ~ prec + Condition(LPI + PCNM7 + PCNM3), data = plot_pcnm_transf))
varpart.abund.herb.ne<- varpart(herb_abund_hell, plot_pcnm_transf[,2], plot_pcnm_transf[,6], plot_pcnm_transf[,c(14, 10)]) #herb.abund.tu ~ PCNM3 + LPI + prec + PCNM1
plot(varpart.abund.herb.ne)

#varpart.abund.herb.ne<- not possible to do a partioning for it is explained only by geographic variation (PCNM7 + PCNM3)

#Lenhosas

varpart.pa.len.rda<- varpart(len.pah, plot_pcnm_transf[,c(3,6)], plot_pcnm_transf[, 8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.len.rda)
varpart.abund.len.rda<- varpart(len_abund_hell, plot_pcnm_transf[,6], plot_pcnm_transf[, 8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.abund.len.rda)

varpart.pa.len.tot<- varpart(len.pa.tot, plot_pcnm_transf[,c(6,3)], plot_pcnm_transf[, 8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.len.tot)
varpart.pa.len.tu<- varpart(len.pa.tu, plot_pcnm_transf[,c(6,3)], plot_pcnm_transf[, 8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.len.tu)

anova(capscale(len.abund.tot ~ Condition(prec) + PCNM1, data = plot_pcnm_transf))
anova(capscale(len.abund.tot ~ prec + Condition(PCNM1), data = plot_pcnm_transf))
varpart.abund.len.tot<- varpart(len_abund_hell, plot_pcnm_transf[,6], plot_pcnm_transf[,8]) #len.abund.tot ~ prec + PCNM1
plot(varpart.abund.len.tot)

anova(capscale(len.abund.tu ~ Condition(WEI + prec) + PCNM1 + PCNM5, data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~ WEI  + Condition(PCNM1 + PCNM5+ prec), data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~  prec + Condition(PCNM1 + PCNM5+ WEI), data = plot_pcnm_transf))
varpart.abund.len.tu<- varpart(len_abund_hell,plot_pcnm_transf[,3] ,plot_pcnm_transf[,6], plot_pcnm_transf[,c(8, 12)]) #len.abund.tu ~ WEI + prec + PCNM1 + PCNM5
plot(varpart.abund.len.tu)

anova(capscale(len.abund.ne ~ Condition(WEI + prec) + PCNM3, data = plot_pcnm_transf))
anova(capscale(len.abund.ne ~ WEI  + Condition(PCNM3+ prec), data = plot_pcnm_transf))
anova(capscale(len.abund.ne ~  prec + Condition(PCNM3 + WEI), data = plot_pcnm_transf))
varpart.abund.len.ne<- varpart(len_abund_hell,plot_pcnm_transf[,3] ,plot_pcnm_transf[,6], plot_pcnm_transf[,10]) #len.abund.tu ~ WEI + prec + PCNM1 + PCNM5
plot(varpart.abund.len.ne)

#varpart.abund.len.ne<- not possible to do a partioning for it is explained only by geographic variation (PCNM3)

#Analysis with presence/absence data############################################################################
#species data transformation####
herb.pa<- read.csv(file = "herb_pa.csv")
len.pa<- read.csv(file="len_pa.csv")
herb.pah<- decostand(herb.pa[,-1], method = "hellinger")
len.pah<- decostand(len.pa[,-1], method = "hellinger")

#RDA with PA data####
#herbaceas - RDA
modT_herb_pa<- rda(herb.pah ~ WEI + LPI + prec + fert_sol + PCNM1 + PCNM3 +
                        PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11,
                      data = plot_pcnm_transf)
vif.cca(modT_herb_pa)
modselec_herb_pa<- ordistep(cca(herb.pah~1, plot_pcnm_transf), scope = formula(modT_herb_pa),
                         direction = 'both', pstep=1000)
mod_herb_pa<- rda(herb.pah ~ PCNM5 + PCNM1, plot_pcnm_transf)
anova(mod_herb_pa)
RsquareAdj(mod_herb_pa)

#Herb?ceas - dbRDA
herb.pair.pa<- beta.pair(herb.pa[,-1])
herb.pa.tot<- herb.pair.pa$beta.sor
herb.pa.tu<- herb.pair.pa$beta.sim
herb.pa.ne<- herb.pair.pa$beta.sne
set.seed(123)
ordistep(capscale(herb.pa.tot~1, plot_pcnm_transf), scope = formula(modT_herb_pa),
         direction = 'both', pstep=1000)
mod.herb.pa.tot<- capscale(herb.pa.tot ~ PCNM1 + PCNM3 + PCNM5 + prec, plot_pcnm_transf)
anova(mod.herb.pa.tot)
RsquareAdj(mod.herb.pa.tot)
plot(mod.herb.pa.tot)

ordistep(capscale(herb.pa.tu~1, plot_pcnm_transf), scope = formula(modT_herb_pa),
         direction = 'both', pstep=1000)
mod.herb.pa.tu<- capscale(herb.pa.tu ~ 1, plot_pcnm_transf)
anova(mod.herb.pa.tu)
RsquareAdj(mod.herb.pa.tu)
plot(mod.herb.abund.tu)

ordistep(capscale(herb.pa.ne~1, plot_pcnm_transf), scope = formula(modT_herb_pa),
         direction = 'both', pstep=1000)
mod.herb.pa.ne<- capscale(herb.pa.ne ~ PCNM3 + PCNM7 + PCNM8, plot_pcnm_transf)
anova(mod.herb.pa.ne)
RsquareAdj(mod.herb.pa.ne)
plot(mod.herb.pa.ne)

#lenhosas - RDA
modT_len_pa<- rda(len.pah ~ LPI + WEI + prec + fert_sol + PCNM1 + PCNM3 +
                     PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + PCNM9 + PCNM10 + PCNM11,
                   data = plot_pcnm_transf)
vif.cca(modT_len_pa)
modselec_len_pa<- ordistep(cca(len.pah~1, plot_pcnm_transf), scope = formula(modT_len_pa),
                            direction = 'both', pstep=1000)
mod_len_pa<- rda(len.pah ~ prec + PCNM1 + WEI, plot_pcnm_transf)
anova(mod_len_pa)
RsquareAdj(mod_len_pa)

#Lenhosas - dbRDA
len.pair.pa<- beta.pair(len.pa[,-1])
len.pa.tot<- len.pair.pa$beta.sor
len.pa.tu<- len.pair.pa$beta.sim
len.pa.ne<- len.pair.pa$beta.sne
set.seed(123)
ordistep(capscale(len.pa.tot~1, plot_pcnm_transf), scope = formula(modT_len_pa),
         direction = 'both', pstep=1000)
mod.len.pa.tot<- capscale(len.pa.tot ~ prec + PCNM1 + WEI, plot_pcnm_transf)
anova(mod.len.pa.tot)
RsquareAdj(mod.len.pa.tot)
plot(mod.len.pa.tot)

ordistep(capscale(len.pa.tu~1, plot_pcnm_transf), scope = formula(modT_len_pa),
         direction = 'both', pstep=1000)
mod.len.pa.tu<- capscale(len.pa.tu ~ WEI + prec + PCNM1 + LPI, plot_pcnm_transf)
anova(mod.len.pa.tu)
RsquareAdj(mod.len.pa.tu)
plot(mod.len.pa.tu)

ordistep(capscale(len.pa.ne~1, plot_pcnm_transf), scope = formula(modT_len_pa),
         direction = 'both', pstep=1000)
mod.len.pa.ne<- capscale(len.pa.ne ~ PCNM11, plot_pcnm_transf)
anova(mod.len.pa.ne)
RsquareAdj(mod.len.pa.ne) 
plot(mod.len.pa.ne)

# Tue Oct 29 11:18:22 2019 ------------------------------
#Variation Partitioning - presence/absence####
#herb?ceas
anova(capscale(herb.pa.tot ~ PCNM1 + PCNM3 + PCNM5 + Condition(prec), data = plot_pcnm_transf))
anova(capscale(herb.pa.tot ~ Condition(PCNM1 + PCNM3 + PCNM5) + prec, data = plot_pcnm_transf))
varpart.pa.herb.tot<- varpart(herb.pah, plot_pcnm_transf[,6], plot_pcnm_transf[,c(8, 10, 12)]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.herb.tot)


#Lenhosas
anova(capscale(len.pa.tot ~ PCNM1 + Condition(prec + WEI), data = plot_pcnm_transf))
anova(capscale(len.pa.tot ~ Condition(PCNM1) + prec + WEI, data = plot_pcnm_transf))
varpart.pa.len.tot<- varpart(len.pah, plot_pcnm_transf[,c(3,6)], plot_pcnm_transf[,8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.len.tot)

anova(capscale(len.pa.tu ~ PCNM1 + Condition(prec + WEI + LPI), data = plot_pcnm_transf))
anova(capscale(len.pa.tu ~ Condition(PCNM1) + prec + WEI + LPI, data = plot_pcnm_transf))
varpart.pa.len.tu<- varpart(len.pah, plot_pcnm_transf[,c(2,3,6)], plot_pcnm_transf[,8]) #herb.abund.tot ~ PCNM3 + LPI + prec + PCNM1 + PCNM10
plot(varpart.pa.len.tu)

# Tue Nov 05 13:40:00 2019 ------------------------------
#Distance matrixes####
dist.LPI<- vegdist(plot_pcnm$LPI, method = "euclidean")
dist.prec<- vegdist(plot_pcnm$prec, method = "euclidean")
dist.wei<- vegdist(plot_pcnm$WEI, method = "euclidean")
dist.gmdi<- vegdist(plot_pcnm$GMDI, method = "euclidean")

dist.LPI.transf<- vegdist(plot_pcnm_transf$LPI, method = "euclidean")
dist.prec.transf<- vegdist(plot_pcnm_transf$prec, method = "euclidean")
dist.wei.transf<- vegdist(plot_pcnm_transf$WEI, method = "euclidean")
dist.gmdi.transf<- vegdist(plot_pcnm_transf$GMDI, method = "euclidean")

#Decay models####
#abundance
decay <- nls()
decay.abund.herb.prec<- decay.model(y = herb.abund.tu, x = dist.prec, model.type = "exponential", y.type = "dissimilarities")
decay.abund.herb.lpi <- decay.model(y = herb.abund.tu, x = dist.LPI,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.herb.wei <- decay.model(y = herb.abund.tu, x = dist.wei,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.herb.gmdi<- decay.model(y = herb.abund.tu, x = dist.gmdi, model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.prec <- decay.model(y = len.abund.tu,  x = dist.prec, model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.lpi  <- decay.model(y = len.abund.tu,  x = dist.LPI,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.wei  <- decay.model(y = len.abund.tu,  x = dist.wei,  model.type = "exponential", y.type = "dissimilarities")
decay.abund.len.gmdi <- decay.model(y = len.abund.tu,  x = dist.gmdi, model.type = "exponential", y.type = "dissimilarities")

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
#table of percentage of explaned variance by each variable
group<- c(rep.int("herb", times = 3),rep.int("len", times = 3))
frac<- c(rep(c("total", "turn", "nest"), times = 2))
disturb<- c(0.08, 0.07, 0, 0, 0.06, 0)
prec<- c(0.05, 0.05, 0, 0.05, 0.06, 0)
space<- c(0.13, 0.09, 1, 0.08, 0.06, 1) 
perc_var<- data.frame(group, frac, disturb, prec, space)
perc_var<- mutate(perc_var, total = disturb + prec + space)

#barplot of percentage of explaned variance by each variable
perc_var2<- pivot_longer(data = perc_var[,-6], cols = c(disturb, prec, space),  names_to= "value",
             names_repair = "unique")
perc_var2$frac<- factor(perc_var2$frac, ordered = T, levels = c("total", "turn", "nest"))
facet_label<- c("total" = "Total Beta-diversity",
                  "turn" = "Turnover",
                  "nest" = "Nestedness")
facet_label<- as_labeller(facet_label)

vp_bars<- ggplot(perc_var2, aes(x=group, y=value...4, fill = value...3))+
  geom_bar(stat = "identity", position = "fill", width = 0.8)+
  facet_grid(cols = vars(frac), switch = "x", labeller = as_labeller(facet_label))+
  scale_y_continuous(name = "% of each variable in R?adj", labels = c("0","25", "50", "75", "100"))+
  scale_x_discrete(name = "", labels = c("Herbaceous", "Wood"))+
  scale_fill_hue(name = "" , labels = c("Disturbance", "Precipitation", "Space"))+
  theme(panel.spacing.x = unit(0, "lines"), panel.background = element_blank(),
        axis.line.y = element_line(colour = "grey"))
print(vp_bars)
ggsave("vp_bars.jpg", plot = vp_bars)        

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
  