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