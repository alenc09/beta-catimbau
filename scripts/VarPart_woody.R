# Tue Nov  1 11:04:05 2022 ------------------------------
#Script for variation partitioning analysis in woody plants dataset

#libraries----
#data----
read_xlsx("data/lenhosas_pp.xlsx") -> len
read_xlsx("data/exp_pp_sem_P7.xlsx")-> plot

#data transformation####
decostand(plot_pcnm[, -c(1:3)], 'standardize') -> plot_pcnm_transf #Faz sentido padronizar as PCNM's? Elas não perderiam a representação da variação espacial dessa forma?
decostand(len[-3,-1], 'hellinger') -> len_abund_hell

##variation part data----
beta.pair.abund(len_abund_hell) ->len.pair.abund
len.pair.abund$beta.bray ->len.abund.tot
len.pair.abund$beta.bray.bal -> len.abund.tu
len.pair.abund$beta.bray.gra ->len.abund.ne

#Analysis----
## Wood ----
### Total beta-diversity ----
anova(capscale(len.abund.tot ~ 
                 PCNM2 +
                 PCNM5 +
                 Condition(WEI + prec),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tot ~ 
                 WEI + 
                 Condition(prec + PCNM2 + PCNM5),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tot ~ 
                 prec + 
                 Condition(WEI + PCNM2 + PCNM5),
               data = plot_pcnm_transf))

varpart(len_abund_hell, 
        plot_pcnm_transf[,3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(9,12)]) -> varpart.abund.len.tot #len.abund.tot ~ WEI + prec + PCNM2 + PCNM5
plot(varpart.abund.len.tot)

### Turnover component ----
anova(capscale(len.abund.tu ~ 
                 PCNM2 +
                 PCNM5 +
                 Condition(WEI + prec),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~ 
                 WEI + 
                 Condition(prec + PCNM2 + PCNM5),
               data = plot_pcnm_transf))
anova(capscale(len.abund.tu ~ 
                 prec + 
                 Condition(WEI + PCNM2 + PCNM5),
               data = plot_pcnm_transf))

varpart(len_abund_hell, 
        plot_pcnm_transf[,3],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(9,12)]) -> varpart.abund.len.tu #len.abund.tot ~ WEI + prec + PCNM2 + PCNM5
plot(varpart.abund.len.tu)

### Nestedness component ----
#not necessary because it was only explained by spatial variables