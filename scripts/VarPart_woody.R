# Tue Nov  1 11:04:05 2022 ------------------------------
#Script for variation partitioning analysis in woody plants dataset

#libraries----
#data----
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