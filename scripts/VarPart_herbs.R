# Tue Nov  1 11:01:36 2022 ------------------------------
#Script for variation partitioning analysis in herbaceus dataset

#libraries----
#data----
#Analysis----

## herbs ----
### total bet-diversity ----
anova(capscale(
  herb.abund.tot ~ 
    PCNM3 +
    PCNM1 + 
    PCNM10 + 
    Condition(LPI + prec), #only spatial effects
  data = plot_pcnm_transf
)) #Analysis to check the significance of each variable group in explaining the beta-diversity components
anova(capscale(
  herb.abund.tot ~ 
    LPI  + 
    Condition(PCNM3 + PCNM1 + PCNM10 + prec), #only disturbance effect
  data = plot_pcnm_transf
))
anova(capscale(
  herb.abund.tot ~ 
    prec + 
    Condition(PCNM3 + PCNM1 + PCNM10 + LPI), #only precipitation effect
  data = plot_pcnm_transf
))
varpart(herb_abund_hell,                    #Proper partitioning analysis
        plot_pcnm_transf[, 2],
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(8, 10, 17)]) -> varpart.abund.herb.tot #herb.abund.tot ~ LPI + prec + PCNM1 + PCNM3 + PCNM10
plot(varpart.abund.herb.tot)

### Turnover component ----
anova(capscale(herb.abund.tu ~ 
                 PCNM3 + 
                 PCNM1 + 
                 PCNM5 + 
                 Condition(prec),
               data = plot_pcnm_transf))
anova(capscale(herb.abund.tu ~ 
                 prec + 
                 Condition(PCNM3 + PCNM1 + PCNM5), data = plot_pcnm_transf))

varpart(herb_abund_hell,
        plot_pcnm_transf[, 6],
        plot_pcnm_transf[, c(8, 10, 12)]) -> varpart.abund.herb.tu #herb.abund.tu ~ prec + PCNM3 + PCNM1 + PCNM5
plot(varpart.abund.herb.tu)

### Nestedness component ----
#not necessary because it was only explained by spatial variables