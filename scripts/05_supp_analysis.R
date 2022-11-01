# Fri Nov 12 18:53:03 2021 ------------------------------
#Script para análises suplementares

# library ----
library(dplyr)

# data ----
read_xlsx("data/herbáceas_pp.xlsx") -> herb
read_xlsx("data/lenhosas_pp.xlsx") -> len

## data table ---- Data came from script VarPart_xxxx
group<- c(rep.int("herb", times = 3),rep.int("len", times = 3))
frac<- c(rep(c("total", "turn", "nest"), times = 2))
LPI<- c(0.08, 0, 0, 0, 0, 0)
WEI<- c(0, 0, 0, 0.12, 0.12, 0)
prec<- c(0.05, 0.06, 0, 0.05, 0.05, 0)
space<- c(0.13, 0.11, 1, 0.08, 0.08, 1) 
unk<- c(0.01, 0, 0, 0.02, 0.02, 0)
perc_var<- data.frame(group, frac, LPI, WEI, prec, space, unk)
perc_var<- mutate(perc_var, total = LPI+ WEI + prec + space + unk)

perc_var%>%
  mutate(prop_LPI = LPI/total,
         prop_WEI = WEI/total,
         prop_prec = prec/total,
         prop_space = space/total
         )%>%
  glimpse

#Proportion of turnover and nestedness----
##herbs----
decostand(herb[,-c(1,59)], 'hellinger')-> herb_abund_hell
beta.pair.abund(herb_abund_hell) -> herb.pair.abund
herb.pair.abund$beta.bray -> herb.abund.tot
herb.pair.abund$beta.bray.bal -> herb.abund.tu
herb.pair.abund$beta.bray.gra -> herb.abund.ne

##wood----
decostand(len[-3,-1], 'hellinger') -> len_abund_hell
beta.pair.abund(len_abund_hell) ->len.pair.abund
len.pair.abund$beta.bray ->len.abund.tot
len.pair.abund$beta.bray.bal -> len.abund.tu
len.pair.abund$beta.bray.gra ->len.abund.ne

##analysis----
###herbs----
cbind(herb.abund.tot, herb.abund.tu, herb.abund.ne)  %>% 
  as_tibble() %>% 
  mutate(turnProp = herb.abund.tu/herb.abund.tot,
         nestProp = herb.abund.ne/herb.abund.tot,
         mean_turnProp = mean(turnProp),
         mean_nestProp = mean(nestProp),
         sd_turnProp = sd(turnProp),
         sd_nestProp = sd(nestProp)) %>% 
  glimpse

###woody----
cbind(len.abund.tot, len.abund.tu, len.abund.ne)  %>% 
  as_tibble() %>% 
  mutate(turnProp = len.abund.tu/len.abund.tot,
         nestProp = len.abund.ne/len.abund.tot,
         mean_turnProp = mean(turnProp),
         mean_nestProp = mean(nestProp),
         sd_turnProp = sd(turnProp),
         sd_nestProp = sd(nestProp)) %>% 
  glimpse
