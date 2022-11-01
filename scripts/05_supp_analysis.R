# Fri Nov 12 18:53:03 2021 ------------------------------
#Script para análises suplementares

# library ----
library(dplyr)

# data ----
perc_var%>%
  mutate(prop_LPI = LPI/total,
         prop_WEI = WEI/total,
         prop_prec = prec/total,
         prop_space = space/total
         )%>%
  glimpse
