# Fri Nov 12 15:54:15 2021 ------------------------------
#script para montar as figuras do artigo

#library ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(reshape2)

# Variation partitioning ----
## data table ---- #A tabela abaixo foi construida baseado nos resultados da partição da variância dos scripts VarPart_
group<- c(rep.int("herb", times = 3),rep.int("len", times = 3))
frac<- c(rep(c("total", "turn", "nest"), times = 2))
LPI<- c(0.08, 0, 0, 0, 0, 0)
WEI<- c(0, 0, 0, 0.12, 0.12, 0)
prec<- c(0.05, 0.06, 0, 0.05, 0.05, 0)
space<- c(0.13, 0.11, 1, 0.08, 0.08, 1) 
unk<- c(0.01, 0, 0, 0.02, 0.02, 0)
perc_var<- data.frame(group, frac, LPI, WEI, prec, space, unk)
perc_var<- mutate(perc_var, total = LPI+ WEI + prec + space + unk)

## figure ----
pivot_longer(
  data = perc_var[, -8],
  cols = c(LPI, WEI, prec, space, unk),
  names_to = "source",
  values_to = "value_rsqr"
) %>%
  mutate(across(.cols = 1:3, .fns = as.factor)) %>%
  mutate(across(.cols = frac, ~ factor(., levels = c(
    "total", "turn", "nest"
  )))) %>%
  glimpse -> perc_var2

facet_label<- c("total" = "Total Beta-diversity",
                "turn" = "Turnover",
                "nest" = "Nestedness")
facet_label<- as_labeller(facet_label)

ggplot(perc_var2, aes(x = group, y = value_rsqr, fill = source)) +
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8) +
  facet_grid(
    cols = vars(frac),
    switch = "x",
    labeller = as_labeller(facet_label)
  ) +
  scale_y_continuous(name = "% of each variable in R²adj",
                     labels = c("0", "25", "50", "75", "100")) +
  scale_x_discrete(name = "", labels = c("Herbaceous", "Wood")) +
  scale_fill_viridis_d(
    name = "",
    labels = c("LPI", "Precipitation", "Space", "Shared", "WEI")
  ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    panel.background = element_blank(),
    axis.line.y = element_line(colour = "grey")
  ) -> vp_bars
                      
#ggsave("/home/lucas/Documentos/Doutorado/projetos_paralelos/artigo_a-b/manuscript/figures/fig1.jpg", plot = vp_bars)

#Exponential decays----
##table for Figures----
### Variables matrices----
source(file = here("scripts/ExpDecay_herb_abund.R"))
source(file = here("scripts/ExpDecay_wood_abund.R"))
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

### species dissimilarity matrices----
herb.mat.tu<-as.matrix(herb.abund.tu)
herb.mat.tu[upper.tri(herb.mat.tu,diag=T)]<-NA
herb.melt.tu<- drop_na(as_tibble(melt(herb.mat.tu)))

wood.mat.tu<-as.matrix(wood.abund.tu)
wood.mat.tu[upper.tri(wood.mat.tu,diag=T)]<-NA
wood.melt.tu<- drop_na(as_tibble(melt(wood.mat.tu)))

###table exponential decay figures----
herb.melt.tu %>% 
  left_join(y = wood.melt.tu, by = c("Var1", "Var2")) %>% 
  rename("herb_tu" = "value.x" ,
         "wood_tu" = "value.y") %>% 
  left_join(y = prec.dist.melt, by = c("Var1", "Var2")) %>% 
  rename("dist_prec" = "value") %>% 
  left_join(y = lpi.dist.melt, by = c("Var1", "Var2")) %>% 
  rename("dist_lpi" = "value") %>%
  left_join(y = wei.dist.melt, by = c("Var1", "Var2")) %>% 
  rename("dist_wei" = "value") %>%
  left_join(y = precL.dist.melt, by = c("Var1", "Var2")) %>% 
  rename("dist_prec.lpi" = "value") %>%
  left_join(y = precW.dist.melt, by = c("Var1", "Var2")) %>% 
  rename("dist_prec.wei" = "value") %>%
  glimpse -> tab_dist_decay

##figure----
