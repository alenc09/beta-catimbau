# Fri Nov 12 15:54:15 2021 ------------------------------
#script para montar as figuras do artigo

#library ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(reshape2)
library(RColorBrewer)
library(here)
library(cowplot)
library(ggpubr)

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
  scale_fill_manual(
    values = c("#738054","#a3a86d", "grey40","#be996e", "#a98e73"),
    name = "",
    labels = c("Precipitation", "Space", "Shared", "LPI", "WEI")
  ) +
  theme(
    panel.spacing.x = unit(0, "lines"),
    panel.background = element_blank(),
    axis.line.y = element_line(colour = "grey80")
  ) -> vp_bars
                      
ggsave(here("figures/varpart.jpg"), plot = vp_bars)

#Exponential decays----
##table for Figures----
### Variables matrices----
source(file = here("scripts/ExpDecay_herb_abund.R"))
source(file = here("scripts/ExpDecay_wood_abund.R"))

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

##figures----
###precipitation distance----
ggplot(data = tab_dist_decay) +
  geom_point(aes(x = dist_prec, y = herb_tu - 1, color = "herb_tu"),
             alpha = 0.4)+
  geom_smooth(aes(x = dist_prec, y = herb_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#738054") +
  geom_point(aes(x = dist_prec, y = wood_tu - 1, color = "wood_tu"),
             alpha = 0.4) +
  geom_smooth(aes(x = dist_prec, y = wood_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#aa722a")+
  scale_y_continuous(breaks = c(-1, -0.8, -0.6,-0.4,-0.2,0), labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+ #modifiquei o label so para o valor de turnover ficar positivo no grafico como ele realmente é
  labs(x = "Precipitation", y = "Abundance-based turnover")+
  scale_color_manual(name = "", 
                     labels = c("Herbaceous", "Wood"),
                     values = c("herb_tu" = "#738054", "wood_tu" = "#aa722a"))+
  theme_classic() -> prec_decay_fig

###LPI distance----
ggplot(data = tab_dist_decay) +
  geom_point(aes(x = dist_lpi, y = herb_tu - 1, color = "herb_tu"),
             alpha = 0.4)+
  geom_smooth(aes(x = dist_lpi, y = herb_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#738054") +
  geom_point(aes(x = dist_lpi, y = wood_tu - 1, color = "wood_tu"),
             alpha = 0.4) +
  geom_smooth(aes(x = dist_lpi, y = wood_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#aa722a")+
  labs(x = "Livestock Pressure Index", y = "")+
  scale_y_continuous(breaks = c(-1, -0.8, -0.6,-0.4,-0.2,0), labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_color_manual(name = "", 
                     labels = c("Herbaceous", "Wood"),
                     values = c("herb_tu" = "#738054", "wood_tu" = "#aa722a"))+
  theme_classic() -> lpi_decay_fig

###WEI distance----
ggplot(data = tab_dist_decay) +
  geom_point(aes(x = dist_wei, y = herb_tu - 1, color = "herb_tu"),
             alpha = 0.4)+
  geom_smooth(aes(x = dist_wei, y = herb_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#738054") +
  geom_point(aes(x = dist_wei, y = wood_tu - 1, color = "wood_tu"),
             alpha = 0.4) +
  geom_smooth(aes(x = dist_wei, y = wood_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "solid",
              lwd = 1,
              color = "#aa722a")+
  labs(x = "Wood Extraction Index", y = "")+
  scale_y_continuous(breaks = c(-1, -0.8, -0.6,-0.4,-0.2,0), labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_color_manual(name = "", 
                     labels = c("Herbaceous", "Wood"),
                     values = c("herb_tu" = "#738054", "wood_tu" = "#aa722a"))+
  theme_classic() -> wei_decay_fig

###prec*lpi distance----
ggplot(data = tab_dist_decay) +
  geom_point(aes(x = dist_prec.lpi, y = herb_tu - 1, color = "herb_tu"),
             alpha = 0.4)+
  geom_smooth(aes(x = dist_prec.lpi, y = herb_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#738054") +
  geom_point(aes(x = dist_prec.lpi, y = wood_tu - 1, color = "wood_tu"),
             alpha = 0.4) +
  geom_smooth(aes(x = dist_prec.lpi, y = wood_tu),
              method = "nls",
              formula = y - 1 ~ a*exp(b*x),
              method.args = list(start = c(a=1, b=0)),
              se = F,
              linetype = "dashed",
              lwd = 1,
              color = "#aa722a")+
  scale_y_continuous(breaks = c(-1, -0.8, -0.6,-0.4,-0.2,0), labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+
  labs(x = "Precipitation and LPI interaction", y = "Abundance-based turnover")+
  scale_color_manual(name = "", 
                     labels = c("Herbaceous", "Wood"),
                     values = c("herb_tu" = "#738054", "wood_tu" = "#aa722a"))+
  theme_classic() -> prec.lpi_decay_fig

###prec*wei distance----
ggplot(data = tab_dist_decay) +
  geom_point(aes(x = dist_prec.wei, y = herb_tu -1, color = "herb_tu"), #reduzi -1 da matriz para os pontos ficarem no mesmo lugar que a curva
           alpha = 0.4)+
  geom_smooth(aes(x = dist_prec.wei, y = herb_tu),
            method = "nls", #Esse é o método utilizado para plotar uma curva exponencial
            formula = y - 1 ~ a*exp(b*x), #reduzi -1 do y no modelo para a curva ficar em funcao da dissimilaridade. Foi o único jeito que encontrei de inverter a direcao da curva
            method.args = list(start = c(a=1, b=0)), #especificacoes do modelo
            se = F, #aparentemente nao da para plotar isso com o standar error do modelo
            linetype = "dashed", #a linha tracejada é nao significativo
            lwd = 1,
            color = "#738054") +
geom_point(aes(x = dist_prec.wei, y = wood_tu -1 , color = "wood_tu"), 
           alpha = 0.4) +
  geom_smooth(aes(x = dist_prec.wei, y = wood_tu),
              method = "nls",
              formula = y - 1 ~ a * exp(b * x), 
              method.args = list(start = c(a = 1, b = 0)),
              se = F,
              linetype = "solid", #linha solida é significativo
              lwd = 1,
              color = "#aa722a")+
  labs(x = "Precipitation and WEI interaction", y = "")+
  scale_color_manual(name = "",
                     labels = c("Herbaceous", "Wood"),
                     values = c("herb_tu" = "#738054", "wood_tu" = "#aa722a"))+
  scale_y_continuous(breaks = c(-1, -0.8, -0.6,-0.4,-0.2,0), labels = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))+ #modifiquei o label so para o valor de turnover ficar positivo no grafico como ele realmente eh
  theme_classic()+
  theme(legend.text = element_text(size = 12)) -> prec.wei_decay_fig

get_legend(prec.wei_decay_fig) %>% 
  as_ggplot() -> lengend_panel

###Panel----
ggarrange(prec_decay_fig, lpi_decay_fig, wei_decay_fig, prec.lpi_decay_fig, prec.wei_decay_fig, lengend_panel,
          labels = c("a", "b", "c", "d", "e", ""),
          common.legend = T,
          legend = "none"
          ) -> fig.decay

ggsave(filename = here("figures/fig.decay.jpg"),
       plot = fig.decay,
       dpi = 300,
       width = 11.5)

#turnover profile----
c(seq(1:19)) -> number
c("P02", "P04", "P08", "P10", "P11", "P14", "P15", "P16", "P17", "P20", "P21", "P22", "P23", "P25", "P26", "P27", "P28", "P29", "P30") -> plot
condicoes2 <- purrr::map2(number, plot, ~quo(Var2 == !!.x ~ !!.y))
condicoes2 <- c(condicoes2, quo(TRUE ~ as.character(Var2)))
condicoes1 <- purrr::map2(number, plot, ~quo(Var1 == !!.x ~ !!.y))
condicoes1 <- c(condicoes1, quo(TRUE ~ as.character(Var1)))

tab_dist_decay %>% 
  mutate(Var2 = case_when(!!!condicoes2),
         Var1 = case_when(!!!condicoes1)) %>%
  glimpse -> tab_dist_plot

