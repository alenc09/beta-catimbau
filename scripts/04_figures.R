# Fri Nov 12 15:54:15 2021 ------------------------------
#script para montar as figuras do artigo

#library ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# Variation partitioning ----
## data table ----
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
perc_var%>%
  pivot_longer(
    data = perc_var[,-8],
    cols = c(LPI, WEI, prec, space, unk),
    names_to = "source",
    values_to = "value_rsqr"
  )%>%
  mutate(across(.cols = 1:3, .fns = as.factor))%>%
  mutate(across(.cols = frac, ~factor(., levels = c("total", "turn", "nest"))))%>%
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
                      
ggsave("/home/lucas/Documentos/Doutorado/projetos_paralelos/artigo_a-b/manuscript/figures/fig1.jpg", plot = vp_bars)
