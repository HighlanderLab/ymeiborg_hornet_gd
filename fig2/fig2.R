library(dplyr)
library(patchwork)

setwd("/scratch/bell/ymeiborg/fig2")

load(file = "../fig2a/Fig2a_plot.Rdata")
load(file = "../fig2b/Fig2b_plot.Rdata")

p <- fig2a / fig2b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
p

ggsave(plot = p, filename = "Fig2.png", height = 25, width = 20, unit = "cm")
