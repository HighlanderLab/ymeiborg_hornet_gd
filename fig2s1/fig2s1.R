library(dplyr)
library(patchwork)

setwd("/scratch/bell/ymeiborg/fig2s1/")

load(file = "../fig2as/Fig2aS1_plot.Rdata")
load(file = "../fig2bs1/Fig2bS1_plot.Rdata")

p <- fig2bs1 / fig2bs1 +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
p

ggsave(plot = p, filename = "Fig2s1.png", height = 25, width = 20, unit = "cm")
