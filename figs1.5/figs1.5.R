##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig1.5")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../figs1a/FigS1.5a.Rdata")
load(file = "../figs1b/FigS1.5b.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs1.5 <- figs1.5a / figs1.5b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14))
figs1.5

ggsave(plot = figs1.5, filename = "FigS1.5.png", height = 18, width = 20, unit = "cm")