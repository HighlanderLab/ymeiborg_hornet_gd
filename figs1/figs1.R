##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs1")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../figs1a/FigS1a.Rdata")
load(file = "../figs1b/FigS1b.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs1 <- figs1a / figs1b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
figs1

ggsave(plot = figs1, filename = "FigS1.png", height = 25, width = 20, unit = "cm")
