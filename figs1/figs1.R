##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig0")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig2a/FigS1a.Rdata")
load(file = "../fig2b/FigS1b.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs1 <- figs1a / figs1b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14))
figs1

ggsave(plot = figs1, filename = "FigS1.png", height = 18, width = 20, unit = "cm")
