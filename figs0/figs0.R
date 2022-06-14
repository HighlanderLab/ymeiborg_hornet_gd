##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig0")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig2a/FigS0a.Rdata")
load(file = "../fig2b/FigS0b.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs0 <- figs0a / figs0b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
figs0

ggsave(plot = figs0, filename = "FigS0.png", height = 25, width = 20, unit = "cm")