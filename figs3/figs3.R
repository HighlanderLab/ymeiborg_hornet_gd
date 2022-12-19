##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs3")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../figs3a/FigS3a.Rdata")
load(file = "../figs3b/FigS3b.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs1 <- figs3a / figs3b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14))
figs1

ggsave(plot = figs1, filename = "FigS3.png", height = 24, width = 20, unit = "cm")
