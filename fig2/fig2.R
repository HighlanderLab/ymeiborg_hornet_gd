##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig2")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig2a/Fig2a.Rdata")
load(file = "../fig2b/Fig2b.Rdata")

#########################################
########## Plot plots ###################
#########################################

fig2 <- fig2a / fig2b +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 17))
fig2

ggsave(plot = fig2, filename = "Fig2.png", height = 25, width = 20, unit = "cm")
