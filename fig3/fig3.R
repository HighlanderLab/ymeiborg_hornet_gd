##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig3")
source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig3g/Fig3g.Rdata")
load(file = "../fig3h/fig3h.Rdata")
load(file = "../fig3e/Fig3e.Rdata")
load(file = "../fig3f/Fig3f.Rdata")

#########################################
########## Plot plots ###################
#########################################

fig3 <- (fig3g | fig3h) /
  (fig3e | fig3f) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14), 
        legend.position='bottom', 
        legend.box = "vertical")
fig3

ggsave(plot = fig3, filename = "Fig3.png", height = 20, width = 20, unit = "cm")
