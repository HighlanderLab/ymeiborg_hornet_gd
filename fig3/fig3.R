##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig3")
source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig3a/Fig3a.Rdata")
load(file = "../fig3b/Fig3b.Rdata")
load(file = "../fig3c/Fig3c.Rdata")
load(file = "../fig3d/Fig3d.Rdata")
load(file = "../fig3e/Fig3e.Rdata")
load(file = "../fig3f/Fig3f.Rdata")

#########################################
########## Plot plots ###################
#########################################

fig3 <- (fig3a | fig3b) /
  (fig3c | fig3d) /
  (fig3e | fig3f) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14), 
        legend.position='bottom', 
        legend.box = "vertical")
fig3

ggsave(plot = fig3, filename = "Fig3.png", height = 25, width = 20, unit = "cm")
