##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs4")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../fig3a/Fig3a.Rdata")
load(file = "../fig3b/Fig3b.Rdata")
load(file = "../fig3c/Fig3c.Rdata")
load(file = "../fig3d/Fig3d.Rdata")
load(file = "../fig3i/Fig3i.Rdata")
load(file = "../fig3j/Fig3j.Rdata")
load(file = "../fig3k/Fig3k.Rdata")
load(file = "../fig3l/Fig3l.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs4 <- (fig3a | fig3b) /
  (fig3c | fig3d) /
  (fig3i | fig3j) /
  (fig3k | fig3l) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14), 
        legend.position='bottom', 
        legend.box = "vertical")
figs4

ggsave(plot = figs4, filename = "FigS4.png", height = 25, width = 20, unit = "cm")
