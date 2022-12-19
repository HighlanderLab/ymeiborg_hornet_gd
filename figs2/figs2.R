##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs4")

source("../model_function.R")

#############################
######## load data ##########
#############################

load(file = "../figs2a/Figs2a.Rdata")
load(file = "../figs2b/Figs2b.Rdata")
load(file = "../figs2c/Figs2c.Rdata")
load(file = "../figs2d/Figs2d.Rdata")
load(file = "../figs2e/Figs2e.Rdata")
load(file = "../figs2f/Figs2f.Rdata")
load(file = "../figs2g/Figs2g.Rdata")
load(file = "../figs2h/Figs2h.Rdata")

#########################################
########## Plot plots ###################
#########################################

figs2 <- (figs2a | figs2b) /
  (figs2c | figs2d) /
  (figs2e | figs2f) /
  (figs2g | fig3s2h) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 14), 
        legend.position='bottom', 
        legend.box = "vertical")
figs2

ggsave(plot = figs2, filename = "FigS2.png", height = 25, width = 20, unit = "cm")
