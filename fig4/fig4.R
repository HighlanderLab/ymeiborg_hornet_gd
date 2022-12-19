

source("../model_function.R")

#############################
######## load data ##########
#############################

load("../fig4ab/Fig4ab.Rdata")
load("../fig4cd/Fig4cd.Rdata")

#############################
######## make plot ##########
#############################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        title=element_text(size=14, hjust=0.5), 
        plot.title=element_text(size=14, hjust=0.5),
        #legend.position="none",
        axis.line = element_line(),
        axis.title=element_text(size=11),
        axis.text = element_text(size = 6.5))

fig4a <- fig4a + PaperTheme
fig4b <- fig4b + PaperTheme
fig4c <- fig4c + PaperTheme
fig4d <- fig4d + PaperTheme

fig4 <- (fig4a | fig4b) /
  (fig4c | fig4d) +
  plot_layout(guides = 'collect') &
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 12), 
        legend.position='bottom', 
        legend.box = "vertical")


ggsave(plot = fig4, filename = "Fig4.png", height = 25, width = 20, unit = "cm")