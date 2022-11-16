##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig3l")
source("../model_function.R")

#############################
######## load data ##########
#############################

filenames <- list.files(pattern="fig3l_[0-9]*_[0-9].Rdata", full.names=TRUE)
load(filenames[1])
allData <- as_tibble(modelOutput)

for (index in 2:length(filenames)){
  load(filenames[index])
  modelOutput <- as_tibble(modelOutput)
  modelOutput$repetitions <- index
  allData <- rbind(allData, modelOutput)
}

modelOutput <- allData

load(file = "../dataGD.Rdata")
#########################################
########## Plot plots ###################
#########################################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        legend.box = "vertical",
        axis.title=element_text(size=12))


modelOutput <- as_tibble(modelOutput)

modelOutput <- mutate(modelOutput,
                        strategy = factor(strategy),
                        N = factor(N),
                        nGD = factor(nGD),
                        k = factor(k),
                        gdSex = factor(gdSex),
                        rmax = factor(rmax),
                        generations = factor(generations),
                        repetitions = factor(repetitions))

heatMapData <- select(modelOutput, generation, repetitions, pnhej, cutRate, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(pnhej, cutRate) %>%
  summarise(suppressionRate = sum(suppressed)/10)

fig3l <- ggplot(data = heatMapData) +
  geom_raster(aes(x = pnhej, y = cutRate, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression rate") +
  geom_point(data = dataGD, aes(x = Pnhej, y = pFunctionalRepair, shape = `Gene drive condition`), fill = "white") +
  scale_shape_manual(values = 21:23) +
  xlab("P(Non-homologous end-joining)") +
  ylab("P(Functional repair)") +
  ggtitle("European paper wasp") +
  PaperTheme
fig3l

ggsave(plot = fig3l, filename = "Fig3l.png", height = 11, width = 10, unit = "cm")

#########################################
########## Save model ###################
#########################################

save(modelOutput, fig3l, file = "Fig3l.Rdata")
