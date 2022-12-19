##########################
########## Setup #########
##########################

setwd("/scratch/ymeiborg/ymeiborg_hornet_gd/figs2e")
source("../model_function.R")

#############################
######## load data ##########
#############################

filenames <- list.files(pattern="FigS2e_[0-9]*_[0-9].Rdata", full.names=TRUE)
load(filenames[1])
allData <- as_tibble(modelOutput)

for (index in 2:length(filenames)){
  load(filenames[index])
  modelOutput <- as_tibble(modelOutput)
  modelOutput$repetitions <- index
  allData <- rbind(allData, modelOutput)
}

modelOutput <- allData

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

modelOutput <- mutate(modelOutput,
                        strategy = factor(strategy),
                        N = factor(N),
                        nGD = factor(nGD),
                        k = factor(k),
                        gdSex = factor(gdSex),
                        rmax = factor(rmax),
                        generations = factor(generations),
                        repetitions = factor(repetitions))

heatMapData <- select(modelOutput, generation, repetitions, cutRate, pFunctionalRepair, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(pFunctionalRepair, cutRate) %>%
  summarise(suppressionRate = sum(suppressed)/10)

figs2e <- ggplot(data = heatMapData) +
  geom_raster(aes(x = pFunctionalRepair, y = cutRate, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression rate") +
  scale_x_continuous(trans='log10') +
  xlab("P(Functional repair)") +
  ylab("P(Cutting)") +
  ggtitle("Asian hornet") +
  PaperTheme
figs2e

save(figs2e, file = "FigS2e.Rdata")
ggsave(plot = figs2e, filename = "FigS2e.png", height = 11, width = 10, unit = "cm")

