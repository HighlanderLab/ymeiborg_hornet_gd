##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/fig3c")

source("../model_function.R")

#############################
######## load data ##########
#############################

################## Fig 3a ##################

load("../fig3a/Fig3a.Rdata")

modelOutput <- as_tibble(modelOutput)

modelOutput <- mutate(modelOutput,
                      strategy = factor(strategy),
                      N = factor(N),
                      nGD = factor(nGD),
                      k = factor(k),
                      gdSex = factor(gdSex),
                      rmax = factor(rmax),
                      generations = factor(generations),
                      repetitions = factor(repetitions),
                      hEffect = factor(hEffect),
                      cutRate = factor(cutRate))

heatMapData <- select(modelOutput, generation, repetitions, pnhej, cutRate, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(pnhej, cutRate) %>%
  summarise(suppressionRate = sum(suppressed)/10)

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5),
        legend.title=element_text(size=12),
        legend.position = "bottom",
        legend.justification = "center",
        axis.title=element_text(size=12))

p1 <- ggplot(data = heatMapData) +
  geom_raster(aes(x = pnhej, y = cutRate, fill = suppressionRate)) +
  scale_fill_viridis(name = "Supression Rate") +
  xlab("P(Nonhomologous endjoining)") +
  ylab("P(Cutting)") +
  PaperTheme
p1

################## Fig 6b ##################

filenames <- list.files(pattern="Fig3c_[0-9]*.Rdata", full.names=TRUE)
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

modelOutput <- mutate(modelOutput,
                        strategy = factor(strategy),
                        N = factor(N),
                        nGD = factor(nGD),
                        k = factor(k),
                        gdSex = factor(gdSex),
                        rmax = factor(rmax),
                        generations = factor(generations),
                        repetitions = factor(repetitions),
                        hEffect = factor(hEffect),
                        pHMort = factor(pHMort))

heatMapData <- select(modelOutput, generation, repetitions, pnhej, pHMort, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(pnhej, pHMort) %>%
  summarise(suppressionRate = sum(suppressed)/10)

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5),
        legend.title=element_text(size=12),
        legend.position = "bottom",
        legend.justification = "center",
        axis.title=element_text(size=12))

p2 <- ggplot(data = heatMapData) +
  geom_raster(aes(x = pnhej, y = pHMort, fill = suppressionRate)) +
  scale_fill_viridis(name = "Supression Rate", limits = c(0,1)) +
  xlab("P(Nonhomologous endjoining)") +
  ylab("P(GD heterozygote mortality)") +
  PaperTheme
p2

ggsave(plot = p2, filename = "fig3b.png", height = 11, width = 10, unit = "cm")


p <- p1 + p2 + 
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  plot_annotation(tag_levels = 'A')
p

ggsave(plot = p, filename = "fig3.png", height = 11, width = 20, unit = "cm")
