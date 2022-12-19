##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs3b")

source("../model_function.R")

#########################################
########## Load data ####################
#########################################

filenames <- list.files(pattern="FigS3b_[0-9]*.Rdata", full.names=TRUE)
load(filenames[1])
allData <- as_tibble(modelOutput)

for (index in 2:length(filenames)){
  load(filenames[index])
  modelOutput <- as_tibble(modelOutput)
  allData <- rbind(allData, modelOutput)
}

modelOutput <- allData %>%
  arrange(gdSex, strategy) %>%
  mutate(repetitions = rep(1:100, each = 26, times = 8))

#########################################
########## Plot plots ###################
#########################################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5),
        legend.title=element_text(size=12),
        legend.position = "bottom",
        legend.justification = "center",
        axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5))

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
                        pHMort = factor(pHMort),
                        cutRate = factor(cutRate))

haploRep <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep$strategy = factor(haploRep$strategy, 
                           levels = c(1,2,3,4),
                           labels = c("Neutral","Male infertility","Female infertility","Both-sex infertility"))
haploRep$`Release` = factor(haploRep$`Release`, 
                           levels = c(1,2),
                           labels = c("Females","Males"))

filter(haploRep, 
       (generation == 25 & Release == "Males" & strategy == "Neutral") & 
         (Allele == "GD" & Frequency == 0)) |> nrow() -> maleFail

figs3b <- ggplot(data = haploRep) +
  facet_grid(
    `Release` ~ strategy,
    scales = "fixed",
    labeller = labeller(.cols = label_value, .rows = label_both)
  ) +
  geom_line(aes(
    x = generation,
    y = Frequency,
    group = interaction(Allele, repetitions),
    colour = Allele
  )) +
  scale_colour_manual(values = alpha(colour = met.brewer("Greek", 5)[-2], 
                                     alpha = 0.15)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlab("Generation") +
  ggtitle("European paper wasp") +
  PaperTheme
figs3b

suppRateData <- select(modelOutput, generation, repetitions, strategy, gdSex, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(strategy, gdSex) %>%
  summarise(suppressionRate = sum(suppressed)/100)

LastGenData <- select(modelOutput, generation, repetitions, strategy, gdSex, popSizeF) %>%
  filter(popSizeF > 0) %>%
  group_by(repetitions, strategy, gdSex) %>%
  filter(generation == max(generation)) %>%
  ungroup() %>%
  group_by(strategy, gdSex) %>%
  summarise(lastGenMean = mean(generation),
            lastGenSD = sd(generation))

#########################################
########## Save model ###################
#########################################

ggsave(plot = figs3b, filename = "FigS3b.png", height = 12, width = 20, unit = "cm")

save(modelOutput, figs3b, file = "FigS3b.Rdata")

