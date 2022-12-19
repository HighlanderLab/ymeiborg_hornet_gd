##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs1a")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:100 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- 300 #average of femlae progeny per queen
input$meanMalProgeny <- 300 #average of male progeny per queen
input$meanFemMatings <- 3.275 #average number of times females mate
input$meanMalMatings <- 0.9 #average number of times males mate
input$maxFemMatings <- 1 #maximum number of times females mate
input$maxMalMatings <- 3 #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- 10 #growth rate of the population
input$N <- 1000 #size of start WT population
input$winterMort <- 0.97 #winter mortality
input$gdSex <- c("F", "M") #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$multiplex <- 1 #how many multiplexes in the gene drives, not used currently
input$strategy <- c(1,2,3,4) #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0.02 #probability of non-homologous end joining, determines the resistance alleles (0.2 in mosquitos)
input$cutRate <- 0.95 #propability CRISPR cuts the opposite DNA strand
input$pHMort <- 0 #mortality of gene drive carriers.
input$pFunctionalRepair <- 0.01 #probability a resistance allele forms after non-homologous end-joining.
inputs <- expand.grid(input)

#########################################
########## Run model ####################
#########################################

# Track population information
generation <- 0:input$generations
var <- c('popSizeF', 'WT', 'GD', 'NF', 'RE', 'homoWT', 'homoGD', 'heteroGD')
modelOutput <- array(NA,
                     dim=c(length(generation),(length(var)+ncol(inputs)+1), nrow(inputs)),
                     dimnames=list(generation=generation, var=c("generation",colnames(inputs),var) , run=1:nrow(inputs)))

for (row in 1:nrow(inputs)){
  results <- modelHornets.comp(inputs[row,])
  modelOutput[,,row] <- unlist(results)
}

modelOutput <- apply(modelOutput, 2, c)

#########################################
########## Save model ###################
#########################################

save(modelOutput, file = "FigS1a.Rdata")

#########################################
########## Plot plots ###################
#########################################

filenames <- list.files(pattern="FigS1a_[0-9]*.Rdata", full.names=TRUE)
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
                        pHMort = factor(pHMort))

haploRep <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep$strategy = factor(haploRep$strategy, 
                           levels = c(1,2,3,4),
                           labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep$`Release` = factor(haploRep$`Release`, 
                           levels = c(1,2),
                           labels = c("Females","Males"))

filter(haploRep, 
       (generation == 25 & Release == "Males" & strategy == "Neutral") & 
         (Allele == "GD" & Frequency == 0)) |> nrow() -> maleFail

figs1a <- ggplot(data = haploRep) +
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
  ggtitle("Asian hornet") +
  PaperTheme
figs1a

popRep <- as_tibble(modelOutput) %>% 
  select(popSizeF, strategy, gdSex, generation, repetitions) %>%
  rename(Release = gdSex) %>%
  mutate(Release = factor(Release, 
                          levels = c(1,2), 
                          labels = c("Females","Males")),
         strategy= factor(strategy,
                          levels = c(1,2,3,4),
                          labels = c("Neutral","Male infertility",
                                     "Female infertility", 
                                     "Both-sex infertility")))

figs1.5a <- ggplot(data = popRep) +
  facet_grid(
    Release ~ strategy,
    scales = "fixed",
    labeller = labeller(.cols = label_value, .rows = label_both)
  ) +
  geom_line(aes(x = generation,
                y = popSizeF,
                group = repetitions),
            alpha = 0.1) +
  ylim(0, 1200) +
  xlab("Generation") +
  ylab("Female population size") +
  ggtitle("Asian hornet") +
  PaperTheme
figs1.5a

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

ggsave(plot = figs1a, filename = "FigS1a.png", height = 12, width = 20, unit = "cm")

save(modelOutput, figs1a, file = "FigS1a.Rdata")

ggsave(plot = figs1.5a, filename = "FigS1.5a.png", height = 12, width = 20, unit = "cm")

save(figs1.5a, file = "FigS1.5a.Rdata")

