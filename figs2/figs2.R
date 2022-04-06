##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/figs2")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- c(10,155,300,445,590) #average of female progeny per queen
input$meanMalProgeny <- c(10,155,300,445,590) #average of male progeny per queen
input$meanFemMatings <- c(1.275,2.275,3.275,4.275,5.275) #average number of times females mate
input$meanMalMatings <- c(0.3,0.6,0.9,1.2,1.5) #average number of times males mate
input$maxFemMatings <- c(2,3,4,5,6) #maximum number of times females mate
input$maxMalMatings <- c(1,2,3,4) #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- c(2,6,10,14,18) #growth rate of the population
input$N <- 1000 #size of start WT population
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- c(0, 0.02) #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- c(1, 0.95, 0.97) #propability CRISPR cuts the opposite DNA strand
input$pHMort <- c(0, 0.1, 0.15) #mortality of gene drive carriers.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter(((meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==meanMalProgeny & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2)) &
            ((pnhej == 0 & cutRate == 1 & pHMort == 0) |
             (pnhej == 0.02 & cutRate == 0.95 & pHMort == 0.1) |
             (pnhej == 0 & cutRate == 0.97 & pHMort == 0.15)))

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

save(modelOutput, file = "FigS2.Rdata")

#################### Prepare data ##################

modelOutput <- as_tibble(modelOutput)

modelOutput <- rowwise(modelOutput) %>% 
  select(-meanMalProgeny) %>%
  rename(meanProgeny = meanFemProgeny) %>%
  mutate(variableRange = case_when(
    meanProgeny != 300 ~ "meanProgeny",
    meanFemMatings != 3.275 ~ "meanFemMatings",
    meanMalMatings != 0.9 ~ "meanMalMatings",
    maxFemMatings != 4 ~ "maxFemMatings",
    maxMalMatings != 2 ~ "maxMalMatings",
    rmax != 10 ~ "rmax",
    TRUE ~ "default"
  ))

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        plot.title=element_text(size=14, hjust=0.5),
        #legend.position="none",
        axis.title=element_text(size=12))

modelOutputOptimal <- filter(modelOutput, cutRate == 1)
modelOutputRealistic <- filter(modelOutput, cutRate == 0.95)
modelOutputMid <- filter(modelOutput, cutRate == 0.97)

year25optimal <- filter(modelOutputOptimal, generation == 25)
year25realistic <- filter(modelOutputRealistic, generation == 25)
year25mid<- filter(modelOutputMid, generation == 25)

ntimes <- nrow(year25optimal)/max(inputs$repetitions)
year25optimal$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25optimal$paramSet <- factor(year25optimal$paramSet)
year25realistic$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25realistic$paramSet <- factor(year25realistic$paramSet)
year25mid$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25mid$paramSet <- factor(year25mid$paramSet)

################## Realistic GD scenario ##################

suppressedRealistic <- rowwise(year25realistic) %>%
  mutate(suppressedRealistic = case_when(popSizeF == 0 ~ 1,
                                         popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressedRealistic)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), condition = "Realistic") %>%
  arrange(name)

suppressedRealistic$name <- factor(suppressedRealistic$name, 
                                   levels = c('rmax','meanFemMatings','meanMalMatings',
                                              'maxFemMatings','maxMalMatings',
                                              'meanProgeny', 'default'), 
                                   labels = c('Max growth rate','Mean female matings','Mean male matings',
                                              'Max female matings','Max male matings',
                                              'Mean progeny', 'Default'))
suppressedRealistic <- suppressedRealistic[order(suppressedRealistic$name),]

modelOutputRealistic$paramSet <- rep(1:ntimes, each = nrow(modelOutputRealistic)/ntimes)
modelOutputRealistic$paramSet <- factor(modelOutputRealistic$paramSet)

lastGenRealistic <- modelOutputRealistic %>%
  filter(popSizeF != 0) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == max(generation))

lastGenStatsRealistic <- group_by(lastGenRealistic, paramSet) %>%
  summarise(mLastGen = mean(generation),
            semLastGen = sd(generation)/sqrt(10)) 

lastGenStatsRealistic <- left_join(x = lastGenStatsRealistic, 
                                   y = lastGenRealistic, by = "paramSet") %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet:semLastGen, meanProgeny:rmax, -k, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name, levels = c('rmax','meanFemMatings','meanMalMatings',
                                        'maxFemMatings','maxMalMatings',
                                        'meanProgeny', 'default'), 
                       labels = c('Max growth rate','Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny', 'Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny', 'default'), 
                                labels = c('Max growth rate','Mean female matings','Mean male matings',
                                           'Max female matings','Max male matings',
                                           'Mean progeny', 'Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStatsRealistic <- lastGenStatsRealistic[order(lastGenStatsRealistic$name),] %>%
  mutate(condition = "Realistic") %>%
  group_by(name) %>%
  mutate(wd = (max(value) - min(value)) * 0.1)

################## Optimal GD scenario ##################

suppressedOptimal <- rowwise(year25optimal) %>%
  mutate(suppressedOptimal = case_when(popSizeF == 0 ~ 1,
                                       popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressedOptimal)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), condition = "Optimal") %>%
  arrange(name)

suppressedOptimal$name <- factor(suppressedOptimal$name, 
                                 levels = c('rmax','meanFemMatings','meanMalMatings',
                                            'maxFemMatings','maxMalMatings',
                                            'meanProgeny', 'default'), 
                                 labels = c('Max growth rate','Mean female matings','Mean male matings',
                                            'Max female matings','Max male matings',
                                            'Mean progeny', 'Default'))
suppressedOptimal <- suppressedOptimal[order(suppressedOptimal$name),]

modelOutputOptimal$paramSet <- rep(1:ntimes, each = nrow(modelOutputOptimal)/ntimes)
modelOutputOptimal$paramSet <- factor(modelOutputOptimal$paramSet)

lastGenOptimal <- modelOutputOptimal %>%
  filter(popSizeF != 0) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == max(generation))

lastGenStatsOptimal <- group_by(lastGenOptimal, paramSet) %>%
  summarise(mLastGen = mean(generation),
            semLastGen = sd(generation)/sqrt(10)) 

lastGenStatsOptimal <- left_join(x = lastGenStatsOptimal, 
                                 y = lastGenOptimal, by = "paramSet") %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet:semLastGen, meanProgeny:rmax, -k, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name, levels = c('rmax','meanFemMatings','meanMalMatings',
                                        'maxFemMatings','maxMalMatings',
                                        'meanProgeny', 'default'), 
                       labels = c('Max growth rate','Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny', 'Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny', 'default'), 
                                labels = c('Max growth rate','Mean female matings','Mean male matings',
                                           'Max female matings','Max male matings',
                                           'Mean progeny', 'Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStatsOptimal <- lastGenStatsOptimal[order(lastGenStatsOptimal$name),] %>%
  mutate(condition = "Optimal") %>%
  group_by(name) %>%
  mutate(wd = (max(value) - min(value)) * 0.1)

################## Intermediate GD scenario ##################

suppressedMid <- rowwise(year25mid) %>%
  mutate(suppressedMid = case_when(popSizeF == 0 ~ 1,
                                   popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressedMid)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), condition = "Intermediate") %>%
  arrange(name)

suppressedMid$name <- factor(suppressedMid$name, 
                             levels = c('rmax','meanFemMatings','meanMalMatings',
                                        'maxFemMatings','maxMalMatings',
                                        'meanProgeny', 'default'), 
                             labels = c('Max growth rate','Mean female matings','Mean male matings',
                                        'Max female matings','Max male matings',
                                        'Mean progeny', 'Default'))
suppressedMid <- suppressedMid[order(suppressedMid$name),]

modelOutputMid$paramSet <- rep(1:ntimes, each = nrow(modelOutputMid)/ntimes)
modelOutputMid$paramSet <- factor(modelOutputMid$paramSet)

lastGenMid <- modelOutputMid %>%
  filter(popSizeF != 0) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == max(generation))

lastGenStatsMid <- group_by(lastGenMid, paramSet) %>%
  summarise(mLastGen = mean(generation),
            semLastGen = sd(generation)/sqrt(10)) 

lastGenStatsMid <- left_join(x = lastGenStatsMid, 
                             y = lastGenMid, by = "paramSet") %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet:semLastGen, meanProgeny:rmax, -k, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name, levels = c('rmax','meanFemMatings','meanMalMatings',
                                        'maxFemMatings','maxMalMatings',
                                        'meanProgeny', 'default'), 
                       labels = c('Max growth rate','Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny', 'Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny', 'default'), 
                                labels = c('Max growth rate','Mean female matings','Mean male matings',
                                           'Max female matings','Max male matings',
                                           'Mean progeny', 'Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStatsMid <- lastGenStatsMid[order(lastGenStatsMid$name),] %>%
  mutate(condition = "Intermediate") %>%
  group_by(name) %>%
  mutate(wd = (max(value) - min(value)) * 0.1)

################## combine data ##################

suppressed <- bind_rows(suppressedRealistic, suppressedOptimal, suppressedMid) %>%
  mutate(condition = factor(condition, levels = c("Optimal", "Intermediate", "Realistic"),
                            labels = c("Optimal", "Intermediate", "Realistic")))

lastGenStats <- bind_rows(lastGenStatsRealistic, lastGenStatsOptimal, lastGenStatsMid) %>%
  mutate(condition = factor(condition, levels = c("Optimal", "Intermediate", "Realistic"),
                            labels = c("Optimal", "Intermediate", "Realistic")))

################## plots ##################

p1 <- ggplot(data = suppressed) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = suppressionRate, col= condition)) +
  geom_point(aes(x = value, y = suppressionRate, col = condition)) +
  scale_colour_manual(values=met.brewer("Greek", n = length(unique(suppressed$condition))),
                      name = "Gene Drive Conditions") +
  ylim(c(0,1)) +
  ylab("Suppression rate") +
  xlab("Parameter value") + 
  PaperTheme

p2 <- ggplot(data = lastGenStats) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = mLastGen, colour = condition)) +
  geom_point(aes(x = value, y = mLastGen, colour = condition)) +
  geom_errorbar(aes(x = value, y = mLastGen, ymin = lower, ymax = upper), 
                size = 0.1, width = lastGenStats$wd) +
  scale_colour_manual(values=met.brewer("Greek", n = length(unique(lastGenStats$condition))),
                      name = "Gene Drive Conditions") +
  ylim(c(0,25)) +
  ylab("Last viable generation") +
  xlab("Parameter value") +
  PaperTheme

p <- p1 + p2 + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = 'collect') & theme(legend.position='bottom')
p

ggsave(plot = p, filename = "FigS2.png", height = 20, width = 25, unit = "cm")

