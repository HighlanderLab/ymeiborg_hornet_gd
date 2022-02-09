##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/figs3")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- c(2,6,10,14,18) #average of female progeny per queen
input$meanMalProgeny <- c(2,6,10,14,18) #average of male progeny per queen
input$meanFemMatings <- c(0.1,0.15,0.2,0.25,0.3) #average number of times females mate
input$meanMalMatings <- c(0.1,0.15,0.2,0.25,0.3) #average number of times males mate
input$maxFemMatings <- c(1,2,3,4) #maximum number of times females mate
input$maxMalMatings <- c(1,2,3,4) #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- c(2,6,10,14,18) #growth rate of the population
input$N <- 1000 #size of start WT population
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$multiplex <- 1 #how many multiplexes in the gene drives, not used currently
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- c(0, 0.02) #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- c(1, 0.95) #propability CRISPR cuts the opposite DNA strand
input$pHMort <- c(0,0.1) #mortality of gene drive carriers.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter(((meanFemProgeny==10 & meanMalProgeny==10 & meanFemMatings==0.2 & meanMalMatings==0.2 & rmax==10 & maxFemMatings==2 & maxMalMatings==2) |
            (meanFemProgeny==10 & meanMalProgeny==10 & meanFemMatings==0.2 & meanMalMatings==0.2 & rmax==10 & maxFemMatings==2) |
            (meanFemProgeny==10 & meanMalProgeny==10 & meanFemMatings==0.2 & meanMalMatings==0.2 & rmax==10 & maxMalMatings==2) |
            (meanFemProgeny==10 & meanMalProgeny==10 & meanFemMatings==0.2 & meanMalMatings==0.2 & maxFemMatings==2 & maxMalMatings==2) |
            (meanFemProgeny==10 & meanMalProgeny==10 & meanFemMatings==0.2 & rmax==10 & maxFemMatings==2 & maxMalMatings==2) |
            (meanFemProgeny==10 & meanMalProgeny==10 & meanMalMatings==0.2 & rmax==10 & maxFemMatings==2 & maxMalMatings==2) |
            (meanFemProgeny==meanMalProgeny & meanFemMatings==0.2 & meanMalMatings==0.2 & rmax==10 & maxFemMatings==2 & maxMalMatings==2)) &
           ((pnhej == 0 & cutRate == 1 & pHMort == 0) |
              (pnhej == 0.02 & cutRate == 0.95 & pHMort == 0.1)))

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

save(modelOutput, file = "FigS3.Rdata")

##################
## ---- Plots ----
##################

modelOutput <- as_tibble(modelOutput)

modelOutput <- rowwise(modelOutput) %>% 
  select(-meanMalProgeny) %>%
  rename(meanProgeny = meanFemProgeny) %>% 
  mutate(variableRange = case_when(
    meanProgeny != 10 ~ "meanProgeny",
    meanFemMatings != 0.2 ~ "meanFemMatings",
    meanMalMatings != 0.2 ~ "meanMalMatings",
    maxFemMatings != 2 ~ "maxFemMatings",
    maxMalMatings != 2 ~ "maxMalMatings",
    rmax != 10 ~ "rmax",
    TRUE ~ "default"
  )) 

modelOutputOptimal <- filter(modelOutput, cutRate == 1)
modelOutputRealistic <- filter(modelOutput, cutRate == 0.95)

year25optimal <- filter(modelOutputOptimal, generation == 25)
year25realistic <- filter(modelOutputRealistic, generation == 25)


ntimes <- nrow(year25optimal)/max(inputs$repetitions)
year25optimal$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25optimal$paramSet <- factor(year25optimal$paramSet)
year25realistic$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25realistic$paramSet <- factor(year25realistic$paramSet)

################## Realistic GD scenario ##################

suppressed <- rowwise(year25realistic) %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressed)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), cond = "opt") %>%
  arrange(name)

suppressed$name <- factor(suppressed$name, 
                          levels = c('rmax','meanFemMatings','meanMalMatings',
                                     'maxFemMatings','maxMalMatings',
                                     'meanProgeny','default'), 
                          labels = c('Max growth rate','Mean female matings','Mean male matings',
                                     'Max female matings','Max male matings',
                                     'Mean progeny','Default'))
suppressed <- suppressed[order(suppressed$name),]

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        plot.title=element_text(size=14, hjust=0.5),
        legend.position="none",
        axis.title=element_text(size=12))

p1 <- ggplot(data = suppressed) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = suppressionRate, colour = name)) +
  geom_point(aes(x = value, y = suppressionRate, colour = name)) +
  scale_colour_manual(values=met.brewer("Greek", 6)) +
  ylim(c(0,1)) +
  ylab("Suppression rate") +
  xlab("Parameter value") + 
  PaperTheme +
  ggtitle("Realistic gene drive")
p1

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
                                        'meanProgeny','default'), 
                       labels = c('Max growth rate','Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny','Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny','default'), 
                                labels = c('Max growth rate','Mean female matings','Mean male matings',
                                           'Max female matings','Max male matings',
                                           'Mean progeny','Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStatsRealistic <- lastGenStatsRealistic[order(lastGenStatsRealistic$name),]

p2 <- ggplot(data = lastGenStatsRealistic) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = mLastGen, colour = name)) +
  geom_point(aes(x = value, y = mLastGen, colour = name)) +
  geom_errorbar(aes(x = value, y = mLastGen, ymin = lower, ymax = upper), 
                size = 0.1, width = 0.1) +
  ylim(c(0,25)) +
  ylab("Last viable generation") +
  xlab("Parameter value") +
  scale_colour_manual(values=met.brewer("Greek", 6)) +
  PaperTheme +
  ggtitle("Realistic gene drive")
p2

################## Optimal GD scenario ##################

suppressed <- rowwise(year25optimal) %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressed)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), cond = "opt") %>%
  arrange(name)

suppressed$name <- factor(suppressed$name, 
                          levels = c('rmax','meanFemMatings','meanMalMatings',
                                     'maxFemMatings','maxMalMatings',
                                     'meanProgeny','default'), 
                          labels = c('Max growth rate','Mean female matings','Mean male matings',
                                     'Max female matings','Max male matings',
                                     'Mean progeny','Default'))
suppressed <- suppressed[order(suppressed$name),]

p3 <- ggplot(data = suppressed) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = suppressionRate, colour = name)) +
  geom_point(aes(x = value, y = suppressionRate, colour = name)) +
  scale_colour_manual(values=met.brewer("Greek", 6)) +
  ylim(c(0,1)) +
  ylab("Suppression rate") +
  xlab("Parameter value") + 
  PaperTheme +
  ggtitle("Optimal gene drive")
p3

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
                                        'meanProgeny','default'), 
                       labels = c('Max growth rate','Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny','Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny','default'), 
                                labels = c('Max growth rate','Mean female matings','Mean male matings',
                                           'Max female matings','Max male matings',
                                           'Mean progeny','Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStatsOptimal <- lastGenStatsOptimal[order(lastGenStatsOptimal$name),]

p4 <- ggplot(data = lastGenStatsOptimal) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = mLastGen, colour = name)) +
  geom_point(aes(x = value, y = mLastGen, colour = name)) +
  geom_errorbar(aes(x = value, y = mLastGen, ymin = lower, ymax = upper), 
                size = 0.1, width = 0.1) +
  ylim(c(0,25)) +
  ylab("Last viable generation") +
  xlab("Parameter value") +
  scale_colour_manual(values=met.brewer("Greek", 6)) +
  PaperTheme +
  ggtitle("Optimal gene drive")
p4

p <- (p1 + p2) / (p3 + p4) + 
  plot_annotation(tag_levels = "A")
p

ggsave(plot = p, filename = "FigS3.png", height = 30, width = 25, unit = "cm")

