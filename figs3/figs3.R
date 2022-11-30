##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs3")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- c(2,11,20,29,38) #average of female progeny per queen
input$meanMalProgeny <- c(2,11,20,29,38) #average of male progeny per queen
input$meanFemMatings <- c(0.1,0.15,0.2,0.25,0.3) #average number of times females mate
input$meanMalMatings <- c(0.3,0.6,0.9,1.2,1.5) #average number of times males mate
input$maxFemMatings <- c(1,2,3,4) #maximum number of times females mate
input$maxMalMatings <- c(1,2,3,4) #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- c(2,6,10,14,18) #growth rate of the population
input$N <- 1000 #size of start WT population
input$winterMort <- c(0.3, 0.4, 0.5, 0.6, 0.7) #winter mortality
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- c(0, 100) #number of gene drive carrying animals to introduce
input$multiplex <- 1 #how many multiplexes in the gene drives, not used currently
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0.02 #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- 0.95 #propability CRISPR cuts the opposite DNA strand
input$pHMort <- 0 #mortality of gene drive carriers.
input$pFunctionalRepair <- 0.01 #probability a resistance allele forms after non-homologous end-joining.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter((meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==2 & maxMalMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==2 & maxMalMatings==2) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & meanMalMatings==0.9 & rmax==10 & maxMalMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & meanMalMatings==0.9 & maxFemMatings==2 & maxMalMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanFemMatings==0.2 & rmax==10 & maxFemMatings==2 & maxMalMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==20 & meanMalProgeny==20 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==2 & maxMalMatings==2 & winterMort == 0.5) |
            (meanFemProgeny==meanMalProgeny & meanFemMatings==0.2 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==2 & maxMalMatings==2 & winterMort == 0.5))

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

#################### Prepare data ##################

modelOutput <- as_tibble(modelOutput)

modelOutput <- rowwise(modelOutput) %>% 
  select(-meanMalProgeny) %>%
  rename(meanProgeny = meanFemProgeny) %>%
  mutate(variableRange = case_when(
    winterMort != 0.5 ~ "winterMort",
    meanProgeny != 20 ~ "meanProgeny",
    meanFemMatings != 0.2 ~ "meanFemMatings",
    meanMalMatings != 0.9 ~ "meanMalMatings",
    maxFemMatings != 2 ~ "maxFemMatings",
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

year25 <- filter(modelOutput, generation == 25)

ntimes <- nrow(year25)/max(inputs$repetitions)
year25$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25$paramSet <- factor(year25$paramSet)

suppressed <- rowwise(year25) %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                         popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressed)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, winterMort, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:winterMort) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name)) %>%
  arrange(name)

suppressed$name <- factor(suppressed$name, 
                                   levels = c('rmax','meanFemMatings','meanMalMatings',
                                              'maxFemMatings','maxMalMatings',
                                              'meanProgeny', 'winterMort', 'default'), 
                                   labels = c('Max growth rate','Mean female matings','Mean male matings',
                                              'Max female matings','Max male matings',
                                              'Mean progeny', 'Winter Mortality', 'Default'))
suppressed <- suppressed[order(suppressed$name),]

modelOutput$paramSet <- rep(1:ntimes, each = nrow(modelOutput)/ntimes)
modelOutput$paramSet <- factor(modelOutput$paramSet)

lastGen <- modelOutput %>%
  filter(popSizeF != 0) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == max(generation))

lastGenStats <- group_by(lastGen, paramSet) %>%
  summarise(mLastGen = mean(generation),
            semLastGen = sd(generation)/sqrt(10)) 

lastGenStats <- left_join(x = lastGenStats, 
                                   y = lastGen, by = "paramSet") %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet:semLastGen, meanProgeny:rmax, winterMort, -k, variableRange) %>%
  pivot_longer(cols = meanProgeny:winterMort) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name, levels = c('rmax','meanFemMatings',
                                        'meanMalMatings', 'maxFemMatings',
                                        'maxMalMatings', 'meanProgeny', 'winterMort', 'default'), 
                       labels = c('Max growth rate',
                                  'Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny', 'Winter Mortality', 'Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny', 'winterMort', 'default'), 
                                labels = c('Max growth rate', 'Mean female matings',
                                           'Mean male matings', 'Max female matings',
                                           'Max male matings', 'Mean progeny', 
                                           'Winter Mortality', 'Default')),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen)
lastGenStats <- lastGenStats[order(lastGenStats$name),] %>%
  group_by(name) %>%
  mutate(wd = (max(value) - min(value)) * 0.1)

test <- lastGen %>%
  select(generation, repetitions, paramSet, meanProgeny:rmax, winterMort, -k, variableRange) %>%
  pivot_longer(cols = meanProgeny:winterMort) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name, levels = c('rmax','meanFemMatings',
                                        'meanMalMatings', 'maxFemMatings',
                                        'maxMalMatings', 'meanProgeny', 'winterMort', 'default'), 
                       labels = c('Max growth rate',
                                  'Mean female matings','Mean male matings',
                                  'Max female matings','Max male matings',
                                  'Mean progeny', 'Winter Mortality', 'Default')), 
         variableRange = factor(variableRange, levels = c('rmax','meanFemMatings','meanMalMatings',
                                                          'maxFemMatings','maxMalMatings',
                                                          'meanProgeny', 'winterMort', 'default'), 
                                labels = c('Max growth rate', 'Mean female matings',
                                           'Mean male matings', 'Max female matings',
                                           'Max male matings', 'Mean progeny', 
                                           'Winter Mortality', 'Default'))) %>%
  group_by(name) %>%
  mutate(wd = (max(value) - min(value)) * 0.1)

################## plots ##################

p1 <- ggplot(data = suppressed) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = suppressionRate)) +
  geom_point(aes(x = value, y = suppressionRate)) +
  ylim(c(0,1)) +
  ylab("Suppression rate") +
  xlab("Parameter value") + 
  PaperTheme

p2 <- ggplot(data = lastGenStats) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_line(aes(x = value, y = mLastGen)) +
  geom_point(aes(x = value, y = mLastGen)) +
  geom_errorbar(aes(x = value, y = mLastGen, ymin = lower, ymax = upper), 
                size = 0.1, width = lastGenStats$wd) +
  ylim(c(0,25)) +
  ylab("Last viable generation") +
  xlab("Parameter value") +
  PaperTheme

p3 <- ggplot(test, aes(y = generation, x = value, group = as.factor(paramSet))) +
  facet_wrap(. ~ name, scales = "free_x") +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE) +
  geom_jitter(width = test$wd) +
  PaperTheme
p3

p <- p1 + p2 + 
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  theme(plot.tag = element_text(size = 14))
p

ggsave(plot = p, filename = "FigS3.png", height = 20, width = 25, unit = "cm")

