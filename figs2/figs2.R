##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/figs2")

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
input$winterMort <- c(0.95, 0.96, 0.97, 0.98, 0.99) #winter mortalty
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0.02 #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- 0.95 #propability CRISPR cuts the opposite DNA strand
input$pHMort <-  0 #mortality of gene drive carriers.
input$pFunctionalRepair <- 0.01 #probability a resistance allele forms after non-homologous end-joining.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter((meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & winterMort == 0.97) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxMalMatings==2 & winterMort == 0.97) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & rmax==10 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97) |
           (meanFemProgeny==300 & meanMalProgeny==300 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97)
           (meanFemProgeny==meanMalProgeny & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97))
 
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
  select(paramSet, meanProgeny:maxMalMatings, rmax, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:rmax) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name)) %>%
  arrange(name)

suppressed$name <- factor(suppressed$name, 
                          levels = c('rmax','meanFemMatings','meanMalMatings',
                                     'maxFemMatings','maxMalMatings',
                                     'meanProgeny', 'default'), 
                          labels = c('Max growth rate','Mean female matings','Mean male matings',
                                     'Max female matings','Max male matings',
                                     'Mean progeny', 'Default'))
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
lastGenStats <- lastGenStats[order(lastGenStats$name),] %>%
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