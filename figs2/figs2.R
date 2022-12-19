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
input$nGD <- c(0, 100) #number of gene drive carrying animals to introduce
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
           (meanFemProgeny==300 & meanMalProgeny==300 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2 & winterMort == 0.97) |
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
    winterMort != 0.97 ~ "winterMort",
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
        panel.border = element_blank(), 
        title=element_text(size=14, hjust=0.5), 
        plot.title=element_text(size=14, hjust=0.5),
        axis.line = element_line(),
        axis.title=element_text(size=12),
        axis.text = element_text(size = 6.5))

modelOutputGD <- filter(modelOutput, nGD == 100)
modelOutputNoGD <- filter(modelOutput, nGD == 0)

year25GD <- filter(modelOutputGD, generation == 25)
year25NoGD <- filter(modelOutputNoGD, generation == 25)

ntimes <- nrow(year25GD)/max(inputs$repetitions)
year25GD$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25GD$paramSet <- factor(year25GD$paramSet)
year25NoGD$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
year25NoGD$paramSet <- factor(year25NoGD$paramSet)

################# Gene Drive data ###################

suppressedGD <- rowwise(year25GD) %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressed)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, winterMort, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:winterMort) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), condition = "Gene Drive") %>%
  arrange(name)

suppressedGD$name <- factor(suppressedGD$name, 
                            levels = c('rmax','meanFemMatings','meanMalMatings',
                                       'maxFemMatings','maxMalMatings',
                                       'meanProgeny', 'winterMort', 'default'), 
                            labels = c('Max growth rate','Mean female matings','Mean male matings',
                                       'Max female matings','Max male matings',
                                       'Mean progeny', 'Winter Mortality', 'Default'))
suppressedGD <- suppressedGD[order(suppressedGD$name),]

modelOutputGD$paramSet <- rep(1:ntimes, each = nrow(modelOutputGD)/ntimes)
modelOutputGD$paramSet <- factor(modelOutputGD$paramSet)

lastGenStatsGD <- modelOutputGD %>%
  filter((popSizeF == 0 & generation != 25) | (popSizeF != 0 & generation == 25)) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == min(generation)) %>%
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
  mutate(wd = (max(value) - min(value)) * 0.075)  %>%
  group_by(paramSet) %>%
  mutate(mLastGen = mean(generation),
         semLastGen = sd(generation)/sqrt(10),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen,
         condition = "Gene Drive")

################# No Gene Drive data ###################

suppressedNoGD <- rowwise(year25NoGD) %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(paramSet) %>%
  mutate(suppressionRate = sum(suppressed)/10) %>%
  distinct(paramSet, .keep_all = TRUE) %>%
  select(paramSet, meanProgeny:maxMalMatings, rmax, winterMort, suppressionRate, variableRange) %>%
  pivot_longer(cols = meanProgeny:winterMort) %>%
  filter(variableRange == name | variableRange == "default") %>%
  mutate(name = factor(name), condition = "No Gene Drive") %>%
  arrange(name)

suppressedNoGD$name <- factor(suppressedNoGD$name, 
                              levels = c('rmax','meanFemMatings','meanMalMatings',
                                         'maxFemMatings','maxMalMatings',
                                         'meanProgeny', 'winterMort', 'default'), 
                              labels = c('Max growth rate','Mean female matings','Mean male matings',
                                         'Max female matings','Max male matings',
                                         'Mean progeny', 'Winter Mortality', 'Default'))
suppressedNoGD <- suppressedNoGD[order(suppressedNoGD$name),]

modelOutputNoGD$paramSet <- rep(1:ntimes, each = nrow(modelOutputNoGD)/ntimes)
modelOutputNoGD$paramSet <- factor(modelOutputNoGD$paramSet)

lastGenStatsNoGD <- modelOutputNoGD %>%
  filter((popSizeF == 0 & generation != 25) | (popSizeF != 0 & generation == 25)) %>%
  group_by(paramSet, repetitions) %>%
  filter(generation == min(generation)) %>%
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
  mutate(wd = (max(value) - min(value)) * 0.075)  %>%
  group_by(paramSet) %>%
  mutate(mLastGen = mean(generation),
         semLastGen = sd(generation)/sqrt(10),
         lower = mLastGen - semLastGen,
         upper = mLastGen + semLastGen,
         condition = 'No Gene Drive')

############## combine data ################

suppressed <- bind_rows(suppressedGD, suppressedNoGD) %>%
  mutate(condition = factor(condition, levels = c("No Gene Drive", "Gene Drive"))) %>%
  arrange(condition)

lastGenStats <- bind_rows(lastGenStatsGD, lastGenStatsNoGD) %>%
  mutate(condition = factor(condition, levels = c("No Gene Drive", "Gene Drive"))) %>%
  arrange(condition)

################## plots ##################
cols <- c(met.brewer("Greek", n = 2)[2], met.brewer("Greek", n = 2)[1])

fig4a <- ggplot(data = suppressed) +
  facet_wrap(. ~ name, scales = "free_x",
             labeller = labeller(name = label_wrap_gen(width = 16))) +
  geom_line(aes(x = value, y = suppressionRate, col = condition)) +
  geom_point(aes(x = value, y = suppressionRate, col = condition)) +
  scale_colour_manual(values=cols,
                      name = "Gene drive conditions") +
  ylim(c(0,1)) +
  ylab("Suppression rate") +
  xlab("Parameter value") + 
  ggtitle("Asian hornet") +
  guides(colour = guide_legend(reverse = TRUE)) +
  PaperTheme
fig4a

fig4b <- ggplot(lastGenStats, aes(x = value, y = generation, col = condition)) +
  facet_wrap(. ~ name, scales = "free_x",
             labeller = labeller(name = label_wrap_gen(width = 16))) +
  stat_summary(fun = mean, geom = "line",lwd = 0.5, aes(col = condition)) +
  stat_summary(fun = mean, geom = "point", size = 1.5, aes(col = condition)) +
  geom_errorbar(aes(x = value, y = mLastGen, ymin = lower, ymax = upper), 
                size = 0.1, width = lastGenStats$wd, colour = "black") +
  scale_colour_manual(values = cols,
                      name = "Gene drive conditions",
                      breaks = c("No Gene Drive", "Gene Drive")) +
  scale_y_continuous(breaks = seq(0, 25, 5), 
                     labels = c(seq(0, 20, 5), "No supp.")) +
  ylab("Time to suppression (generation)") +
  xlab("Parameter value") +
  guides(colour = "none") +
  ggtitle("Asian hornet") +
  geom_point(shape = 1, alpha = 0.8) +
  PaperTheme
fig4b

p <- fig4a + fig4b + 
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = 'collect') & theme(legend.position='bottom') &
  theme(plot.tag = element_text(size = 14))
p


ggsave(plot = p, filename = "FigS2.png", height = 15, width = 25, unit = "cm")

save(modelOutput, fig4a, fig4b, file = "Fig4ab.Rdata")

