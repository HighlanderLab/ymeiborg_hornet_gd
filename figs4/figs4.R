##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/figs4")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- c(10, 20, 40, 80, 160, 320, 640) #average of female progeny per queen
input$meanMalProgeny <- c(10, 20, 40, 80, 160, 320, 640) #average of male progeny per queen
input$meanFemMatings <- c(3.275, 0.2) #average number of times females mate
input$meanMalMatings <- 0.9 #average number of times males mate
input$maxFemMatings <- c(4, 2) #maximum number of times females mate
input$maxMalMatings <- 3 #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- 10 #growth rate of the population
input$N <- 1000 #size of start WT population
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$multiplex <- 1 #how many multiplexes in the gene drives, not used currently
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- c(0, 0.02) #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- c(1, 0.95, 0.97) #propability CRISPR cuts the opposite DNA strand
input$pHMort <- c(0, 0.1, 0.15) #mortality of gene drive carriers.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter((meanFemProgeny==meanMalProgeny & meanFemMatings==3.275 & maxFemMatings==4 | 
           meanFemProgeny==meanMalProgeny & meanFemMatings==0.2 & maxFemMatings==2) &
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

save(modelOutput, file = "FigS4.Rdata")

#################### Prepare data ##################

# modelOutput <- as_tibble(modelOutput)
# 
# modelOutput <- rowwise(modelOutput) %>% 
#   select(-meanMalProgeny) %>%
#   rename(meanProgeny = meanFemProgeny) %>%
#   mutate(condition = case_when(
#     maxFemMatings == 2 ~ "European Paper Wasp",
#     maxFemMatings == 4 ~ "Asian Hornet",
#   ))
# 
# PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
#   theme(strip.background = element_blank(),
#         panel.grid = element_blank(),
#         title=element_text(size=14, hjust=0.5), 
#         plot.title=element_text(size=14, hjust=0.5),
#         axis.title=element_text(size=12))
# 
# modelOutputOptimal <- filter(modelOutput, cutRate == 1)
# modelOutputRealistic <- filter(modelOutput, cutRate == 0.95)
# modelOutputMid <- filter(modelOutput, cutRate == 0.97)
# 
# year25optimal <- filter(modelOutputOptimal, generation == 25)
# year25realistic <- filter(modelOutputRealistic, generation == 25)
# year25mid<- filter(modelOutputMid, generation == 25)
# 
# ntimes <- nrow(year25optimal)/max(inputs$repetitions)
# year25optimal$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
# year25optimal$paramSet <- factor(year25optimal$paramSet)
# year25realistic$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
# year25realistic$paramSet <- factor(year25realistic$paramSet)
# year25mid$paramSet <- rep(1:ntimes, each = max(inputs$repetitions))
# year25mid$paramSet <- factor(year25mid$paramSet)

################## Realistic GD scenario ##################

