##########################
###########setup #########
##########################

setwd("/scratch/bell/ymeiborg/fig3f")
source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1 # number of reps
input$generations <- 50 #runtime of the simulation
input$meanFemProgeny <- 10 #average of female progeny per queen
input$meanMalProgeny <- 10 #average of male progeny per queen
input$meanFemMatings <- 0.2 #average number of times females mate
input$meanMalMatings <- 0.2 #average number of times males mate
input$maxFemMatings <- 2 #maximum number of times females mate
input$maxMalMatings <- 2 #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- 10 #growth rate of the population
input$N <- 1000 #size of start WT population
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0 #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- seq(0.8, 1, 0.01) #propability CRISPR cuts the opposite DNA strand
input$hEffect <- TRUE #logical, determines whether there is a fitness cost associated with the gene drive
input$pHMort <- seq(0, 0.5, 0.025) #only if hEffect == TRUE, mortality of gene drive carriers.
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

save(modelOutput, file = "Fig3f_1.Rdata")
