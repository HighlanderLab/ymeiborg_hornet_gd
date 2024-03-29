##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/ymeiborg_hornet_gd/fig3c")
source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1 # number of reps
input$generations <- 50 #runtime of the simulation
input$meanFemProgeny <- 300 #average of female progeny per queen
input$meanMalProgeny <- 300 #average of male progeny per queen
input$meanFemMatings <- 3.275 #average number of times females mate
input$meanMalMatings <- 0.9 #average number of times males mate
input$maxFemMatings <- 4 #maximum number of times females mate
input$maxMalMatings <- 2 #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- 10 #growth rate of the population
input$N <- 1000 #size of start WT population
input$winterMort <- 0.97
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0 #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- seq(0.91, 1, 0.01) #probability CRISPR cuts the opposite DNA strand
input$pHMort <- seq(0, 0.25, 0.025) #mortality of gene drive carriers.
input$pFunctionalRepair <- 0 #probability a resistance allele forms after non-homologous end-joining.
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

save(modelOutput, file = "Fig3c_1_3.Rdata")

