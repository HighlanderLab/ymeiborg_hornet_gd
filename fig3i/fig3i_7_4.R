##########################
########## Setup #########
##########################

setwd("/scratch/ymeiborg/ymeiborg_hornet_gd/fig3i")
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
input$pnhej <- 0.02 #probability of non-homologous end joining, determines the resistance alleles (0.2 in mosquitos)
input$cutRate <- seq(0.8, 0.9, 0.01) #probability CRISPR cuts the opposite DNA strand
input$pHMort <- 0 #mortality of gene drive carriers.
input$pFunctionalRepair <- lseq(0.000001, 0.1, 21)[12:21] #probability a resistance allele forms after non-homologous end-joining.
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

save(modelOutput, file = "Fig3i_7_4.Rdata")
