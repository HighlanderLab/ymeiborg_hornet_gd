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
input$gdSex <- "F" #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$strategy <- 3 #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0 #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- 1 #propability CRISPR cuts the opposite DNA strand
input$pHMort <- 0 #mortality of gene drive carriers.
inputs <- expand.grid(input)

inputs <- inputs %>%
  filter((meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & meanMalMatings==0.9 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanFemMatings==3.275 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==300 & meanMalProgeny==300 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2) |
            (meanFemProgeny==meanMalProgeny & meanFemMatings==3.275 & meanMalMatings==0.9 & rmax==10 & maxFemMatings==4 & maxMalMatings==2)) 

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

save(modelOutput, file = "FigS2_1.Rdata")

