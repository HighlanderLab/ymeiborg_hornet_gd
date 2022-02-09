##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/fig3d")
source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 50 #runtime of the simulation
input$meanFemProgeny <- 20 #average of female progeny per queen
input$meanMalProgeny <- 20 #average of male progeny per queen
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
input$pnhej <- seq(0, 0.02, 0.001) #probability of non-homologous end joining, determines the resistance alleles (0.02 in mosquitos)
input$cutRate <- 1 #probability CRISPR cuts the opposite DNA strand
input$pHMort <- seq(0, 0.5, 0.025) #mortality of gene drive carriers.
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

save(modelOutput, file = "Fig3d.Rdata")

#########################################
########## Plot plots ###################
#########################################

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12))

modelOutput <- as_tibble(modelOutput)

modelOutput <- mutate(modelOutput,
                      strategy = factor(strategy),
                      N = factor(N),
                      nGD = factor(nGD),
                      k = factor(k),
                      gdSex = factor(gdSex),
                      rmax = factor(rmax),
                      generations = factor(generations),
                      repetitions = factor(repetitions),
                      pHMort = factor(pHMort),
                      cutRate = factor(cutRate))

heatMapData <- select(modelOutput, generation, repetitions, pnhej, pHMort, popSizeF) %>%
  filter(generation == max(generation)) %>%
  rowwise() %>%
  mutate(suppressed = case_when(popSizeF == 0 ~ 1,
                                popSizeF > 0 ~ 0)) %>%
  group_by(pnhej, pHMort) %>%
  summarise(suppressionRate = sum(suppressed)/10)

fig3d <- ggplot(data = heatMapData) +
  geom_raster(aes(x = pnhej, y = pHMort, fill = suppressionRate)) +
  scale_fill_gradientn(colors=met.brewer("Greek"), limits = c(0,1), name = "Suppression rate") +
  xlab("P(Nonhomologous endjoining)") +
  ylab("P(GD heterozygote mortality)") +
  ggtitle("Asian hornet") +
  PaperTheme
fig3d

ggsave(plot = fig3d, filename = "Fig3d.png", height = 11, width = 10, unit = "cm")

#########################################
########## Save model ###################
#########################################

save(modelOutput, fig3d, file = "Fig3d.Rdata")
