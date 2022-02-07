##########################
########## Setup #########
##########################

setwd("/scratch/bell/ymeiborg/fig2as1")

source("../model_function.R")

#########################################
########## Set input parameters #########
#########################################

input <- list()

input$repetitions <- 1:10 # number of reps
input$generations <- 25 #runtime of the simulation

input$meanFemProgeny <- 300 #average of femlae progeny per queen
input$meanMalProgeny <- 300 #average of male progeny per queen
input$meanFemMatings <- 3.275 #average number of times females mate
input$meanMalMatings <- 0.9 #average number of times males mate
input$maxFemMatings <- 1 #maximum number of times females mate
input$maxMalMatings <- 1 #maximum number of times males mate
input$k <- 1000 #simulated k --> carrying capacity (Carrying capacity K equals to 10.6/km^2)
input$rmax <- 10 #growth rate of the population
input$N <- 1000 #size of start WT population
input$gdSex <- c("F", "M") #which sex carries the gene drive F or M
input$nGD <- 100 #number of gene drive carrying animals to introduce
input$multiplex <- 1 #how many multiplexes in the gene drives, not used currently
input$strategy <- c(1,2,3,4) #what targeting strategy to use 1 = neutral, 2 = male, 3 = female
input$pnhej <- 0.02 #probability of non-homologous end joining, determines the resistance alleles (0.2 in mosquitos)
input$cutRate <- 0.95 #propability CRISPR cuts the opposite DNA strand
input$hEffect <- FALSE #logical, determines whether there is a fitness cost associated with the gene drive
input$pHMort <- 0 #only if hEffect == TRUE, mortality of gene drive carriers.
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

save(modelOutput, file = "Fig2aS1.Rdata")

#########################################
########## Plot plots ###################
#########################################

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
                        hEffect = factor(hEffect),
                        pHMort = factor(pHMort))

haploRep <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Introduction sex" = gdSex, Allele = name, Frequency = value)

haploRep$strategy = factor(haploRep$strategy, 
                           levels = c(1,2,3,4),
                           labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep$`Introduction sex` = factor(haploRep$`Introduction sex`, 
                           levels = c(1,2),
                           labels = c("Female","Male"))

p1 <- ggplot(data = haploRep) +
  facet_grid(
    `Introduction sex` ~ strategy,
    scales = "fixed",
    labeller = labeller(.cols = label_value, .rows = label_both)
  ) +
  geom_line(aes(
    x = generation,
    y = Frequency,
    group = interaction(Allele, repetitions),
    colour = Allele
  )) +
  scale_colour_manual(values=met.brewer("Greek", 4)) +
  xlab("Generation") +
  labs(title = "Asian hornet")
p1

PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=12, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12),
        plot.title = element_text(hjust = 0.5))

ggsave(plot = p1 + PaperTheme, filename = "Fig2aS1.png", height = 12, width = 20, unit = "cm")

fig2as1 <- p1 + PaperTheme
save(fig2as1, file = "Fig2aS1_plot.Rdata")

###################
##### Get SD ######
###################

haploRepSD <- filter(haploRep, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Introduction sex`, Allele) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")
