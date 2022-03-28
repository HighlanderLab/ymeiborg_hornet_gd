setwd("/scratch/bell/ymeiborg")
source("./model_function.R")

###############################
############ get data #########
###############################

load("fig2a/Fig2a.Rdata")

haploRep_2a <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep_2a$strategy = factor(haploRep_2a$strategy, 
                           levels = c(1,2,3,4),
                           labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep_2a$`Release` = factor(haploRep_2a$`Release`, 
                            levels = c(1,2),
                            labels = c("Females","Males"))

load("fig2b/Fig2b.Rdata")

haploRep_2b <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep_2b$strategy = factor(haploRep_2b$strategy, 
                              levels = c(1,2,3,4),
                              labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep_2b$`Release` = factor(haploRep_2b$`Release`, 
                               levels = c(1,2),
                               labels = c("Females","Males"))

load("figs1a/FigS1a.Rdata")

haploRep_s1a <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep_s1a$strategy = factor(haploRep_s1a$strategy, 
                              levels = c(1,2,3,4),
                              labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep_s1a$`Release` = factor(haploRep_s1a$`Release`, 
                               levels = c(1,2),
                               labels = c("Females","Males"))


load("figs1b/FigS1b.Rdata")

haploRep_s1b <- as_tibble(modelOutput) %>% 
  select(strategy, gdSex, generation, repetitions, WT:RE, N, nGD, k) %>% 
  pivot_longer(cols = WT:RE) %>% 
  rename("Release" = gdSex, Allele = name, Frequency = value)

haploRep_s1b$strategy = factor(haploRep_s1b$strategy, 
                               levels = c(1,2,3,4),
                               labels = c("Neutral","Male infertility","Female infertility", "Both-sex infertility"))
haploRep_s1b$`Release` = factor(haploRep_s1b$`Release`, 
                                levels = c(1,2),
                                labels = c("Females","Males"))

haploRep <- bind_rows(haploRep_2a, haploRep_2b, haploRep_s1a, haploRep_s1b) %>%
  mutate(dataset = rep(c("2a", "2b", "s1a", "s1b"), each = nrow(haploRep_2a)))

##################################################
############## SD at generation 7 ################
##################################################

haploRepSD_g7 <- filter(haploRep, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Release`, Allele, dataset) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")

##################################################
############## SD at generation 10 ###############
##################################################

haploRepSD_g10 <- filter(haploRep, strategy == "Neutral" | strategy == "Female infertility", generation == 10) %>%
  group_by(strategy, `Release`, Allele, dataset) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")

#####################################
############## max SD ###############
#####################################

haploRepSD_max <- filter(haploRep, strategy == "Neutral" | strategy == "Female infertility") %>%
  group_by(strategy, `Release`, Allele, dataset, generation) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(sd == max(sd), Allele == "GD")