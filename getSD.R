setwd(/scratch/bell/ymeiborg)
source("./model_function.R")

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

haploRepSD_2a <- filter(haploRep_2a, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Release`, Allele) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")

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

haploRepSD_2b <- filter(haploRep_2b, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Release`, Allele) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")


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

haploRepSD_s1a <- filter(haploRep_s1a, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Release`, Allele) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")

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

haploRepSD_s1b <- filter(haploRep_s1b, strategy == "Neutral" | strategy == "Female infertility", generation == 7) %>%
  group_by(strategy, `Release`, Allele) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")