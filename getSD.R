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
                           labels = c("Neutral","Male infertility",
                                      "Female infertility", 
                                      "Both-sex infertility"))
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
                              labels = c("Neutral","Male infertility",
                                         "Female infertility", 
                                         "Both-sex infertility"))
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
                              labels = c("Neutral","Male infertility",
                                         "Female infertility", 
                                         "Both-sex infertility"))
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
                               labels = c("Neutral",
                                          "Male infertility",
                                          "Female infertility", 
                                          "Both-sex infertility"))
haploRep_s1b$`Release` = factor(haploRep_s1b$`Release`, 
                                levels = c(1,2),
                                labels = c("Females","Males"))

haploRep <- bind_rows(haploRep_2a, haploRep_2b, haploRep_s1a, haploRep_s1b) %>%
  mutate(dataset = rep(c("2a", "2b", "s1a", "s1b"), each = nrow(haploRep_2a)))

##################################################
############## SD at generation 7 ################
##################################################

haploRepSD_g7 <- filter(haploRep, 
                        strategy == "Neutral" | strategy == "Female infertility", 
                        generation == 7) %>%
  group_by(strategy, `Release`, Allele, dataset) %>%
  summarise(mean = mean(Frequency), sd = sd(Frequency)) %>%
  filter(Allele == "GD")

##################################################
############## SD at generation 10 ###############
##################################################

haploRepSD_g10 <- filter(haploRep, 
                         strategy == "Neutral" | strategy == "Female infertility", 
                         generation == 10) %>%
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

#######################################################
############## generation spread at 0.5 ###############
#######################################################
PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        title=element_text(size=14, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        legend.box = "vertical",
        axis.title=element_text(size=12))

haploRepSD25 <- filter(haploRep, strategy == "Neutral" | strategy == "Female infertility") %>%
  group_by(strategy, `Release`, Allele, dataset, repetitions) %>%
  filter(Allele == "GD" & Frequency >= 0.25) %>%
  #filter(Frequency == min(Frequency)) %>%
  filter(generation == min(generation)) %>%
  group_by(strategy, `Release`, Allele, dataset) %>%
  #summarise(minGen = min(generation), maxGen = max(generation))
  
p <- ggplot(haploRepSD25) +
  facet_grid(Release ~ strategy) +
  geom_violin(aes(x = dataset, y = generation), draw_quantiles = 0.5) +
  geom_jitter(aes(x = dataset, y = generation), width = 0.1) +
  PaperTheme
p
  
  
  