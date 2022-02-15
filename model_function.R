###########
## Setup ##
###########

rm(list = ls())

packages <- c('AlphaSimR','tidyr', 'dplyr', 'reshape2', 'data.table',
              'doParallel', 'compiler', 'ggplot2', 'extraDistr',
              'patchwork','MetBrewer')

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos='http://cran.us.r-project.org')
      library(x, character.only = TRUE)
    }
  }
)

####################################
## Gene Drive simulation function ##
####################################

modelHornets <- function(input){
  
  ##########################################
  ## ---- Assign all variables in input ----
  ##########################################
  
  for (i in 1:ncol(input)) {assign(names(input)[i], input[,i])}
  
  ######################
  ## Define functions ##
  ######################
  
  source("../AltphaSimR.R")
  
  gdStartPop <- function(founderPop, N, nGD, gdSex, simParam = SP) {
    
    pop <- newPop(founderPop, simParam = simParam)
    pop <- editGenomeFix(pop, ind = 1:pop@nInd, chr = rep(1, pop@nLoci), 
                         segSites = 1:pop@nLoci, allele = 0, 
                         simParam = simParam)
    
    pop_f <- pop[1:(N+nGD)]
    pop_f@sex <- rep("F", pop_f@nInd)
    
    pop_m <- pop[(N+nGD+1):((N+nGD)*2)]
    pop_m@sex <- rep("M", pop_m@nInd)
    
    if (gdSex == "F"){
      pop_f <- editHaplo(pop_f, ind = 1:nGD, chr = rep(1, pop_f@nLoci), 
                          segSites = 1:pop_f@nLoci, allele = 1, 
                          haplotype = rep(1, pop_f@nLoci), simParam = simParam)
    } else if (gdSex == "M") {
      pop_m <- editGenomeFix(pop_m, ind = 1:nGD, chr = rep(1, pop_m@nLoci), 
                              segSites = 1:pop_m@nLoci, allele = 1, 
                              simParam = simParam)
    }
    
    pop <- mergePops(list(pop_f, pop_m))
    
    return(pop)
  }
  
  fertility <- function(pop, strategy = 1){
    # If there are no individuals of either sex, the model stops later on
    if ("M" %in% pop@sex) {
      pop_m <- pop[pop@sex == "M"]
    } else {
      return(pop)}
    if ("F" %in% pop@sex) {
      pop_f <- pop[pop@sex == "F"]
    } else {
      return(pop)}
    
    if (strategy == 1){ # Neutral gene drive
      return(pop)
    }
    
    if (strategy == 2){ # Male infertility
      haploDF <- createHaploDF(pop_m)
      haplo1 <- haploDF[haplo == 1 & (locus1 == "GD" | locus1 == "NF"), id]
      haplo2 <- haploDF[haplo == 2 & (locus1 == "GD" | locus1 == "NF"), id]
      if (length(haplo1) > 0 & length(haplo2) > 0){
        pop_m <- pop_m[!(pop_m@id %in% intersect(haplo1, haplo2))]
        pop <- mergePops(list(pop_m, pop_f))
      }
      
    } else if (strategy == 3) { # Female infertility
      haploDF <- createHaploDF(pop_f)
      haplo1 <- haploDF[haplo == 1 & (locus1 == "GD" | locus1 == "NF"), id]
      haplo2 <- haploDF[haplo == 2 & (locus1 == "GD" | locus1 == "NF"), id]
      if (length(haplo1) > 0 & length(haplo2) > 0){
        pop_f <- pop_f[!(pop_f@id %in% intersect(haplo1, haplo2))]
        pop <- mergePops(list(pop_m, pop_f))
      }
    } else if (strategy == 4) { # Female infertility
      haploDF <- createHaploDF(pop)
      haplo1 <- haploDF[haplo == 1 & (locus1 == "GD" | locus1 == "NF"), id]
      haplo2 <- haploDF[haplo == 2 & (locus1 == "GD" | locus1 == "NF"), id]
      if (length(haplo1) > 0 & length(haplo2) > 0){
        pop <- pop[!(pop@id %in% intersect(haplo1, haplo2))]
      }
    }
    
    return(pop)
  }
  
  heterozygousMortality <- function(pop, pHeteroMort) {
    if (pHeteroMort > 0) { 
      haploDF <- createHaploDF(pop)
      GDs <- haploDF[locus1 == "GD", id]
      if (length(GDs) > 0) {
        pop_GDs <- pop[pop@id %in% unique(GDs)]
        pop <- pop[!(pop@id %in% pop_GDs@id)]
        pop_GDs <- mortality(pop_GDs, mortality = pHeteroMort)
        pop <- mergePops(list(pop, pop_GDs))
      }
    }
    return(pop)
  }
  
  createHaploDF <- function(pop, simParam = SP) {
    haplotype <- setDT(bind_cols(name = row.names(pullQtlHaplo(pop, 
                                                          simParam = simParam)), 
                                 QTL_1 = pullQtlHaplo(pop, 
                                                      simParam = simParam)[,1], 
                                 QTL_2 = pullQtlHaplo(pop, 
                                                      simParam = simParam)[,2]))
    haplotype[, `:=`(id = strsplit(name, split = "_")[[1]][1], 
                     haplo = strsplit(name, split = "_")[[1]][2]), 
              by = name] [, name := NULL]
    
    haplotype[, locus1 := fifelse(test = QTL_1 == 1, 
                                  yes = fifelse(test = QTL_2 == 1, 
                                                yes = "GD", no = "RE"),
                                  no = fifelse(test = QTL_2 == 1, 
                                               yes = "NF", no = "WT"))]
    return(haplotype[ , .(id, haplo, locus1)])
  }
  
  maleDH <- function(pop, meanMaleProgeny, simParam = SP){
    pop <- homing(pop, p_nhej = pnhej, cut_rate = cutRate, simParam = simParam)
    tmp <- list()
    for(ind in 1:pop@nInd){
      tmp <- c(tmp, pop[ind])
    }
    tmp <- lapply(tmp, FUN = function(x) 
      makeDH(x, nDH = actuar::rztpois(n = 1, lambda = meanMaleProgeny), 
             keepParents = FALSE, simParam = simParam))
    pop <- mergePops(tmp)
    pop@sex <- rep("M", pop@nInd)
    return(pop)
  }
  
    maleOffspring <- function(females, meanMaleProgeny, simParam = SP){
      females <- homing(females, p_nhej = pnhej, cut_rate = cutRate, simParam = simParam)
    tmp <- list()
    for(ind in 1:females@nInd){
      tmp <- c(tmp, females[ind])
    }
    tmp <- lapply(tmp, FUN = function(x) 
      makeDH(x, nDH = actuar::rztpois(n = 1, lambda = meanMaleProgeny), 
             keepParents = FALSE, simParam = simParam))
    offspring <- mergePops(tmp)
    offspring@sex <- rep("M", offspring@nInd)
    return(offspring)
  }
  
  polyandryCross <- function(pop,
                             femaleMatings, maleMatings, 
                             maxFem, maxMal, simParam = SP) {
    pop_f <- pop[pop@sex == "F"]
    pop_m <- pop[pop@sex == "M"]
    
    nMates <- rtpois(pop_f@nInd, femaleMatings, 0, maxFem)
    nMatesMale <- rtpois(pop_m@nInd, maleMatings, -Inf, maxMal)
    if (sum(nMates) < 1 | sum(nMatesMale) < 1) {
      crossPlan <- as.data.table(matrix(NA, ncol = 2, nrow = 0))
      return(crossPlan) # Return an empty crossPlan
    }
    
    crossPlan <- as.data.table(matrix(NA, ncol = 2, nrow = sum(nMates)))
    setnames(crossPlan, c("Mothers", "Fathers"))
    
    crossPlan[, Mothers := rep(pop_f@id, nMates)]
    
    maleIDs <- rep(pop_m@id, nMatesMale)
    if (length(maleIDs) < sum(nMates)) {
      crossPlan <- crossPlan[match(sample(Mothers, 
                                          length(maleIDs)),
                                   Mothers),][,Fathers := sample(maleIDs, 
                                                              length(maleIDs))]
    } else {
      crossPlan <- crossPlan[, Fathers := sample(maleIDs, nrow(crossPlan))]
    }

    return(crossPlan)
  }
  
  removeUnsuccessfulCrosses <- function(crosses, fertileFemales, fertileMales){
    # Keep only fertile crosses
    crosses <- crosses[Mothers %in% fertileFemales & Fathers %in% fertileMales]
    
    return(crosses)
  }
  
  femaleOffspring <- function(females, males, crossPlan, meanProgeny, 
                              simParam = SP){
    females <- homing(females, p_nhej = pnhej, cut_rate = cutRate, 
                      simParam = simParam)
    
    offspring <- makeCross2_pois(females = females, 
                                 males = males, crossPlan = as.matrix(crossPlan), 
                                 meanProgeny = meanProgeny, simParam = simParam)
    offspring@sex <- rep("F", offspring@nInd) 
    
    return(offspring)
  }
 
  homing <- function(pop, p_nhej, cut_rate, simParam = SP) {
    haplotype <- as_tibble(createHaploDF(pop))
    
    for (row in seq(1, nrow(haplotype), 2)) {
      for (column in 3:ncol(haplotype)){
        tmpHaplo <- haplotype[row:(row+1), c(1:2, column)]
        if (!any(tmpHaplo[,3] == "GD")) {next}
        qtl <- which(3:ncol(haplotype) == column)*2
        if (any(tmpHaplo[, 3] == "WT") & 
            any(tmpHaplo[, 3] == "GD")){
          if (runif(1) < cut_rate) { # Cutting
            if (runif(1) < 1-p_nhej) { # Homing
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = 1, # 1 because it's the gene drive
                               haplotype = which(tmpHaplo[,3] != "GD"),
                               simParam = simParam)
            } else if (runif(1) < 2/3) { # NHEJ, non-functional repair
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = c(0, 1), # 0,1 for NF
                               haplotype = which(tmpHaplo[,3] != "GD"), 
                               simParam = simParam)
            }  else { # NHEJ, functional repair
              pop <- editHaplo(pop = pop, 
                               ind = which(seq(1, nrow(haplotype), 2) == row), 
                               chr = rep(1, 2),
                               segSites = (qtl-1):qtl,
                               allele = c(1, 0), # 1,0 for RE
                               haplotype = which(tmpHaplo[,3] != "GD"), 
                               simParam = simParam)
            }
          }
        }
      }
    }
    
    return(pop)
  }
  
  mortality <- function(pop, mortality){
    if (mortality == 0){
      return(pop)
    }
    deaths <- which(runif(pop@nInd) < mortality)
    pop <- pop[-deaths]
    pop
  }
  
  ###########################
  ## Initialise population ##
  ###########################
  
  founderPop <- quickHaplo(nInd=((N+nGD)*2), nChr=1, segSites=2, genLen = 0)
  
  SP <- SimParam$new(founderPop)
  SP$addTraitA(nQtlPerChr=2)
  SP$setSexes("yes_sys")

  pop <- gdStartPop(founderPop = founderPop, N = N, nGD = nGD, gdSex = gdSex, simParam = SP)
  crossPlan <- as.data.table(matrix(NA, ncol = 2, nrow = (N+nGD)))
  setnames(crossPlan, c("Mothers", "Fathers"))
  crossPlan[, Mothers := pop@id[1:(N+nGD)] ]
  crossPlan[, Fathers := pop@id[(N+nGD+1):pop@nInd]]
  
  if (gdSex == "M") {
    crossPlan <- crossPlan[((nGD + 1):(nGD + N)), .(Mothers, Fathers)]
    popMaleGD <- pop[pop@sex == "M"][1:nGD]
  }
  
  queens <- pop[pop@id %in% crossPlan[,Mothers]]
  drones <- pop[pop@id %in% crossPlan[,Fathers]]
  
  #################################
  ## Setup data collection table ##
  #################################
  
  generation <- 0:generations
  var <- c('popSizeF', 'WT', 'GD', 'NF', 'RE', 'homoWT', 'homoGD', 'heteroGD')
  results <- array(NA,
                   dim=c(length(generation),length(var)),
                   dimnames=list(generation=generation, var=var))
  N0 <- pop[pop@sex == "F"]@nInd
  results[1, var] <- c(N0,
                       ((N*2)+ifelse(gdSex == "F", nGD, 0))/(N0 * 2), 
                       ifelse(gdSex == "F", nGD/(N0 * 2), 0), 
                       0, 
                       0, 
                       N, 
                       0, 
                       ifelse(gdSex == "F", nGD, 0))
  
  #####################
  ## Start modelling ##
  #####################
  
  for(generation in 1:generations){
    
    cat("\r",sprintf("Generation %i/%i...", generation, generations))
    
    #########################
    ## Population dynamics ##
    #########################
    
    # Save current population numbers
    Nt <- queens@nInd
    
    ##### Offspring generation ##### 
    # Haploid male offspring
    offspring_m <- maleOffspring(females = queens, 
                                 meanMaleProgeny = meanMalProgeny, 
                                 simParam = SP)
    
    # In case of a male GD, adds the GD carrying males to the WT population
    if (generation == 1 & gdSex == "M") {
      offspring_m <- mergePops(list(offspring_m, popMaleGD))
    } 
    
    # Diploid female offspring
    offspring_f <- femaleOffspring(females = queens, 
                                   males = drones, 
                                   crossPlan = crossPlan,
                                   meanProgeny = meanFemProgeny,
                                   simParam = SP)
    pop <- mergePops(list(offspring_m, offspring_f))
    
    ##### Mating ##### 
    # Gene drive fitness cost
    pop <- heterozygousMortality(pop = pop, pHeteroMort = pHMort)
    
    # Makes cross plan based on whole population
    crossPlan <- polyandryCross(pop = pop,
                                femaleMatings = meanFemMatings,
                                maleMatings = meanMalMatings,
                                maxFem = maxFemMatings,
                                maxMal = maxMalMatings,
                                simParam = SP)
    
    if (nrow(crossPlan) == 0) {  
      results[(generation + 1):nrow(results), var[1]] <- rep(0, 
                                                             length((generation + 1):nrow(results)))
      break
    }
    
    # Keep only reproducing individuals
    fertiles <- fertility(pop = pop, strategy = strategy)
    queens <- fertiles[fertiles@sex == "F"]
    drones <- fertiles[fertiles@sex == "M"]
    
    if (queens@nInd == 0 | drones@nInd == 0) {  
      results[(generation + 1):nrow(results), var[1]] <- rep(0, 
                                                             length((generation + 1):nrow(results)))
      break
    }
    
    ##### Mortality ##### 
    # Density dependent female mortality with logistic function
    maxPopSize <- rpois(1, k/(1+((k-Nt)/Nt)*rmax^-1))
    mortalityRate <- 1-maxPopSize/queens@nInd
    
    queens <- mortality(queens, mortalityRate)
    
    if (queens@nInd == 0 | drones@nInd == 0) {
      results[(generation + 1):nrow(results), var[1]] <- rep(0, 
                                               length((generation + 1):nrow(results)))
      break
    }
    
    # Remove infertile crosses
    crossPlan <- removeUnsuccessfulCrosses(crosses = crossPlan,
                                           fertileFemales = queens@id,
                                           fertileMales = drones@id)
    
    if (nrow(crossPlan) == 0) {  
      results[(generation + 1):nrow(results), var[1]] <- rep(0, 
                                                             length((generation + 1):nrow(results)))
      break
    }
    
    #############################
    ## Track population values ##
    #############################
    
    haploDF       <- createHaploDF(queens)
    countHaplo    <- haploDF[, .N, by = locus1]
    haploDFWider  <- setDT(pivot_wider(haploDF, id_cols = id, 
                                       names_from = haplo, 
                                       names_prefix = "haplo", 
                                       values_from = locus1))
    popSizeF      <- queens@nInd
    WT            <- ifelse(length(countHaplo[locus1 == "WT", N]) > 0, 
                            yes = countHaplo[locus1 == "WT", N]/(popSizeF * 2),
                            no = 0)
    GD            <- ifelse(length(countHaplo[locus1 == "GD", N]) > 0, 
                            yes = countHaplo[locus1 == "GD", N]/(popSizeF * 2),
                            no = 0)
    NF            <- ifelse(length(countHaplo[locus1 == "NF", N]) > 0, 
                            yes = countHaplo[locus1 == "NF", N]/(popSizeF * 2),
                            no = 0)
    RE            <- ifelse(length(countHaplo[locus1 == "RE", N]) > 0, 
                            yes = countHaplo[locus1 == "RE", N]/(popSizeF * 2),
                            no = 0)
    homoWT        <- haploDFWider[haplo1 == "WT" & haplo2 == "WT", .N]
    homoGD        <- haploDFWider[haplo1 == "GD" & haplo2 == "GD", .N]
    heteroGD      <- haploDFWider[xor(haplo1 == "GD", haplo2 == "GD"), .N]
    
    
    for (v in var) {
      results[generation+1,v] <- get(x=v)
    }
  }
  
  ###########################
  ## ---- Process output ----
  ###########################
  
  results <- cbind(generation=0:generations, input, results, row.names = NULL) 
  
  return(results)
}

######################
## Compile function ##
######################

modelHornets.comp <- cmpfun(modelHornets)
