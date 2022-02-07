editGenomeFix <- function (pop, ind, chr, segSites, allele, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr) == length(segSites))
  allele = as.integer(allele)
  stopifnot(all(allele == 0L | allele == 1L))
  allele = as.raw(allele)
  if(length(allele) == 1L){
    allele = rep(allele, length(segSites))
  }
  stopifnot(length(allele) == length(segSites))
  for (selChr in unique(chr)) {
    sel = which(chr == selChr)
    for (i in sel) {
      BYTE = (segSites[i] - 1L)%/%8L + 1L
      BIT = (segSites[i] - 1L)%%8L + 1L
      for (selInd in ind) {
        for (j in 1:pop@ploidy) {
          TMP = pop@geno[[selChr]][BYTE, j, selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = allele[i]
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE, j, selInd] = TMP
        }
      }
    }
  }
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop = pop, simParam = simParam)
  pop@pheno = PHENO
  pop@ebv = EBV
  return(pop)
}

editHaplo <- function (pop, ind, chr, segSites, allele, haplotype = 1:pop@ploidy ,simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  ind = unique(as.integer(ind))
  stopifnot(all(ind %in% (1:pop@nInd)))
  chr = as.integer(chr)
  segSites = as.integer(segSites)
  stopifnot(length(chr) == length(segSites))
  allele = as.integer(allele)
  stopifnot(all(allele == 0L | allele == 1L))
  allele = as.raw(allele)
  if(length(allele) == 1L){
    allele = rep(allele, length(segSites))
  }
  stopifnot(length(allele) == length(segSites))
  for (selChr in unique(chr)) {
    sel = which(chr == selChr)
    for (i in sel) {
      BYTE = (segSites[i] - 1L)%/%8L + 1L
      BIT = (segSites[i] - 1L)%%8L + 1L
      for (selInd in ind) {
        for (j in haplotype) {
          TMP = pop@geno[[selChr]][BYTE, j, selInd]
          TMP = rawToBits(TMP)
          TMP[BIT] = allele[i]
          TMP = packBits(TMP)
          pop@geno[[selChr]][BYTE, j, selInd] = TMP
        }
      }
    }
  }
  PHENO = pop@pheno
  EBV = pop@ebv
  pop = resetPop(pop = pop, simParam = simParam)
  pop@pheno = PHENO
  pop@ebv = EBV
  return(pop)
}

makeCross2_pois <- function (females, males, crossPlan, meanProgeny = 1, simParam = NULL) 
{
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if ((females@ploidy%%2L != 0L) | (males@ploidy%%2L != 0L)) {
    stop("You can not cross aneuploids")
  }
  if (is.character(crossPlan)) {
    crossPlan = cbind(match(crossPlan[, 1], females@id), 
                      match(crossPlan[, 2], males@id))
    if (any(is.na(crossPlan))) {
      stop("Failed to match supplied IDs")
    }
  }
  if ((max(crossPlan[, 1]) > nInd(females)) | (max(crossPlan[, 2]) > nInd(males)) | (min(crossPlan) < 1L)) {
    stop("Invalid crossPlan")
  }
  nProgeny <- rpois(n = length(unique(crossPlan[, 1])), lambda = meanProgeny)
  if (sum(nProgeny) == 0) {
    return(females[-1:-females@nInd]) #return an empty population
  } else {
    tmpCrossPlan <- matrix(ncol = 2, nrow = 0)
    for (mother in unique(crossPlan[, 1])){
      tmpMother <- rep(mother, nProgeny[which(unique(crossPlan[, 1]) == mother)]) #generate the number of offspring per mother
      fathers <- crossPlan[crossPlan[,1] == mother, 2]
      if (length(fathers) == 1){
        tmpFather <- rep(fathers, length(tmpMother))
      } else {
        tmpFather <- sample(x = fathers, size = length(tmpMother), replace = TRUE) #distribute fathers randomly over the number of offspring
      }
      tmpCrossPlan <- rbind(tmpCrossPlan, cbind(tmpMother, tmpFather))
    }
    crossPlan <- tmpCrossPlan
  }
  tmp = AlphaSimR:::cross(females@geno, crossPlan[, 1], males@geno, crossPlan[,
      2], simParam$femaleMap, simParam$maleMap, simParam$isTrackRec,
      females@ploidy, males@ploidy, simParam$v, simParam$p,
      simParam$femaleCentromere, simParam$maleCentromere, simParam$quadProb,
      simParam$nThreads)
  rPop = new("RawPop", nInd = nrow(crossPlan), nChr = females@nChr, 
             ploidy = as.integer((females@ploidy + males@ploidy)/2), 
             nLoci = females@nLoci, geno = tmp$geno)
  if (simParam$isTrackRec) {
    simParam$addToRec(tmp$recHist)
  }
  return(newPop(rawPop = rPop, mother = females@id[crossPlan[,
        1]], father = males@id[crossPlan[, 2]], simParam = simParam,
        iMother = females@iid[crossPlan[, 1]], iFather = males@iid[crossPlan[,
            2]], femaleParentPop = females, maleParentPop = males,
        hist = hist))
}
