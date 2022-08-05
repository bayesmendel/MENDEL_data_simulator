# This script contains the functions to simulated a specified number of test 
# families for MENDEL option 14. It will save a single .in file with the all 
# pedigrees and a .csv file which contains the parameters used to simulate the 
# test families. Seethe function doc strings for more details.

#' Simulates a new person for a pedigree
#' 
#' Adds a simulated proband or relative to a pedigree data frame with a hypothetical
#' gene status and cancer status based on a pathogenic allele frequency and a cancer
#' probability distribution.
#' 
#' @param ped. data frame of the pedigree. To create a new pedigree with a 
#' single founding father set ped to NULL and ensure founding.father = T.
#' @param motherid. Numeric Mother ID of the new person. Omit if this is a founder or
#' a parent that is related to the family only through their children (ie married-in).
#' @param fatherid. Numeric Father ID of the new person. Omit if this is a founder or 
#' a parent that is related to the family only through their children (ie married-in)
#' @param proband. Logical operator indicating if the person is a proband. Default
#' is FALSE.
#' @param pb.alleles. Numeric vector of length two where values can only be 0 
#' (not pathogenic) or 1 (pathogenic). Only provide a value is proband = TRUE.
#' The first element is the proband's allele 1 pathogenic status and the second
#' element is the proband's allele 2 pathogenic status.
#' @param founding.father. Logical operator indicating if the person to be created is
#' the founding father or not. If the argument is TRUE, a new pedigree is created.
#' Default is FALSE.
#' @param ff.allele. Numeric value representing the pathogenic allele status which
#' thee proband inherited from the founding mother. Must be 0 (not pathogenic)
#' or 1 (pathogenic). Only provide a value if founding.father = TRUE.
#' @param fm.allele. Numeric value representing the pathogenic allele status the
#' proband inherited from the founding mother. Must be 0 (not pathogenic)
#' or 1 (pathogenic). Only provide a value if the simulated person is intended to
#' be the pedigree's founding mother.
#' @param partner. Logical operator indicating if the person to be created is to be a
#' partner of the last person created in the pedigree. If TRUE, partner.sex and
#' partner.age must have values. Default is FALSE.
#' @param partner.sex. String of either c("M","F") representing the sex of partner
#' that already exists in the pedigree. Optional, only required if partner = TRUE.
#' @param partner.age. Numeric value representing the age of the partner that
#' already exists in the pedigree. Optional, only required if partner = TRUE.
#' @param path.allele.prevalence. Prevalence of the pathogenic allele of study.
#' @param distribution. Character string of the parametric distribution type, 
#' one of c("binomial", "exponential", "normal").
#' @param c.cancer.dist.mean. The mean value of the distribution of the cancer
#' penetrance for carriers.
#' @param nc.cancer.dist.mean. The mean value of the distribution of the cancer
#' penetrance for non-carriers.
#' @param c.cancer.dist.var. The variance of the distribution of the cancer
#' penetrance for carriers. Optional for distributions where the variance is 
#' just a transformation of the mean.
#' @param nc.cancer.dist.var. The variance of the distribution of the cancer 
#' penetrance for non-carriers. Optional for distributions where the variance is 
#' just a transformation of the mean.
#' @param doc. data frame that contains death penetrances for death from non-cancer
#' causes with columns c("age","sex","deathPen").
#' @details 
#'  - The founding father's age is randomly selected between 70 and 95 and the 
#'  founding mother is +/- 10 years of their partner's age, with a max age of 95.
#'  - Children are always between 15 years younger than the youngest parent and
#'  40 years younger than their oldest parent. 
#'  - Partners' ages that "marry-in" to the family are +/- 10 years of their 
#'  spouse but cannot be younger than 16.
#'  - Generally, the age rules above result in a pedigree with three levels.
#' @return a pedigree data frame with columns:
#'  - ID: unique numeric identified of the person within the family
#'  - MotherID: identifier of the person's mother corresponding to the ID column
#'  - FatherID: identifier of the person's father corresponding to the ID column
#'  - Sex: person's sex, one of c("M","F") for male or female
#'  - Allele1: numeric value, one of c(0,1) corresponding to absence or presence
#'  of a pathogenic allele.
#'  - Allele2: numeric value, one of c(0,1) corresponding to absence or presence
#'  of a pathogenic allele.
#'  - Gene: numeric value, one of c(0,1) corresponding to absence of presence of
#'  at least one pathogenic allele. This column identifies carriers of the pathogenic
#'  mutation.
#'  - Proband: numeric value, 1 if proband, 0 otherwise.
#'  - Cancer: numeric value, one of c(0,1) corresponding to no affection or 
#'  affection with cancer due to carrier status.
#'  - AgeCancer: numeric value, age of cancer diagnosis.
#'  - Age: numeric value, subject's current age if alive or death age if dead.
#'  - Death: numeric value, 1 if dead, 0 is alive.
create.person <- function(ped, motherid = NULL, fatherid = NULL, 
                          proband = FALSE, pb.alleles = NULL,
                          founding.father = F, ff.allele = NULL, fm.allele = NULL,
                          partner = F, 
                          partner.sex = NULL, partner.age = NULL, 
                          path.allele.prevalence, 
                          distribution,
                          c.cancer.dist.mean,
                          nc.cancer.dist.mean,
                          c.cancer.dist.var = NULL,
                          nc.cancer.dist.var = NULL,
                          doc){
  
  # check input combinations
  if(founding.father & is.null(ff.allele)){
    stop("No founding father allele")
  } 
  if(proband & all(is.na(pb.alleles))){
    stop("No proband alleles")
  } 
  if(sum(proband, founding.father, partner) > 1){
    stop("Cannon specify person as multiple types of either proband, founding father, and/or partner.")
  }
  if(founding.father & (!is.null(motherid) | !is.null(fatherid))){
    stop("founding father should not have a motherid or fatherid specified.")
  }
  
  # ids
  id <- ifelse(founding.father, 1, nrow(ped) + 1)
  if(founding.father | partner){ 
    motherid <- NA
    fatherid <- NA
  }
  
  # proband
  proband <- ifelse(proband == TRUE, 1, 0)
  
  # sex and age limits
  if(founding.father){
    sex <- "M"
    age.max <- 95
    age.min <- 70
  } else if(partner){
    sex <- ifelse(partner.sex == "M", "F", "M")
    age.max <- partner.age + 10
    age.min <- partner.age - 10
  } else {
    sex <- sample(c("F","M"), 1)
    mother.age <- ped$Age[ped$ID == motherid]
    father.age <- ped$Age[ped$ID == fatherid]
    age.max <- ifelse(min(c(mother.age, father.age)) - 15 < 1, 1, 
                      min(c(mother.age, father.age)) - 15)
    age.min <- ifelse(max(c(mother.age, father.age)) - 40 < 1, 1, 
                      max(c(mother.age, father.age)) - 40)
  }
  
  # age
  age <- sample(age.min:age.max, 1)
  if(partner){ 
    age <- ifelse(age < 16, 16, ifelse(age > 95, 95, age)) 
  }
  
  # death and death age
  death.pens <- doc$deathPen[which(doc$sex == sex & doc$age >= age.min & doc$age <= age)]
  death <- ifelse(proband, 0, rbinom(1, 1, prob = sum(death.pens)))
  # make the current age the death age, if the person is dead
  age <- ifelse(death == 1, sample(age.min:age, size = 1, prob = death.pens), age)
  
  ## alleles
  # proband: both provided (from create.test.fam())
  if(proband){
    pb.alleles <- sample(pb.alleles, size = 2, replace = F)
    allele1 <- pb.alleles[1]
    allele2 <- pb.alleles[2]
    
    # founder: 1 allele from the proband, 1 random
  } else if(founding.father | (partner & !is.null(fm.allele))){
    inherited.allele <- ifelse(!is.null(ff.allele), ff.allele, fm.allele)
    allele.position <- sample(c(1,2), size = 1)
    allele1 <- ifelse(allele.position == 1, inherited.allele, 
                      rbinom(1, 1, prob = path.allele.prevalence))
    allele2 <- ifelse(allele.position == 2, inherited.allele, 
                      rbinom(1, 1, prob = path.allele.prevalence))
    
    # partner that "married-in": two random alleles
  } else if(partner){
    allele1 <- rbinom(1, 1, prob = path.allele.prevalence)
    allele2 <- rbinom(1, 1, prob = path.allele.prevalence)
    
    # children of mothers and fathers in the pedigree, excluding the proband,
    # determined by Mendelian inheritance
  } else {
    allele1 <- sample(c(ped$Allele1[ped$ID == motherid], 
                        ped$Allele2[ped$ID == motherid]), size = 1)
    allele2 <- sample(c(ped$Allele1[ped$ID == fatherid], 
                        ped$Allele2[ped$ID == fatherid]), size = 1)
  }
  
  # carrier status
  gene <- ifelse(allele1 == 1 | allele2 == 1, 1, 0)
  
  # cancer and cancer age
  cancer.dist.mean <- ifelse(gene == 1, c.cancer.dist.mean, nc.cancer.dist.mean)
  if(!is.null(c.cancer.dist.var) & !is.null(nc.cancer.dist.var) &
     !distribution %in% c("binomial","negative binomial","poisson","exponential")){
    cancer.dist.var <- ifelse(gene == 1, c.cancer.dist.var, nc.cancer.dist.var)
  }
  if(distribution == "binomial"){
    cancer.prob <- cancer.dist.mean
    cancer.age.prob <- rep(1 / age, age)
  } else if(distribution == "exponential"){
    cancer.prob <- pexp(q = age, rate = 1/cancer.dist.mean)
    cancer.age.prob <- dexp(x = 1:age, rate = 1/cancer.dist.mean)
  } else if(distribution == "normal"){
    cancer.prob <- pnorm(q = age, mean = cancer.dist.mean, sd = sqrt(cancer.dist.var))
    cancer.age.prob <- dnorm(x = age, mean = cancer.dist.mean, sd = sqrt(cancer.dist.var))
  }
  cancer <- sample(x = c(1,0), size = 1, prob = c(cancer.prob, 1-cancer.prob))
  cancer.age <- ifelse(cancer == 1, 
                       sample(x = 1:age, size = 1, prob = cancer.age.prob),
                       NA)
  
  # new row
  row <- data.frame(ID = c(id),
                    MotherID = c(motherid),
                    FatherID = c(fatherid),
                    Sex = c(sex),
                    Allele1 = c(allele1),
                    Allele2 = c(allele2),
                    Gene = c(gene),
                    Proband = c(proband),
                    Cancer = c(cancer),
                    AgeCancer = c(cancer.age),
                    Age = c(age),
                    Death = c(death))
  if(is.null(ped)){
    ped <- row
  } else {
    ped <- rbind(ped, row)
  }
  
  ped
}


#' Simulated a pedigree
#' 
#' Creates a simulated pedigree with a hypothetical pathogenic gene and cancer.
#' Pedigrees are usually three levels, but could possibly be only be two levels. 
#' The proband is always on the 2nd level and can be specified to either have or 
#' not have a pathogenic gene, or randomly assign pathogenic gene status according 
#' to the specified pathogenic allele frequency. Founders can have between one and
#' five children, and children of founders have have between zero and five children 
#' each. 
#' 
#' @param proband.carriers. See `create.test.fams()`.
#' @param path.allele.prevalence. See `create.person()`.
#' @param c.cancer.dist.mean. See `create.person()`.
#' @param nc.cancer.dist.mean. See `create.person()`.
#' @param c.cancer.dist.var. See `create.person()`.
#' @param nc.cancer.dist.var. See `create.person()`.
#' @param doc. See `create.person()`.
#' @return a data frame containing a single test family with columns that
#' match those created by create.person().
create.test.fam <- function(proband.carriers,
                            path.allele.prevalence, 
                            distribution,
                            c.cancer.dist.mean,
                            nc.cancer.dist.mean,
                            c.cancer.dist.var = NULL,
                            nc.cancer.dist.var = NULL,
                            doc){
  
  # initiate pathogenic allele status of proband and founders
  if(!is.null(proband.carriers)){
    
    # proband must be a carrier of at least 1 allele
    if(proband.carriers){
      pb.allele1 <- 1
      pb.allele2 <- rbinom(n = 1, size = 1, prob = path.allele.prevalence)
      
      # proband must be a non-carrier
    } else if(!proband.carriers){
      pb.allele1 <- pb.allele2 <- 0
    }
    
    # proband.carriers are NULL, so allele statuses generated using prevalence
  } else {
    pb.allele1 <- rbinom(n = 1, size = 1, prob = path.allele.prevalence)
    pb.allele2 <- rbinom(n = 1, size = 1, prob = path.allele.prevalence)
  }
  
  father.allele.contrib <- sample(c("PBallele1","PBallele2"), size = 1)
  if(father.allele.contrib == "PBallele1"){
    father.allele <- pb.allele1
    mother.allele <- pb.allele2
  } else if(father.allele.contrib == "PBallele2"){
    father.allele <- pb.allele2
    mother.allele <- pb.allele1
  }
  
  # founding father
  ped <- create.person(ped = NULL, 
                       founding.father = T, ff.allele = father.allele,
                       path.allele.prevalence = path.allele.prevalence, 
                       distribution = distribution,
                       c.cancer.dist.mean = c.cancer.dist.mean, 
                       nc.cancer.dist.mean = nc.cancer.dist.mean, 
                       c.cancer.dist.var = c.cancer.dist.var, 
                       nc.cancer.dist.var = nc.cancer.dist.var,
                       doc = doc)
  
  # founding mother
  ped <- create.person(ped = ped, 
                       fm.allele = mother.allele, partner = T, 
                       partner.sex = "M", 
                       partner.age = ped$Age[1], 
                       path.allele.prevalence = path.allele.prevalence, 
                       distribution = distribution,
                       c.cancer.dist.mean = c.cancer.dist.mean, 
                       nc.cancer.dist.mean = nc.cancer.dist.mean, 
                       c.cancer.dist.var = c.cancer.dist.var, 
                       nc.cancer.dist.var = nc.cancer.dist.var,
                       doc = doc)
  
  # children of founders
  num.child <- sample(1:5, 1, prob = c(0.2, 0.4, 0.2, 0.1, 0.1))
  proband.num <- sample(1:num.child, 1)
  for(nc in 1:num.child){
    pb <- ifelse(nc == proband.num, TRUE, FALSE)
    pb.alleles <- NULL
    if(pb){ pb.alleles <- c(pb.allele1, pb.allele2) }
    ped <- create.person(ped = ped, motherid = 2, fatherid = 1, 
                         proband = pb, pb.alleles = pb.alleles,
                         path.allele.prevalence = path.allele.prevalence, 
                         distribution = distribution,
                         c.cancer.dist.mean = c.cancer.dist.mean, 
                         nc.cancer.dist.mean = nc.cancer.dist.mean, 
                         c.cancer.dist.var = c.cancer.dist.var, 
                         nc.cancer.dist.var = nc.cancer.dist.var,
                         doc = doc)
    
    # no grandchildren if dead before age 16
    if(ped$Death[nrow(ped)] == 1 & ped$Age[nrow(ped)] < 16){ next }
    
    # grandchildren of founders
    num.gchild <- sample(0:5, 1, prob = c(0.1, 0.18, 0.38, 0.18, 0.08, 0.08))
    if(num.gchild > 0){
      
      # create partners of founders' children
      ped <- create.person(ped = ped, partner = T, 
                           partner.sex = ped$Sex[nrow(ped)], 
                           partner.age = ped$Age[nrow(ped)], 
                           path.allele.prevalence = path.allele.prevalence, 
                           distribution = distribution,
                           c.cancer.dist.mean = c.cancer.dist.mean, 
                           nc.cancer.dist.mean = nc.cancer.dist.mean, 
                           c.cancer.dist.var = c.cancer.dist.var, 
                           nc.cancer.dist.var = nc.cancer.dist.var,
                           doc = doc)
      
      # record mother and father ids for creating grandchildren
      if(ped$Sex[nrow(ped)] == "F"){
        m.id <- ped$ID[nrow(ped)]
        f.id <- ped$ID[(nrow(ped)-1)]
      } else {
        m.id <- ped$ID[(nrow(ped)-1)]
        f.id <- ped$ID[nrow(ped)]
      }
      
      # populate grandchildren of founders
      for(ngc in 1:num.gchild){
        ped <- create.person(ped = ped, motherid = m.id, fatherid = f.id, 
                             path.allele.prevalence = path.allele.prevalence, 
                             distribution = distribution,
                             c.cancer.dist.mean = c.cancer.dist.mean, 
                             nc.cancer.dist.mean = nc.cancer.dist.mean, 
                             c.cancer.dist.var = c.cancer.dist.var, 
                             nc.cancer.dist.var = nc.cancer.dist.var,
                             doc = doc)
      }
    }
  }
  
  ped
}


#' Create Test Pedigrees
#' 
#' Creates a specified number of test pedigrees with a hypothetical pathogenic
#' gene and cancer, formatted for MENDEL.
#' 
#' @param num.fams numeric value, the number of test families to create
#' @param proband.carriers logical operator indicating if all test families
#' should have a proband with a pathogenic mutation or not. If NULL, then families
#' will be generated with probands with and without pathogenic mutations based on
#' the allele frequency provided to create.
#' @param path.allele.prevalence. See `create.person()`.
#' @param distribution. See `create.person()`.
#' @param c.cancer.dist.mean. See `create.person()`.
#' @param nc.cancer.dist.mean. See `create.person()`.
#' @param c.cancer.dist.var. See `create.person()`.
#' @param nc.cancer.dist.var. See `create.person()`.
#' @param doc. See `create.person()`.
#' @returns a list with elements:
#'  - `Data`: a data frame of pedigrees with the same columns as created by
#'  `create.person()` however `PedigreeID` and `Twins` has been added, both allele
#'  columns have been removed, and the contents of `Gene` have been reformatted to
#'  contain the allele status information for both alleles. `PedigreeID` contains
#'  a unique numeric pedigree identifier, `Twins` is currently populated with all
#'  blanks, and `Gene` contains strings of the format X1/X2 where X1 is the 1st
#'  allele and X2 is the 2nd allele. Alleles are 1 for non-pathogenic and 2 for
#'  pathogenic. All NAs have been replaced with blanks.
#'  - `Desc`: a one row data frame of descriptors of how the pedigrees were created 
#'  with columns:
#'   - `Version`: a string of the format "YYYY-MM-DD_X" where X is a version number
#'   unique to the version date.
#'   - `NumPeds`: number of pedigrees in the data.
#'   - `PbAlleleFilter`: A string containing the allele statuses for all of the 
#'   probands.
#'   - `PropAlleleMasking`: a number between 0 and 1 indicating the proportion of
#'   gene information masked.
#'   - `Distribution`: the cancer distribution type used.
#'   - `CarrierMean`: the cancer distribution mean for carriers.
#'   - `NoncarrierMean`: the cancer distribution mean for non-carriers.
#'   - `CarrierVar`: the cancer distribution variance for carriers.
#'   - `NoncarrierVar`: the cancer distribution variance for non-carriers.
#'  - `DataFile`: a character string of the path and file name of the pedigree data 
#'  frame.
#'  - `DescFile`: a character string of the path and file name of the descriptors 
#'  data frame.
#'  - `Version`: a character string containing the date and a number that indicates
#'  the file version used in the file's name.
create.test.fams <- function(num.fams, proband.carriers = NULL, 
                             path.allele.prevalence, 
                             distribution,
                             c.cancer.dist.mean,
                             nc.cancer.dist.mean,
                             c.cancer.dist.var = NULL,
                             nc.cancer.dist.var = NULL,
                             doc){
  
  all.peds <- NULL
  for(tf in 1:num.fams){
    tmp.ped <- create.test.fam(proband.carriers = proband.carriers,
                               path.allele.prevalence = path.allele.prevalence, 
                               distribution = distribution,
                               c.cancer.dist.mean = c.cancer.dist.mean,
                               nc.cancer.dist.mean = nc.cancer.dist.mean,
                               c.cancer.dist.var = c.cancer.dist.var,
                               nc.cancer.dist.var = nc.cancer.dist.var,
                               doc = doc)
    tmp.ped <- mutate(tmp.ped, PedigreeID = tf, .before = "ID")
    if(is.null(all.peds)){
      all.peds <- tmp.ped
    } else {
      all.peds <- rbind(all.peds, tmp.ped)
    }
  }
  
  # convert to MENDEL pedigree .in file format
  all.peds <- 
    all.peds %>%
    mutate(Twins = NA, .after = "Sex") %>%
    mutate(across(.cols = starts_with("Allele"), ~ .+1)) %>%
    mutate(Gene = paste0(Allele1,"/",Allele2)) %>%
    select(-c(Allele1,Allele2)) %>%
    mutate(across(.cols = c(MotherID, FatherID, Twins, AgeCancer), 
                  ~ replace_na(as.character(.), "")))
  
  # create descriptors table for the simulated data
  pb.statuses <- ifelse(is.null(proband.carriers), "1/1, 1/2, 2/1, 2/2", 
                        ifelse(proband.carriers, "1/2, 2/1, 2/2", "1/1"))
  f.names <- make.file.names()
  if(distribution == "binomial"){
    c.cancer.dist.var <- c.cancer.dist.mean * (1-c.cancer.dist.mean)
    nc.cancer.dist.var <- nc.cancer.dist.mean * (1-nc.cancer.dist.mean)
  } else if(distribution == "exponential"){
    c.cancer.dist.var <- c.cancer.dist.mean^2
    nc.cancer.dist.var <- nc.cancer.dist.mean^2
  }
  desc <- data.frame(Version = f.names$Version,
                     NumPeds = num.fams,
                     PbAlleleFilter = "1/1, 1/2, 2/1, 2/2",
                     PropAlleleMasking = 0,
                     PathAllelePrevalence = path.allele.prevalence,
                     Distribution = distribution,
                     CarrierMean = c.cancer.dist.mean,
                     NoncarrierMean = nc.cancer.dist.mean,
                     CarrierVar = c.cancer.dist.var,
                     NoncarrierVar = nc.cancer.dist.var)
  
  # save pedigree as .in file and descriptors file as .csv
  
  write.table(all.peds, file = f.names$Pedigree, quote = FALSE, sep = ",", 
              row.names = FALSE, col.names = FALSE)
  write.csv(desc, file = f.names$Desc, row.names = F)
  
  list(Data = all.peds, 
       Desc = desc, 
       DataFile = f.names$Pedigree, 
       DescFile = f.names$Desc,
       Version = f.names$Version)
}


#' Create Pedigree Data Set File Name with Path
#' 
#' Creates a standardized but unique file name (and path) based on the date plus 
#' an appended number.
#' 
#' @details 
#'  - The file path directory is always "./simulated-data"
#'  - The file naming convention is "test_fam_[YYYY-MM-DD]_[VERSION]". Where
#'  the date is from `Sys.Date()` and VERSION is 1 if there are no other files
#'  in the directory from today's date, or, if other files are dated today, then
#'  one greater than the highest version number in the directory.
#' @returns two file names: 
#'  - `Pedigree`: pedigree file name as a .in file
#'  - `Desc`: decriptors file name as a .csv file
make.file.names <- function(){
  
  other.files <- list.files("./simulated-data", 
                            pattern = paste0("^test_fams_", Sys.Date(),"_"), 
                            recursive = T)
  if(length(other.files) > 0){
    last_ <- sapply(str_locate_all(other.files, pattern = "_"), 
                    function(x){ as.numeric(x[nrow(x),"start"])})
    last. <- sapply(str_locate_all(other.files, pattern = "\\."), 
                    function(x){ as.numeric(x[nrow(x),"start"])})
    highest <- max(as.numeric(str_sub(other.files, start = last_+1, end = last.-1)))
    v <- paste0(Sys.Date(), "_", highest+1)
  } else {
    v <- paste0(Sys.Date(), "_1")
  }
  p.name <- paste0("./simulated-data/test_fams_", v, ".in")
  d.name <- paste0("./simulated-data/test_fams_", v, ".csv")
  
  list(Pedigree = p.name, Desc = d.name, Version = v)
}


#' Modify Simulated Pedigrees
#' 
#' Filters pedigrees by probands with specified allele statuses and masks
#' specified proportion of gene information in each pedigree. The function can 
#' also randomly sample a subset of pedigrees to include. Saves the updated
#' pedigrees and descriptors file by a naming convention.
#' 
#' @param df data frame of pedigree data with at least columns `PedigreeID`, 
#' `Proband`, and `Gene`.
#' @param selected.alleles a character vector specifying which the allele status(es)
#' of which to retain. Any of c("1/1", "1/2", "2/1", "2/2") which correspond to
#' no pathogenic alleles for 1/1, heterozygous pathogenic alleles status for either
#' 1/2 or 2/1, and homozygous pathogenic allele status for 2/2. Default is for 
#' all four options with would apply no filtering.
#' @param mask.proportion a numeric value between 0 and 1 indicating the proportion
#' of family members that will have masked gene information. Default is `0` which
#' would provide no masking.
#' @param sample.peds a numeric value of the number of pedigrees to retain. 
#' Default is `NULL` which results in all pedigrees being retained after filtering
#' by alleles.
#' @param desc.file a character string of the path and file name where the 
#' descriptor file for the `df` is saved.
#' @returns a list with elements:
#'  - `Data`: a data frame of modified pedigrees. See `create.test.fams()` for a 
#'  more complete description.
#'  - `Desc`: a data frame of descriptors of how the pedigrees were created/modified.
#'  See `create.test.fams()` for a more complete description.
#'  - `DataFile`: a character string of the path and file name of the pedigree data 
#'  frame.
#'  - `DescFile`: a character string of the path and file name of the descriptors 
#'  data frame.
#'  - `Version`: a character string containing the date and a number that indicates
#'  the file version used in the file's name.
modify.pedigrees <- function(df, 
                             selected.alleles = c("1/1","1/2","2/1","2/2"), 
                             mask.proportion = 0,
                             sample.peds = NULL,
                             desc.file){
  
  # filter by proband gene status
  selected.peds <-
    df %>%
    filter(Proband == 1) %>%
    mutate(Selected = ifelse(Gene %in% selected.alleles, "Keep", "Drop")) %>%
    filter(Selected == "Keep") %>%
    select(PedigreeID) %>%
    as.matrix() %>%
    as.vector()
  df2 <- filter(df, PedigreeID %in% selected.peds)
  
  # mask gene information
  ped.ids <- unique(df2$PedigreeID)
  for(pid in ped.ids){
    tmp.ped <- filter(df2, PedigreeID == pid & Proband == 0)
    mask.cnt <- floor(mask.proportion * nrow(tmp.ped))
    keep.ids <- sample(tmp.ped$ID, size = mask.cnt)
    df2 <- mutate(df2, 
                  Gene = ifelse(PedigreeID == pid & !ID %in% keep.ids, "", Gene))
  }
  
  # randomly sample specified number of pedigrees
  if(!is.null(sample.peds)){
    if(sample.peds > length(ped.ids)){ 
      stop("There are not enough pedigrees remaining to sample.")
    } 
    sample.peds <- sample(ped.ids, size = sample.peds)
    df2 <- filter(df2, PedigreeID %in% sample.peds)
  }
  
  # modify descriptors file
  desc <- read.csv(desc.file)
  fnames <- make.file.names()
  desc$Version <- fnames$Version
  desc$NumPeds = length(unique(df2$PedigreeID))
  desc$PbAlleleFilter <- paste0(selected.alleles, collapse = ", ")
  desc$PropAlleleMasking <- mask.proportion
  
  # save modified files
  write.table(df2, file = fnames$Pedigree, quote = FALSE, sep = ",", 
              row.names = FALSE, col.names = FALSE)
  write.csv(desc, file = fnames$Desc, row.names = F)
  
  list(Data = df2, 
       Desc = desc, 
       DataFile = fnames$Pedigree, 
       DescFile = fnames$Desc,
       Version = fnames$Version)
}


