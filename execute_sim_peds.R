# Set-up
library(tidyverse)
source("./ped_creator_functions.R") # import simulation functions
DOC <- readRDS("./DOC.rds") # death penetrances for death from non-cancer causes

# Generate test families 
# includes all probands/pedigrees and 100% of gene information is present
sim.data <- create.test.fams(num.fams = 6000, 
                             proband.carriers = NULL, 
                             path.allele.prevalence = 0.2, 
                             distribution = "exponential",
                             c.cancer.dist.mean = 50,
                             nc.cancer.dist.mean = 100,
                             c.cancer.dist.var = NULL,
                             nc.cancer.dist.var = NULL,
                             doc = DOC)
test.peds <- sim.data$Data
desc <- sim.data$Desc
peds.file.name <- sim.data$DataFile
desc.file.name <- sim.data$DescFile
version <- sim.data$Version

