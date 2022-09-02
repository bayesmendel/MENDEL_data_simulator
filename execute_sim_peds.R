# Set-up
library(tidyverse)
source("./ped_creator_functions.R") # import simulation functions
DOC <- readRDS("./DOC.rds") # death penetrances for death from non-cancer causes

# Generate test families 
# includes all probands/pedigrees and 100% of gene information is present
sim.data <- create.test.fams(num.fams = 6000, 
                             proband.carriers = NULL, 
                             path.allele.prevalence = 0.2, 
                             distribution = "normal",
                             c.cancer.dist.mean = 35,
                             nc.cancer.dist.mean = 70,
                             c.cancer.dist.var = 144,
                             nc.cancer.dist.var = 400,
                             doc = DOC)
test.peds <- sim.data$Data
desc <- sim.data$Desc
peds.file.name <- sim.data$DataFile
desc.file.name <- sim.data$DescFile
version <- sim.data$Version