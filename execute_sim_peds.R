# Set-up
library(tidyverse)
source("./ped_creator_functions.R") # import simulation functions
DOC <- readRDS("./DOC.rds") # death penetrances for death from non-cancer causes

# Generate test families 
# includes all probands/pedigrees and 100% of gene information is present
sim.data <- create.test.fams(num.fams = 100, 
                             proband.carriers = NULL, 
                             path.allele.prevalence = 0.1, 
                             distribution = "exponential",
                             c.cancer.dist.mean = 40,
                             nc.cancer.dist.mean = 80,
                             c.cancer.dist.var = NULL,
                             nc.cancer.dist.var = NULL,
                             doc = DOC)
test.peds <- sim.data$Data
desc <- sim.data$Desc
peds.file.name <- sim.data$DataFile
desc.file.name <- sim.data$DescFile
version <- sim.data$Version

# Pre-MENDEL Run Modifications (Optional)
# optionally remove pedigrees with probands that don't have specified allele
# statuses and/or mask a proportion of gene information from each pedigree
mod.data <- modify.pedigrees(df = test.peds,
                             selected.alleles = c("1/2","2/1","2/2"), 
                             mask.proportion = 0.5,
                             desc.file = desc.file.name)
mod.peds <- mod.data$Data
mod.desc <- mod.data$Desc
mod.peds.file.name <- mod.data$DataFile
mod.desc.file.name <- mod.data$DescFile
mod.version <- mod.data$Version
