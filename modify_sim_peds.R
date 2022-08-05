# Set-up
library(tidyverse)
source("./ped_creator_functions.R") # import simulation functions
# change the file argument below as needed
test.peds <- read.table(file = "./simulated-data/test_fams_2022-08-05_1.in",
                        sep = ",",
                        col.names = c("PedigreeID", "ID", "MotherID", "FatherID", 
                                      "Sex", "Twins", "Gene", "Proband",
                                      "Cancer", "AgeCancer", "Age", "Death"))

# Pre-MENDEL Run Modifications (Optional)
# optionally remove pedigrees with probands that don't have specified allele
# statuses and/or mask a proportion of gene information from each pedigree and/or
# sample a subset of pedigrees.
mod.data <- modify.pedigrees(df = test.peds,
                             selected.alleles = c("1/2","2/1","2/2"), 
                             mask.proportion = 0.5,
                             sample.peds = NULL,
                             desc.file = desc.file.name)
mod.peds <- mod.data$Data
mod.desc <- mod.data$Desc
mod.peds.file.name <- mod.data$DataFile
mod.desc.file.name <- mod.data$DescFile
mod.version <- mod.data$Version
