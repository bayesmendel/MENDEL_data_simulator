# MENDEL_data_simulator
Simulates test pedigrees for MENDEL penetrance estimation.

**Author**: Stephen Knapp (sknapp@ds.dfci.harvard.edu)

Use `execute_sim_peds.R` to generate simulated pedigrees according to various parameters explained in the `create.test.fams()`, `create.test.fam()`, and `create.person()` doc strings.

Use `modify_sim_peds.R` to modify the generated data by any of the following methods: filter pedigrees by proband allele status, mask a portion of allele information, and/or sample a subset of pedigrees. See the `modify.pedigrees()` doc string.

Both `execute_sim_peds.R` and `modify_sim_peds.R` will save their output data to the `/simulated-data` sub-directory as a .in and will also save a .csv file that contains the simulation parameters. See `make.file.names()` for the file naming convention.

`output_evaluator_functions.R` contains a function called `test.MENDEL.output()` which compares the true versus MENDEL estimated distribution information. Import this function into another .R or .Rmd file to evaluate your experiment's results.

Note `IR.SEER.RData`, `DOC.R`, and `DOC.rds` are all realted to death penetrances and should not need to be modified.

View the detailed doc strings in `ped_creator_functions.R` and `output_evaluator_functions.R` for more information.
