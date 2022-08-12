# Run all tests
rm(list = ls())
devtools::document('.')
devtools::load_all('.')
devtools::test()
devtools::check()

# library(covr, testthat)
# get a shiny to find which line is not yet tested. Very helpfull
report(x = package_coverage(line_exclusions = list("R/Sim_NonDem.R")))
report(x = package_coverage())
# covr::package_coverage()

# rerun all coverage for dataframe analysis.
# cov <- package_coverage()
# view results as a data.frame
# df_cov <- as.data.frame(cov)

# df_cov[df_cov$functions=='tank',]


x <- file_coverage('R/Sim_Deter.R', 'tests/testthat/test-Sim_Deter.R')
x <- file_coverage('R/class_ipm.R', 'tests/testthat/test-class_ipm.R')
x <- file_coverage('R/class_species.R', 'tests/testthat/test-class_species.R')
x <- file_coverage('R/class_forest.R', 'tests/testthat/test-class_forest.R')
x <- file_coverage('R/generic_methods.R', 'tests/testthat/test-gen_meth.R')
report(x)
# zero_coverage() shows only uncovered lines.
# If run within RStudio, `zero_coverage()` will open a marker pane with the
# uncovered lines.
# zero_coverage(cov)


# Update the doc part of the github with pkgdown
# pkgdown::build_site()

# Link between functions
library(mvbutils)
foodweb(where = asNamespace( "treeforce"), cex = 0.8, color.lines = F)

# Update the coverage to codecov, don't forget the token
# covr::codecov()
