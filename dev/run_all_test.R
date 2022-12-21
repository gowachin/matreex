# Run all tests
rm(list = ls())
devtools::load_all('.')
devtools::document('.')
devtools::test()
devtools::check()


# change name
# grep "treeforce" -r ./*

library(treeforce)
?treeforce

# library(covr, testthat)
# get a shiny to find which line is not yet tested. Very helpfull
# report(x = package_coverage(line_exclusions = list("R/Sim_NonDem.R")))
x <- package_coverage(
    line_exclusions = list(
        "R/harvest.R" = 225,
        "R/class_species.R" = 263,
        "R/dev_in.R"
    )
)
report(x)

# view results as a data.frame
# df_cov <- as.data.frame(cov)

# df_cov[df_cov$functions=='tank',]

x <- file_coverage('R/class_fit_sgr.R', 'tests/testthat/test-class_fit_sgr.R')
x <- file_coverage('R/class_ipm.R', 'tests/testthat/test-class_ipm.R')
x <- file_coverage('R/class_mu_gr.R', 'tests/testthat/test-class_mu_gr.R')
x <- file_coverage('R/class_species.R', 'tests/testthat/test-class_species.R',
                   line_exclusions = list("R/class_species.R" = 263))
x <- file_coverage('R/class_forest.R', 'tests/testthat/test-class_forest.R')
x <- file_coverage('R/generic_methods.R', 'tests/testthat/test-gen_meth.R')
x <- file_coverage('R/load_oldIPM_rec.R', 'tests/testthat/test-oldIPM_rec.R')
x <- file_coverage('R/harvest.R', 'tests/testthat/test-harvest.R',
                   line_exclusions = list("R/harvest.R" = 225))
x <- file_coverage('R/make_ipm.R', 'tests/testthat/test-make_ipm.R')
x <- file_coverage('R/step_IPM.R', 'tests/testthat/test-step_IPM.R')
x <- file_coverage('R/Sim_Deter.R', 'tests/testthat/test-Sim_Deter.R')
x <- file_coverage('R/formulas.R', 'tests/testthat/test-formula.R')
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
