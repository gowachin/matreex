# Run all tests
devtools::load_all('.')
devtools::document('.')
devtools::test()
devtools::check()

pkgdown::init_site()
pkgdown::build_home()
pkgdown::build_site()

devtools::build_vignettes()


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

