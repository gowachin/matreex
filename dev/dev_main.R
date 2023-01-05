#' Script for development procedure.
#'
#' Developping matreex as a package depend on {devtools} package.
#' It also include gitlab or github CI and {pkgdown} package
#' for documentation presentation.
#'
#' Below are some of the procedure to run once an edit has been
#' done in the package.
#'
#' Most of them need to be used before git commit.
#'
#' Run this script interactively and only select some chunks !
#'
#' maxime Jaunatre 31/03/2022

library(usethis)
library(devtools)
library(pkgdown)
library(fs)
library(beepr)

## Building Binary package for Releases ####
# FROM SO https://stackoverflow.com/questions/54634056/how-to-include-an-html-vignette-in-a-binary-r-package
# build_vignettes_to_inst <- function() {
#     build_vignettes() # Builds vignettes to 'doc' and 'Meta'. Updates '.gitignore'.
#     cat("Builded\n")
#     unlink(c("inst/doc", "inst/Meta"), recursive = TRUE) # Remove the directories if they exist
#     dir.create("inst/doc"); dir.create("inst/Meta") # Create empty directories
#     has_worked <- c( # Copy files to 'inst' subfolders
#         file.copy(list.files("doc", full.names = TRUE), to = "inst/doc"),
#         file.copy(list.files("Meta", full.names = TRUE), to = "inst/Meta")
#     )
#     cat("Moved\n")
#     unlink(c("doc", "Meta"), recursive = TRUE) # Optional: Remove unwanted directories
#     return(all(has_worked)) # Returns TRUE if everything worked OK
# }

# Command for pkg compil ####

## Reload matreex in empty env. ####
devtools::unload("matreex")
.rs.restartR()
devtools::load_all()

## Force C++ compil without edition ####
devtools::unload("matreex")
.rs.restartR()
devtools::clean_dll()
devtools::load_all()

## Reinstall package ####
devtools::document()
devtools::unload("matreex")
.rs.restartR()
remove.packages("matreex")
devtools::install(upgrade = "never",
                  build_vignettes = TRUE)

## Check the full package ####
devtools::check() # takes few minutes

## Compile documentation site ####
devtools::load_all()
rmarkdown::render("vignettes/multisp_deter_sim.Rmd")
rmarkdown::render("vignettes/Basic_functions.Rmd")
rmarkdown::render("README.rmd", clean = TRUE, quiet = TRUE)
fs::file_delete("README.html")
pkgdown::build_site(preview = TRUE, devel = TRUE) # will open website
pkgdown::build_home(preview = TRUE) # will open website

## Clear objects in env ####
rm(list = ls())
gc()
# Restart R
.rs.restartR()

## Building Binary package for Releases ####
build_vignettes_to_inst() # Call the function
devtools::build(binary = TRUE, args = c('--preclean'))
beepr::beep(5)
# we can push releases to gitlab with curl. See this SO post:
# https://stackoverflow.com/questions/29013457/how-to-store-releases-binaries-in-gitlab

# ## Compile the full notice ####
# build_vignettes_to_inst() # This will recompile all vignettes to pdf and move them int/doc/
# setwd("inst/notice")
# tinytex::pdflatex("Notice_matreex.tex")
# tinytex::pdflatex("Notice_matreex.tex")
# setwd("../../")

## Pre compile intensive compu vignettes ##
# sourcing this script will kill R, but running it interactively work.
# also, very long script to run.
# source("vignettes/pre_compile.R")

#' I tested using some options but it's poorly used in the package.
# options(dev = TRUE)
# options(matreex.dev = TRUE) # TRUE show verbose, FALSE not. NULL equal FALSE

# Index of files ####


# Notes ####
#' I can't put simulation in data/ because maximum recommended is 5MB
