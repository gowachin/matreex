test_that("new_ipm works", {

    place <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                  "testthat", "testdata")
    file <- here(place, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- c(sgdd = 2000, wai = -0.2)
    clim_lab <- "paradis"
    compress <- TRUE
    delay <- 0

    expect_identical(new_ipm(IPM, BA, mesh, species, climatic, clim_lab, compress),
                     structure(list(
                         IPM = IPM, BA = BA, mesh = mesh, climatic = climatic,
                         info = c(species = "Yggdrasil", clim_lab = "paradis",
                                  compress = "TRUE", delay = "0")), class = "ipm"))
})


test_that("validate_ipm works", {

    place <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                  "testthat", "testdata")
    file <- here(place, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- 1
    climatic <- c(sgdd = 2000, wai = -0.2)
    clim_lab <- "paradis"
    compress <- TRUE
    x <- new_ipm(IPM, BA, mesh, species, climatic, clim_lab, compress)


    expect_identical(x, validate_ipm(x))
    tmp <- x
    names(tmp) <- c("IPM", "BA", "meshpts", "climatic", "info")
    expect_error(
        validate_ipm(tmp),
        "IPM class must be composed of elements IPM, BA, mesh, climatic and info"
    )
    tmp <- x
    names(tmp$info) <- c("species", "clim_lab", "compressed")
    expect_error(
        validate_ipm(tmp),
        "IPM class must have info of elements species, climatic, compress and delay"
    )
})


test_that("old_ipm2ipm works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")
    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- c(sgdd = 2688.2488978521455464943, wai = -0.1870090464753778603, sgddb = 0.0003719893648236884,
                  waib = 1.23002598696163, wai2 = 0.0349723834636300399, sgdd2 = 7226682.1368032749742269516,
                  PC1 = -0.3447177811918831214, PC2 = 0.1229983549657221314, N = 1, SDM = 0.04901825789273488
    )
    clim_lab <- "1"
    compress <- TRUE
    delay <- 0

    expect_identical(
        old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_ipm(IPM, BA, mesh, species, climatic, clim_lab = clim_lab,
                compress = compress, delay)
    )

    expect_identical(
        old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1, delay = 2),
        delay(new_ipm(IPM, BA, mesh, species, climatic, clim_lab, compress, delay), delay = 2)
    )

})

