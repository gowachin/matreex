sv_params <- gr_params <- c(intercept = 2, size = .5)
rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
sv_family <- list(linkinv = function(eta){
    pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
})
gr_sigma = .6
max_dbh <- 42

fit <- new_fit_sgr(sv_params, sv_family,
                             gr_params, gr_sigma, rec_params,
                             "Yggdrasil", 42)


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

    expect_identical(new_ipm(IPM = IPM, BA = BA, mesh = mesh, species = species,
                             correction = "constant", climatic = climatic,
                             clim_lab = clim_lab, fit = fit, compress = compress),
                     structure(list(
                         IPM = IPM, BA = BA, mesh = mesh, climatic = climatic,
                         fit = fit,
                         info = c(species = "Yggdrasil", correction = "constant",
                                  clim_lab = "paradis", surv = "TRUE",
                                  compress = "TRUE", delay = "0"),
                         int = c(gl1 = 0, gl2 = 0, gl_tresh = 0, gl_min = 0,
                                  mb_tresh = 0, mid_level = 0, mb_max = 0,
                                  year_delta = 0, max_error = 0)),
                         class = "ipm"))
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
    x <- new_ipm(IPM = IPM, BA = BA, mesh = mesh, species = species,
                 correction = "constant", climatic = climatic,
                 clim_lab = clim_lab, fit = fit, compress = compress)


    expect_identical(x, validate_ipm(x))
    tmp <- x
    names(tmp) <- c("IPM", "BA", "meshpts", "climatic", "info", "int")
    expect_error(
        validate_ipm(tmp),
        "IPM class must be composed of elements IPM, BA, mesh, climatic, fit, info and int"
    )

    tmp <- x
    names(tmp$info) <- c("species", "clim_lab", "compressed")
    expect_error(
        validate_ipm(tmp),
        "IPM class must have info of elements species, correction, clim_lab, surv, compress and delay"
    )


})


test_that("old_ipm2ipm works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")
    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]
    LIPM <- IPM$LIPM
    BA <- 1:length(LIPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- c(sgdd = 2688.2488978521455464943, wai = -0.1870090464753778603, sgddb = 0.0003719893648236884,
                  waib = 1.23002598696163, wai2 = 0.0349723834636300399, sgdd2 = 7226682.1368032749742269516,
                  PC1 = -0.3447177811918831214, PC2 = 0.1229983549657221314, N = 1, SDM = 0.04901825789273488
    )
    clim_lab <- "1"
    fit <- old_fit2fit(species, path = path, replicat = 1, mean = FALSE)
    compress <- TRUE
    delay <- 0


    expect_equal(
        old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_ipm(IPM = LIPM, BA = BA, mesh = mesh, species = species,
                correction = "constant", climatic = climatic,
                clim_lab = clim_lab, fit = fit, compress = compress)
    )

    expect_equal(
        old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1, delay = 2),
        delay(new_ipm(IPM = LIPM, BA = BA, mesh = mesh, species = species,
                      correction = "constant", climatic = climatic,
                      clim_lab = clim_lab, fit = fit, compress = compress), delay = 2)
    )

})

