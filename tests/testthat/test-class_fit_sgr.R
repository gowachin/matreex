test_that("new_fit_sgr works", {

    sv_params <- gr_params <- c(intercept = 2, size = .5)
    rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    sv_family <- list(linkinv = function(eta){
        pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
    })
    gr_sigma = .6
    species <- "shadok"
    max_dbh <- 42

    expect_identical(new_fit_sgr(sv_params, sv_family,
                                 gr_params, gr_sigma, rec_params,
                                 "shadok", 42),
                     structure(list(
                         sv = list(params_m = sv_params, family = sv_family),
                         gr = list(params_m = gr_params, sigma = gr_sigma),
                         rec = list(params_m = rec_params),
                         info = c(species = "shadok", max_dbh = "42")
                     ),
                     class = "fit_sgr"))

    expect_identical(fit_sgr(sv_params, sv_family,
                                 gr_params, gr_sigma,
                             rec_params,
                                 "shadok", 42),
                     new_fit_sgr(sv_params, sv_family,
                             gr_params, gr_sigma,
                             rec_params,
                             "shadok", 42))
})


test_that("validate_fit_sgr works", {

    sv_params <- gr_params <- c(intercept = 2, size = .5)
    rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    sv_family <- list(linkinv = function(eta){
        pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
    })
    gr_sigma = .6

    x <- new_fit_sgr(sv_params, sv_family, gr_params, gr_sigma, rec_params, "shadok", 42)


    expect_identical(x, validate_fit_sgr(x))
    tmp <- x
    names(tmp) <- c("survival", "gr", "info")
    expect_error(
        validate_fit_sgr(tmp),
        "fit_sgr class must be composed of elements sv, gr, rec and info"
    )
    tmp <- x
    tmp$sv$params_m <- c(Intercept = 2, size = .5)
    expect_error(
        validate_fit_sgr(tmp),
        "Survival model should at least depend on an intercept and the size variable."
    )
    tmp <- x
    tmp$gr$params_m <- c(Intercept = 2, size = .5)
    expect_error(
        validate_fit_sgr(tmp),
        "Growth model should at least depend on an intercept and the size variable."
    )

    tmp <- x
    tmp$rec$params_m <- c(intercept = 2, BATOTSP = 0, BATOTnonSP = 1)
    expect_error(
        validate_fit_sgr(tmp),
        "Recruitment model should at least depend on an intercept, BATOTSP and BATOTNonSP variable."
    )

    tmp <- x
    names(tmp$info) <- c("sp", "max_dbh")
    expect_error(
        validate_fit_sgr(tmp),
        "fit_sgr class must have info of elements species and max_dbh"
    )
})


test_that("mean_oldfit works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    species = "Yggdrasil"
    max_dbh <- 1600
    file <-  here(path, "output", species, "fit_sgr_all.Rds")
    fit <- readRDS(file)

    exp <- structure(list(
        sv = list(params_m = c(BATOTcomp = 0.00788664089885266,
            intercept = 11.9013764892525, logsize = -3.93693226264533, `logsize:sgdd` = 0.000706566602363076,
            `logsize:sgddb` = -300.240902250853, `logsize:wai` = -0.254715109352006,
            `logsize:waib` = 0.58900275676364, sgdd = -0.00270369231603175,
            sgdd2 = -1.49517976896259e-07, sgddb = 1280.8317883853, size = 0.00883370086149075,
            `size:sgdd` = -1.37614662736779e-06, `size:sgddb` = 0.654994130494874,
            `size:wai` = 0.000525819849898456, `size:waib` = -0.00127077266533475,
            wai = 0.634372719503291, wai2 = 0.393938802895601, waib = -2.61369245577149
        ), family = structure(list(
            family = "binomial", link = "cloglog",
            linkinv = function (eta) {
                pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
            }), class = "family")),
        gr = list(params_m = c(BATOTcomp = -0.0153566533095071,
            intercept = -4.96490345545489, logsize = 0.969879429895294, `logsize:sgdd` = 0.000101550631573779,
            `logsize:wai` = 0.117461947720243, sgdd = 0.000512040990260767,
            sgdd2 = -1.2621957152595e-07, size = -0.00146753201555766, `size:sgdd` = -3.43810044991891e-07,
            `size:wai` = 6.95635485702585e-05, wai = -0.3283117658971, wai2 = -0.184076301010191
        ), sigma = 0.6229243935538504),
        rec = list(params_m = c(BATOTNonSP = -0.0248842776930938, BATOTSP = -0.0190822994852905,
                                intercept = -0.908737368267839, logBATOTSP = -0.252604353458172,
                                sgddb = 321.429377000321, wai = -0.0614645745262896, wai2 = 0.290906922523664
        )),
        info = c(species = "Yggdrasil", max_dbh = "1600")), class = "fit_sgr")


    res <- mean_oldfit(fit, species, max_dbh)
    res$sv$family <- structure(list(
        family = "binomial", link = "cloglog",
        linkinv = function (eta) {
            pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
        }), class = "family")

    expect_equal(res, exp)
})



test_that("old_fit2fit works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    species = "Yggdrasil"
    max_dbh <- 1649
    file <-  here(path, "output", species, "fit_sgr_all.Rds")
    fit <- readRDS(file)
    replicat <- 42

    res_fit <- fit[[replicat]]

    exp <- fit_sgr(res_fit$sv$params_m, res_fit$sv$family,
                   res_fit$gr$params_m, res_fit$gr$sigma,
                   res_fit$rec$params_m,
                   species = species, max_dbh = 1544)

    expect_equal(old_fit2fit(species, path, replicat = 42, mean = FALSE), exp)

    expect_equal(old_fit2fit(species, path, mean = TRUE),
                 mean_oldfit(fit, species, max_dbh))

})
