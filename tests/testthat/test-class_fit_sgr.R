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
                                 "shadok", 42, 0),
                     structure(list(
                         sv = list(params_m = sv_params, family = sv_family),
                         gr = list(params_m = gr_params, sigma = gr_sigma),
                         rec = list(params_m = rec_params),
                         info = c(species = "shadok", max_dbh = "42", delay = "0")
                     ),
                     class = "fit_sgr"))

    expect_identical(fit_sgr(sv_params, sv_family,
                                 gr_params, gr_sigma,
                             rec_params,
                                 "shadok", 42, 0),
                     new_fit_sgr(sv_params, sv_family,
                             gr_params, gr_sigma,
                             rec_params,
                             "shadok", 42, 0))
})


test_that("validate_fit_sgr works", {

    sv_params <- gr_params <- c(intercept = 2, size = .5)
    rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    sv_family <- list(linkinv = function(eta){
        pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
    })
    gr_sigma = .6

    x <- new_fit_sgr(sv_params, sv_family, gr_params, gr_sigma, rec_params, "shadok", 42, 0)


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
    names(tmp$info) <- c("sp", "max_dbh", "delay")
    expect_error(
        validate_fit_sgr(tmp),
        "fit_sgr class must have info of elements species, max_dbh and delay"
    )
})


test_that("mean_oldfit works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    species = "Yggdrasil"
    max_dbh <- 1600
    file <-  here(path, "output", species, "fit_sgr_all.Rds")
    fit <- readRDS(file)

    exp <- structure(
        list(
            sv = list(
                params_m = c(
                    BATOTcomp = 0.0080549554722945, intercept = 10.9945773188906,
                    logsize = -3.75929970819068, `logsize:sgdd` = 0.000417191231113778,
                    `logsize:sgddb` = -126.484281755702, `logsize:wai` = -0.250138019156562,
                    `logsize:waib` = 0.896224691499341, sgdd = -0.0013050149281753,
                    sgdd2 = -1.55165230858647e-07, sgddb = 309.189963118227, size = 0.0085064699761369,
                    `size:sgdd` = -7.62488191975794e-07, `size:sgddb` = 0.385320322502985,
                    `size:wai` = 0.000516298446070642, `size:waib` = -0.00198689723046645,
                    wai = 0.753649327992577, wai2 = 0.296526935867094, waib = -3.96730108212054
                ), family = structure(
                    list(family = "binomial", link = "cloglog", linkinv = function (eta)
                        pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)), class = "family")),
            gr = list(params_m = c(BATOTcomp = -0.0152773871016923, intercept = -4.94352065003243,
                                   logsize = 0.963767044418592, `logsize:sgdd` = 0.000105410889632953,
                                   `logsize:wai` = 0.101311873777565, sgdd = 0.000514978757323322,
                                   sgdd2 = -1.34024298547651e-07, size = -0.00154789437793654,
                                   `size:sgdd` = -3.10421883815199e-07, `size:wai` = 0.000140274923123379,
                                   wai = -0.263098590772547, wai2 = -0.181252473633295), sigma = 0.621825127187286),
            rec = list(params_m = c(BATOTNonSP = -0.024793354587699,
                                    BATOTSP = -0.0187862744077507, intercept = -0.918591942470068,
                                    logBATOTSP = -0.257347121425412, sgddb = 295.822304118005,
                                    wai = -0.0643689062773369, wai2 = 0.303847830350663)),
            info = c(species = "Yggdrasil", max_dbh = "1600", delay = "0")), class = "fit_sgr")


    res <- mean_oldfit(fit, species, max_dbh, delay = 0)
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
    replicat <- 2

    res_fit <- fit[[replicat]]

    exp <- fit_sgr(res_fit$sv$params_m, res_fit$sv$family,
                   res_fit$gr$params_m, res_fit$gr$sigma,
                   res_fit$rec$params_m,
                   species = species, max_dbh = 1649, delay = 0)

    expect_equal(old_fit2fit(species, path, replicat = 2, mean = FALSE), exp)

    expect_equal(old_fit2fit(species, path, mean = TRUE),
                 mean_oldfit(fit, species, max_dbh, 0))

})
