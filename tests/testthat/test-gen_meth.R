
# sp_name ####
test_that("sp_name works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                 "testthat", "testdata")

    x <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(sp_name(x), "Yggdrasil")
    x <- old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(sp_name(x), "Yggdrasil")
    x <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                    mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                    level = c(3, 10), midbin_tresh = 2)
    expect_identical(sp_name(x), "Picea_abies")
})

test_that("climatic works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                 "testthat", "testdata")

    x <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(climatic(x), "1")
    x <- old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(climatic(x), "1")
    x <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                    mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                    level = c(3, 10), midbin_tresh = 2)
    expect_identical(climatic(x), "mu_gr")
})

sv_params <- gr_params <- c(intercept = 2, size = .5)
rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
sv_family <- list(linkinv = function(eta){
    pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
})
gr_sigma = .6
max_dbh <- 42

fit <- new_fit_sgr(sv_params, sv_family,
                   gr_params, gr_sigma, rec_params,
                   "Yggdrasil", 42, 0)

# delay ####
test_that("delay dtCMatrix works", {
    x <- new("dtCMatrix",
             i = c(0L, 1L, 2L, 1L, 2L, 2L),
             p = c(0L, 3L,  5L, 6L),
             Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")

    expect_identical(
        delay(x, 2),
        new("dgCMatrix",
            i = c(1L, 2L, 2L, 3L, 4L, 3L, 4L, 4L),
            p = c(0L, 1L, 2L, 5L,  7L, 8L),
            Dim = c(5L, 5L),
            x = c(1, 1, 1, 2, 3, 1, 2, 1))
    )

    expect_identical(
        delay(delay(x, 2), -1),
        new("dtCMatrix",
            i = c(1L, 1L, 2L, 3L, 2L, 3L, 3L),
            p = c(0L, 1L, 4L,  6L, 7L),
            Dim = c(4L, 4L), x = c(1, 1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    )

    expect_identical( delay(x, 0), x )
})


test_that("delay numeric works", {
    x <- as.numeric(1:10)

    expect_identical( delay(x, 2), c(0,0, x) )
    expect_identical( delay(delay(x, 2), -1), c(0, x) )
    expect_identical( delay(x, 0), x )
})


test_that("delay ipm works", {

    rec <-  fit
    mat<- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
              p = c(0L, 3L,  5L, 6L),
              Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", correction = "constant",
                 climatic = 1, clim_lab = "1", fit = rec,
                 compress = FALSE)

    exp <- new_ipm(IPM = list( delay(mat, 2) ), BA = 1, mesh = c(0,0,1:3),
                   species = "darwin", correction = "constant",
                   climatic = 1, clim_lab = "1", fit = rec,
                   compress = FALSE, delay = 2)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( delay(x, 2), exp )
    expect_identical( delay(x, 0), x )

    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", correction = "constant",
                 climatic = 1, clim_lab = "1", fit = rec,
                 compress = TRUE)

    exp <- new_ipm(IPM = list( delay(mat * 1e-7, 2) ), BA = 1, mesh = c(0,0,1:3),
                   species = "darwin", correction = "constant",
                   climatic = 1, clim_lab = "1", fit = rec,
                   compress = FALSE, delay = 2)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( delay(x, 2), exp )

    expect_error( delay(delay(x, 2), -3),
                  "Negative delay is not possible for ipm objects. Minimal value here is 2")
})


test_that("delay species works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = 1, clim_lab = "1",
    fit = fit,
    compress = FALSE)
    validate_ipm(IPM)

    x <- new_species(IPM = IPM, init_pop = def_init,
                     harvest_fun = def_harv,
                     disturb_fun = def_disturb)

    exp <- new_species(delay(IPM, 2), init_pop = def_init,
                       harvest_fun = def_harv,
                       disturb_fun = def_disturb)
    # validate_species(x)
    # validate_species(exp)

    expect_identical( delay(x, 2), exp )
    expect_identical( delay(x, 0), x )

})


test_that("delay forest works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = c(sgdd = 1), clim_lab = "1",
    fit = fit,
    compress = FALSE)
    validate_ipm(IPM)

    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv,
                      disturb_fun = def_disturb)

    x <- new_forest(list(darwin = sp), favoured_sp = c(darwin = FALSE))
    exp <- new_forest(list(darwin = delay(sp, 2)), favoured_sp = c(darwin = FALSE))
    # hack ne pas nommer les elements ici fout la merde.

    # validate_forest(x)
    # validate_forest(exp)

    expect_identical( delay.forest(x, 2), exp )
    expect_identical( delay(x, 0), x )

})


# correction ####
test_that("correction ipm works", {

    rec <-  fit
    mat<- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
              p = c(0L, 3L,  5L, 6L),
              Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    new_mat<- new("dtCMatrix", i = c(0L, 1L, 1L),
                  p = c(0L, 2L,  3L, 3L),
                  Dim = c(3L, 3L), x = c(1, 2, 1), uplo = "L", diag = "N")
    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin",  correction = "constant",
                 climatic = 1, clim_lab = "1", fit = rec, compress = FALSE)

    exp <- new_ipm(IPM = list( new_mat  ), BA = 1, mesh = 1:3,
                   species = "darwin",  correction = "constant",
                   climatic = 1, clim_lab = "1", fit = rec, compress = FALSE)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( correction(x, "cut"), exp )
    expect_identical( correction(x), x )

    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin",  correction = "constant",
                 climatic = 1, clim_lab = "1", fit = rec, compress = TRUE)

    exp <- new_ipm(IPM = list( new_mat * 1e-7 ), BA = 1, mesh = 1:3,
                   species = "darwin",  correction = "constant",
                   climatic = 1, clim_lab = "1", fit = rec, compress = FALSE)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( correction(x, "cut"), exp )
})

test_that("correction mu_sgr works", {

    x <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                    mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                    level = c(3, 10), midbin_tresh = 2)
    expect_identical(x, correction(x))

})

test_that("correction species works", {

    rec <-  c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = 1, clim_lab = "1", fit = fit, compress = FALSE)
    validate_ipm(IPM)


    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv,
                      disturb_fun = def_disturb)

    exp <- new_species(correction(IPM, "cut"), init_pop = def_init,
                       harvest_fun = def_harv,
                       disturb_fun = def_disturb)
    # validate_species(sp)
    # validate_species(exp)

    expect_identical( correction(sp, "cut"), exp )
    expect_identical( correction(sp), sp )

})


test_that("correction forest works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = 1, clim_lab = "1",
    fit = fit, compress = FALSE)
    validate_ipm(IPM)

    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv,
                      disturb_fun = def_disturb)

    x <- new_forest(list(darwin = sp))
    exp <- new_forest(list(darwin = correction(sp, "cut")))
    # hack ne pas nommer les elements ici fout la merde.

    expect_identical( correction(x, "cut"), exp )
    expect_identical( correction(x), x )

})


# sp_rec ####
test_that("sp_rec mu_gr works", {

    x <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                    mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                    level = c(3, 10), midbin_tresh = 2)
    climate <- c(sgdd = 1444.66662815716, wai = 0.451938692369727, sgddb = 0.000692201218266957,
                 waib = 0.688734314510131, wai2 = 0.204248581660859, sgdd2 = 2087061.66651097,
                 PC1 = 1.67149830316836, PC2 = 0.0260206360117464, N = 2, SDM = 0.676055555555556
    )

    expect_identical(
        exp_recFun(params = x$fit$rec$params_m, list_covs = climate),
        sp_rec(x, climate)
    )

})


sv_params <- gr_params <- c(intercept = 2, size = .5)
rec_params <- c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
sv_family <- list(linkinv = function(eta){
    pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
})
gr_sigma = .6
max_dbh <- 42

fit <- new_fit_sgr(sv_params, sv_family,
                   gr_params, gr_sigma, rec_params,
                   "Yggdrasil", 42, 0)

test_that("sp_rec species works", {

    climate <- c(sgdd = 1444.66662815716, wai = 0.451938692369727, sgddb = 0.000692201218266957,
                 waib = 0.688734314510131, wai2 = 0.204248581660859, sgdd2 = 2087061.66651097,
                 PC1 = 1.67149830316836, PC2 = 0.0260206360117464, N = 2, SDM = 0.676055555555556
    )

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = climate, clim_lab = "1", fit = fit, compress = FALSE)
    validate_ipm(IPM)

    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv,
                      disturb_fun = def_disturb)

    expect_identical(
        sp$recruit_fun,
        sp_rec(sp, climate)
    )

    x <- make_mu_gr(species = "Picea_abies", fit_Picea_abies,
                    mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                    level = c(3, 10), midbin_tresh = 2)
    sp <- new_species(IPM = x, init_pop = def_init,
                      harvest_fun = def_harv,
                      disturb_fun = def_disturb)

    expect_identical(
        sp_rec(x, climate),
        sp_rec(sp, climate)
    )

})


# Max dbh ####
test_that("sp_rec species works", {

    expect_identical(
        get_maxdbh(fit_Abies_alba),
        1420
    )
})

# Summary
test_that("summary IPM works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = 1, clim_lab = "1",
    fit = fit, compress = FALSE)
    validate_ipm(IPM)



    res<-evaluate_promise({
        summary(IPM)
    })

    expect_equal(
        res$messages[1],
        paste0("IPM object for species darwin at climate '1' \n\n",
               "Integation was done on a mesh from 1.00 to 3.00 with 3 cells for BA between 1 and 1. \n",
               "Gauss-Legendre was used on 0 cells with 0 x 0 levels (size at t * t+1) \n",
               "Midbin was used on 0 cells with 0 levels. \n",
               "The correction was constant, the IPM does contain survival and the recruitment delay is 0. \n\n"
        )
    )
})

test_that("summary IPM works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin",  correction = "constant",
    climatic = 1, clim_lab = "1",
    fit = fit, compress = FALSE)

    sp <- species(IPM)
    validate_species(sp)



    res<-evaluate_promise({
        summary(sp)
    })

    expect_equal(
        res$messages[1],
        paste0("IPM object for species darwin at climate '1' \n\n",
               "For Uneven harvest, dth = 175.0, dha = 575.0 and hmax = 1.00.\n",
               "rdi_coef are not defined for this species. Even harvest not possible\n",
               "disturb_coef are not defined for this species. Disturbance models not possible\n\n"
        )
    )

    expect_equal(
        res$messages[2],
        paste0("Integation was done on a mesh from 1.00 to 3.00 with 3 cells for BA between 1 and 1. \n",
               "Gauss-Legendre was used on 0 cells with 0 x 0 levels (size at t * t+1) \n",
               "Midbin was used on 0 cells with 0 levels. \n",
               "The correction was constant, the IPM does contain survival and the recruitment delay is 0. \n\n"
        )
    )

})

