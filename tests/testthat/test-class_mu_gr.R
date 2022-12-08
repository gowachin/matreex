
test_that("getRangemu works", {

    species<- "Abies_alba"
    climate <- subset(treeforce::climate_species,
                      sp == species, select = -sp)
    fit <- fit_Abies_alba

    expect_equal(
        getRangemu(climate, fit),
                     c(min= -3.015334875, max = 2.277579061, sig = 0.596726204)
    )
})


test_that("make_mu_gr works", {

    f <- fit_Picea_abies

    expect_equal(
        make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                   mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                   level = c(3, 10), midbin_tresh = 2),
        structure(list(
            mu_gr = structure(c(3.71708104286679e-29, 2.82449363960135e-21,
                                1.11935798433023e-14, 2.31758162258592e-09, 2.52288935185623e-05,
                                0.0147976387634345, 0.510163925556425, 1.3383334444985, 0.764519947327727,
                                1.24618232395601e-08, 6.2383412414785e-05, 0.0179203037693807,
                                0.3404417294973, 0.539601019103744, 0.0779642652306831, 0.000853912101669262,
                                5.86099137759088e-07, 2.29177829522159e-11, 2.91148618333743e-13,
                                2.10040492675858e-08, 8.35372197595973e-05, 0.0188966476335774,
                                0.252024945238924, 0.203411215879399, 0.00993355105319125, 2.86012123207744e-05,
                                4.68802794462808e-09), dim = c(9L, 3L)),
            BA = 0:200,
            mesh = c(145.5, 256.5, 367.5, 478.5, 589.5, 700.5, 811.5, 922.5, 1033.5, 1144.5),
            mu_tab = c(-5, -4, -3, -2, -1, 0, 1, 2, 3), fit = f,
            info = c(species = "Picea_abies", correction = "constant",
                     clim_lab = "mu_gr", step = "1", surv = "TRUE"),
            int = c(gl1 = 3, gl2 = 10, gl_tresh.U = 1, mb_tresh = 2, mid_level = 5, year_delta = 1)
            ), class = "mu_gr")
    )

    expect_error(
        make_mu_gr(species = "Pseudotsuga_menziesii", fit = fit_Picea_abies,
                   mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                   level = c(3, 10), midbin_tresh = 2),
        "This species is not listed in species for which treeforce package has climate."
    )

    false_fit <- fit_Picea_abies
    false_fit$gr$params_m <- c(false_fit$gr$params_m, temp = 32)

    expect_error(
        make_mu_gr(species = "Picea_abies", fit = false_fit,
                   mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                   level = c(3, 10), midbin_tresh = 2),
        "Missing climate variables used in fit model.
Missing : temp"
    )

    expect_error(
        make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                   mesh = c(mesh = 10, L = 90, U = 1200), stepMu = 1,
                   level = c(3, 10), midbin_tresh = 2),
        "mesh must be consitued of m, L and U"
    )
})


test_that("validate_mu_gr works", {

    x <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
               mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
               level = c(3, 10), midbin_tresh = 2)
    expect_identical(x, validate_mu_gr(x))
    tmp <- x
    names(tmp) <- c("mu_gr", "BA", "mu_mesh", "mu_tab", "fit",
                    "info", "int")
    expect_error(
        validate_mu_gr(tmp),
        "mu_gr class must be composed of elements mu_gr, BA, mesh, mu_tab, fit, info and int"
    )
    tmp <- x
    names(tmp$info) <- c("sp", "correction", "clim_lab",
                         "step", "surv")
    expect_error(
        validate_mu_gr(tmp),
        "species class must have info of elements species, correction, clim_lab, step and surv"
    )
})
