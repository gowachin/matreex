test_that("get_step_IPM.mu_gr works", {

    species<- "Picea_abies"
    climate <- subset(treeforce::climate_species,
                      sp == species & N == 2, select = -sp)
    climate <- drop(as.matrix(climate))

    x <- make_mu_gr(species = "Picea_abies", fit_Picea_abies,
               mesh = c(m = 3, L = 90, U = 1200), stepMu = 1,
               level = c(3, 10), midbin_tresh = 2)

    expect_equal(
        get_step_IPM.mu_gr(x = x, BA = 20, climate = climate, sim_corr = "cut"),
        new("dtCMatrix", i = c(0L, 1L, 1L), p = c(0L, 2L, 3L, 3L),
            Dim = c(3L, 3L), Dimnames = list(NULL, NULL),
            x = c(0.0123028767143828, 0.34158918921637, 0.034689615536287),
            uplo = "L", diag = "N")
    )

})

test_that("get_step_IPM.ipm works", {

    species<- "Picea_abies"
    climate <- subset(treeforce::climate_species,
                      sp == species & N == 2, select = -sp)
    climate <- drop(as.matrix(climate))
    mesh = c(m = 10, L = 90, U = 1500)
    fit =  fit_Picea_abies

    x <- make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1:10,
                  correction = "none", level = 420, diag_tresh = -1,
                  midbin_tresh = 0)

    expect_equal(
        get_step_IPM(x = x, BA = 2, climate = climate, sim_corr = "cut"),
        new("ddiMatrix", diag = "N", Dim = c(10L, 10L), Dimnames = list(
            NULL, NULL), x = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    )

})
