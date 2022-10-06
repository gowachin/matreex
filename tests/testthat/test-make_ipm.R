test_that("gaussQuadInt works", {

    expect_equal(gaussQuadInt(10, 32, 3),
                 list(weights = c(6.11111111, 9.7777778, 6.11111111),
                      nodes = c(12.4794366, 21.00000, 29.5205634)))
})


test_that("build_weight_matrix works", {

    w <- c(6.11111111, 9.7777778, 6.11111111)

    exp <- structure(c(w, rep(0, 6), w), dim = c(6L, 2L))
    expect_equal(build_weight_matrix(w, 2), exp)
})


test_that("fun_mid_int works", {
    mesh <- 1:3
    h <- 1
    gr <- function(size){
        size * .2
    }
    sig_gr <- 0.5
    N_ini <- 1
    N_int <- 4

    exp <- matrix(c(0.628783402782724, 0.246197111495438, 0.0540056118837384,
                    0.0121485419484611, 0.50811935990912, 0.327619717560027,
                    0.100285491892948, 0.0285841489726821, 0.365150947810569,
                    0.376154577492603, 0.15955742006563, 0.0574874304830683
    ), nrow = 4)

    expect_equal(fun_mid_int(mesh, h, gr, sig_gr, N_ini, N_int, Level = 10),
                 exp)

    N_ini <- 0

    exp <- matrix(c(0.0546026335466382, 0.628783402782724, 0.246197111495438,
                    0.0540056118837384, 0.0225332638056529, 0.50811935990912,
                    0.327619717560027, 0.100285491892948, 0.00805767982538893,
                    0.36515094781057, 0.376154577492603, 0.15955742006563), nrow = 4)

    expect_equal(fun_mid_int(mesh, h, gr, sig_gr, N_ini, N_int, Level = 10),
                 exp)
})

# Common stuff for make_ipm function !
species <- "Yggdrasil"
path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
             "testthat", "testdata")
clim <- readRDS(here(path, "output", "Yggdrasil", "plots_pred.Rds"))
climate <- unlist(clim[1,])
mesh = c(m = 10, L = 90, U = 1500)
fit <- list(
    sv = list(
        params_m = c(intercept = 10.247, size = 0.00819,
                     logsize = -3.72, BATOTcomp = 0.00811,
                     sgdd = -0.00261, wai = 1.507,
                     sgdd2 = -3.10-07, wai2 = 0.593,
                     `size:wai` = 0.00103, `logsize:wai` = -0.500,
                     `size:sgdd` = -1.52e-06, `logsize:sgdd` = 0.000834),
        family = list(linkinv = function (eta) {
            pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
        })
    ),
    gr = list(params_m = c(intercept = -5.039, size = -0.00154,
                           logsize = 0.985, BATOTcomp = -0.0152,
                           sgdd = 0.000572, wai = -0.218, sgdd2 = -1.29e-07,
                           wai2 = -0.181, `size:wai` = 0.000143,
                           `logsize:wai` = 0.0927, `size:sgdd` = -2.92e-07,
                           `logsize:sgdd` = 9.08e-05),
              sigma = 0.622),
    rec = list(params_m = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2)))



test_that("make_IPM works : fonctions communes", {

    ex <- Matrix(matrix(0, ncol = 10, nrow = 10), sparse = TRUE)

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
                IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
                climatic = climate, clim_lab = "test",
                rec = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2),
                species = species, compress = FALSE
            ))
    )

    # survival is FALSE
    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0, IsSurv = FALSE),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
            rec = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2),
            species = species, compress = FALSE
        ))
    )
})


test_that("make_IPM works : mid_bin", {

    ex <- new("dtCMatrix",
              i = c(0L, 1L,
                    2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                    9L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
                    4L, 5L, 6L, 7L, 8L, 9L, 5L, 6L, 7L, 8L, 9L, 6L, 7L, 8L, 9L, 7L,
                    8L, 9L, 8L, 9L, 9L),
              p = c(0L, 10L, 19L, 27L, 34L, 40L, 45L, 49L, 52L, 54L, 55L),
              Dim = c(10L, 10L), Dimnames = list(NULL, NULL),
              x = c(1.01333114220533, 1.61045146216677e-08, 1.61320196334519e-13,
                    2.64333452316302e-16, 2.65234143094034e-18, 6.980403200971e-20,
                    3.37267989949e-21, 2.48371086225448e-22, 2.49690542628745e-23,
                    3.19176831738741e-24, 1.00210128362939, 1.61160301667125e-05,
                    1.67627200513434e-09, 8.21181222733903e-12, 1.69306363131613e-13,
                    7.62453850172581e-15, 5.65627501055709e-16, 5.95319771207903e-17,
                    8.12869240486978e-18, 1.00022515448486, 7.48758019845865e-05,
                    1.44430673305361e-08, 9.45222829274529e-11, 2.35811721904441e-12,
                    1.22444702546171e-13, 1.01781440097791e-14, 1.17781037605766e-15,
                    1.00010722924696, 9.03813351108587e-05, 1.88698207220419e-08,
                    1.28166158976529e-10, 3.27663167124276e-12, 1.7327618116984e-13,
                    1.46154286405327e-14, 1.00042131597956, 5.7567712496093e-05,
                    9.95581924782902e-09, 6.19016825678057e-11, 1.49308009368092e-12,
                    7.55985272647216e-14, 1.00141196795186, 2.41785907927703e-05,
                    2.94691265419789e-09, 1.55500331121483e-11, 3.36668075217783e-13,
                    1.00393820357338, 7.31035708195537e-06, 5.63018923561546e-10,
                    2.39594441241346e-12, 1.00914606999184, 1.65299717183477e-06,
                    7.43244156964355e-11, 1.01640114121826, 2.84286807834341e-07,
                    1.01804165011446), uplo = "L", diag = "N")

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 12),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
            rec = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2),
            species = species, compress = FALSE
        ))
    )

    res<-evaluate_promise({
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 12, verbose = TRUE)
    })

    expect_equal(res$messages[1], "Launching integration loop\n")
    expect_equal(res$messages[2], "GL integration won't occur because of negative treshold\n")
    expect_equal(res$messages[3], "midbin integration occur on 10 cells\n")
    expect_equal(res$messages[4], "Loop done.\n")

})


test_that("make_IPM works : gauss-legendre", {

    ex <- new("dtCMatrix",
              i = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
                    1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                    9L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 4L, 5L, 6L, 7L, 8L, 9L, 5L, 6L,
                    7L, 8L, 9L, 6L, 7L, 8L, 9L, 7L, 8L, 9L, 8L, 9L, 9L),
              p = c(0L, 10L, 19L, 27L, 34L, 40L, 45L, 49L, 52L, 54L, 55L),
              Dim = c(10L, 10L),
              Dimnames = list(NULL, NULL),
              x = c(0.987873947547578, 0.00810661999982615,
                    3.47635662166081e-09, 5.54583300226408e-12, 6.72529101190695e-14,
                    2.19520784744741e-15, 1.30519172219056e-16, 1.16565513329924e-17,
                    1.40006631794013e-18, 2.10913283764152e-19, 0.97839383949078,
                    0.0206743766398513, 4.34949552351638e-08, 1.15214394053169e-10,
                    1.91526908493958e-12, 7.8458489611432e-14, 5.56426609267383e-15,
                    5.73582104510589e-16, 7.77177256219818e-17, 0.973930230276361,
                    0.0256545390341093, 8.43080602926875e-08, 2.6907314125609e-10,
                    5.09003930211168e-12, 2.29850332856189e-13, 1.76093619277022e-14,
                    1.93399669311651e-15, 0.97600394888106, 0.0235788815457929, 6.84027400173429e-08,
                    2.21107963731384e-10, 4.28622847986053e-12, 1.97996733591586e-13,
                    1.54649395916916e-14, 0.981363468520974, 0.0179422878510356,
                    3.24602456928403e-08, 9.52781209727371e-11, 1.76988672702588e-12,
                    7.96614663705882e-14, 0.986926197921166, 0.0116988880646031,
                    1.04185542456469e-08, 2.57805731886707e-11, 4.38539732042285e-13,
                    0.990734397864969, 0.00659128451113505, 2.43087957418674e-09,
                    4.80576731717906e-12, 0.992316908839443, 0.00320329048782211,
                    4.28583052972224e-10, 0.993374287906479, 0.00133498605576727,
                    0.998860122517894), uplo = "L", diag = "N")

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = 1500,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
            rec = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2),
            species = species, compress = FALSE
        ))
    )

    res<-evaluate_promise({
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = 1500,
                 midbin_tresh = 0, verbose = TRUE)
    })

    expect_equal(res$messages[1], "Launching integration loop\n")
    expect_equal(res$messages[2], "GL integration occur on 10 cells\n")
    expect_equal(res$messages[3], "midbin integration won't occur because of treshold at 0\n")
    expect_equal(res$messages[4], "Loop done.\n")

})



test_that("make_IPM works : corrections", {

    ex <- Matrix(matrix(0, ncol = 10, nrow = 10), sparse = TRUE)

    # expect_equal(
    #     make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
    #              correction = "constant", level = 420, diag_tresh = -1,
    #              midbin_tresh = 0),
    #     validate_ipm(new_ipm(
    #         IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
    #         climatic = climate, clim_lab = "test",
    #         species = species, compress = FALSE
    #     ))
    # )

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "sizeExtremes", level = 420, diag_tresh = -1,
                 midbin_tresh = 0, IsSurv = FALSE),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
            rec = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2),
            species = species, compress = FALSE
        ))
    )

    # expect_equal(
    #     make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
    #              correction = "ceiling", level = 420, diag_tresh = -1,
    #              midbin_tresh = 0, IsSurv = FALSE),
    #     validate_ipm(new_ipm(
    #         IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
    #         climatic = climate, clim_lab = "test",
    #         species = species, compress = FALSE
    #     ))
    # )
})


test_that("make_IPM stops", {

    expect_error(
        make_IPM(species, c(temp = 666), clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0),
        paste0("Missing climate variables used in fit model.\n",
               "Missing : sgdd wai sgdd2 wai2")
    )

    expect_error(
        make_IPM(species, climate, clim_lab = "test", fit,
                 mesh = c(length = 10, L = 90, U = 1500), BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0),
        paste0("mesh must be consitued of m, L and U")
    )

})

