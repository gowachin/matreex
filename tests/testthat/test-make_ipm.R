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
    rec = list(params_m = c(intercept = 1, BATOTSP = 1, BATOTNonSP = 2)),
    info = c(species = "Yggdrasil", max_dbh = 42, delay = "0"))
class(fit) <- "fit_sgr"
int_log <- c(year_delta = 1, MaxError = 1,
             GL_Nint = 0, GL_level = 420, GL_min = 0,
             MB_Nint = 0, MB_level = 5, MB_max = 0)
int <- c(gl1 = 3, gl2 = 140, gl_tresh = 0, gl_min = 0,
         mb_tresh = 0, mid_level = 5, mb_max = 0,
         year_delta = 1, max_error = 1)

test_that("make_IPM works : fonctions communes", {

    ex <- Matrix(matrix(0, ncol = 10, nrow = 10), sparse = TRUE)

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
                IPM = list(`1` = ex), BA = 1, mesh = seq(160.5, 1429.5, length.out = 10),
                climatic = climate, clim_lab = "test", fit = fit,
                species = species, correction = "none",
                compress = FALSE, int = int
            ))
    )

    # survival is FALSE
    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0, IsSurv = FALSE),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(160.5, 1429.5, length.out = 10),
            climatic = climate, clim_lab = "test", fit = fit,
            species = species, correction = "none",
            compress = FALSE, int = int, survival = FALSE
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
              p = c(0L, 10L, 19L, 27L, 34L, 40L, 45L,
                    49L, 52L, 54L, 55L),
              Dim = c(10L, 10L), Dimnames = list(NULL, NULL),
              x = c(0.0177923861010391, 2.15897414025614e-09, 7.01713402329326e-14,
                    1.53910645070346e-16, 1.76631330180612e-18, 5.02602662431653e-20,
                    2.55656515818775e-21, 1.95296951797064e-22, 2.01815614912049e-23,
                    2.6358842010398e-24, 0.385658620814465, 3.59489520443886e-06,
                    8.57403445923304e-10, 5.24752349068835e-12, 1.20278467647037e-13,
                    5.76784829653512e-15, 4.46231924651537e-16, 4.84064118259471e-17,
                    6.7626138992482e-18, 0.639261189470722, 1.90190906415751e-05,
                    7.70412906772501e-09, 6.18849303706583e-11, 1.70386650365037e-12,
                    9.38382294686819e-14, 8.1146437521555e-15, 9.66207448149019e-16,
                    0.67573893109161, 2.33395471607757e-05, 1.01194487834026e-08,
                    8.41726166950068e-11, 2.37268045758997e-12, 1.33015148795711e-13,
                    1.16680504470854e-14, 0.590105265258337, 1.42927941101703e-05,
                    5.27144561392864e-09, 4.03546846602714e-11, 1.07560806686453e-12,
                    5.78038011370428e-14, 0.444479914720633, 5.57730476889259e-06,
                    1.52369885296707e-09, 9.99897245411065e-12, 2.40218348645655e-13,
                    0.287922687051759, 1.52969813125378e-06, 2.8212775910506e-10,
                    1.51303497763449e-12, 0.158131275591236, 3.0812990531071e-07,
                    3.58919090657476e-11, 0.0724572291183668, 4.65251304814939e-08,
                    0.0272823103923344), uplo = "L", diag = "N")
    int["max_error"] <- 0.98220761
    int["mb_max"] <- 0.67573893
    int["mb_tresh"] <- 10

    res<-evaluate_promise({
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 12)
    })

    expect_equal(
        res$result,
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(160.5, 1429.5, length.out = 10),
            climatic = climate, clim_lab = "test", fit = fit,
            species = species, correction = "none",
            compress = FALSE, int = int
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
    expect_equal(res$warnings[1], "At least one mid_bin integration has values above 1e-2. This is linked with insufficient Gauss-Legendre integration treshold. value : -1 mm.")

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
    int["max_error"] <- 4.019428970623240005011211906e-03
    int["gl_min"] <- 8.7694376e-23
    int["gl_tresh"] <- 10
    int["mb_max"] <- 0
    int["mb_tresh"] <- 0

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = 1500,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(160.5, 1429.5, length.out = 10),
            climatic = climate, clim_lab = "test", fit = fit,
            species = species, correction = "none",
            compress = FALSE, int = int
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
    int["max_error"] <- 1
    int["gl_min"] <- 0
    int["gl_tresh"] <- 0

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
            IPM = list(`1` = ex), BA = 1, mesh = seq(160.5, 1429.5, length.out = 10),
            climatic = climate, clim_lab = "test", fit = fit,
            species = species, correction = "sizeExtremes",
            compress = FALSE, int = int, surv = FALSE
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

