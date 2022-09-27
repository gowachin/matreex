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

    exp <- matrix(c(0.313617744964713, 0.26270420944295, 0.0579194449143623, 0.0129767617813288,
                    0.29226166432228, 0.343017667025785, 0.106182904046211, 0.0302430842263083,
                    0.234061337356817, 0.386723531757648, 0.166827434994055, 0.0602530149602582),
                  nrow = 4)

    expect_equal(fun_mid_int(mesh, h, gr, sig_gr, N_ini, N_int, Level = 100),
                 exp)
})

# Common stuff for make_ipm function !
species <- "Yggdrasil"
climate <- c(sgdd = 2688.248, wai = -0.187,
             sgddb = 0.0003719895, waib = -5.347594,
             wai2 = 0.034969, sgdd2 = 7226677)
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
              sigma = 0.622))

test_that("make_IPM works : fonctions communes", {

    ex <- Matrix(matrix(0, ncol = 10, nrow = 10), sparse = TRUE)

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
                IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
                climatic = climate, clim_lab = "test",
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
            species = species, compress = FALSE
        ))
    )
})


test_that("make_IPM works : mid_bin", {

    ex <- Matrix(matrix(c(
        1.60412995318835e-11, 1.61324689407932e-13, 2.64341611026272e-16,
        2.65242854955873e-18, 6.98064278509081e-20, 3.37279963023509e-21,
        2.48380146935195e-22, 2.49699861241952e-23, 3.1918897832306e-24,
        4.94168889504289e-25, 0, 6.9049138561878e-08, 1.67632668729453e-09,
        8.21211513803251e-12, 1.693130836345e-13, 7.62485713181455e-15,
        5.65652085690709e-16, 5.95346477055914e-17, 8.1290667843527e-18,
        1.35861358263381e-18, 0, 0, 4.70272297886378e-07, 1.4443592433973e-08,
        9.45261951071426e-11, 2.35822263343473e-12, 1.22450479298458e-13,
        1.01786443399612e-14, 1.17787021721043e-15, 1.74462929581719e-16,
        0, 0, 0, 5.96108753996694e-07, 1.88705858588181e-08, 1.28172080625983e-10,
        3.27679526396673e-12, 1.73285314494444e-13, 1.46162315234124e-14,
        1.71211267503423e-15, 0, 0, 0, 0, 3.37964922093222e-07, 9.95627051291696e-09,
        6.19048727530764e-11, 1.49316314605996e-12, 7.56029633062582e-14,
        6.1591516955936e-15, 0, 0, 0, 0, 0, 1.14296558855166e-07, 2.9470620556131e-09,
        1.5550926035931e-11, 3.36688897704768e-13, 1.57262919864518e-14,
        0, 0, 0, 0, 0, 0, 2.59909606071119e-08, 5.63050816771989e-10,
        2.39609742012653e-12, 4.50311893474017e-14, 0, 0, 0, 0, 0, 0,
        0, 4.20490074524273e-09, 7.43291099620453e-11, 2.45856059994029e-13,
        0, 0, 0, 0, 0, 0, 0, 0, 4.9795458121722e-10, 7.02657705613401e-12,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4.38033095281154e-11),
        ncol = 10, nrow = 10), sparse = TRUE)

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = -1,
                 midbin_tresh = 12),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
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


test_that("make_IPM works : mid_bin", {

    ex <- new(
        "dtCMatrix",
        i = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
              1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
              9L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 4L, 5L, 6L, 7L, 8L, 9L, 5L, 6L,
              7L, 8L, 9L, 6L, 7L, 8L, 9L, 7L, 8L, 9L, 8L, 9L, 9L),
        p = c(0L, 10L, 19L, 27L, 34L, 40L, 45L, 49L, 52L, 54L, 55L),
        Dim = c(10L, 10L), Dimnames = list(NULL, NULL),
        x = c(0.987873888381782, 0.00810671476198199,
              3.47645904163991e-09, 5.54602476538689e-12, 6.72554426783967e-14,
              2.19529536203141e-15, 1.30524600743269e-16, 1.16570526768275e-17,
              1.40012822063037e-18, 2.10922830083143e-19, 0.978393620548668,
              0.0206746151850283, 4.34963654797835e-08, 1.15218819586787e-10,
              1.91534946496806e-12, 7.84619812699488e-14, 5.56452471238397e-15,
              5.73609694155996e-16, 7.7721570993954e-17, 0.973929927604643,
              0.0256548548997383, 8.43110755096776e-08, 2.6908455692255e-10,
              5.09027510741877e-12, 2.29861615470281e-13, 1.76102640475857e-14,
              1.93409916768572e-15, 0.976003638440618, 0.0235792066477467,
              6.84054478307206e-08, 2.21118320635311e-10, 4.28644724986862e-12,
              1.98007424483321e-13, 1.54658100358733e-14, 0.981363206203981,
              0.0179425725107829, 3.24616703937358e-08, 9.52830468207385e-11,
              1.76998615392203e-12, 7.96661915292784e-14, 0.986926018814668,
              0.0116991045130483, 1.04190612864791e-08, 2.57820429195632e-11,
              4.38566813030585e-13, 0.9907343114035, 0.00659142737623629, 2.43101061063417e-09,
              4.80606895808284e-12, 0.992316884209617, 0.00320337175432724,
              4.28608594233547e-10, 0.993374221687891, 0.00133502553507573,
              0.998859808152198), uplo = "L", diag = "N")

    expect_equal(
        make_IPM(species, climate, clim_lab = "test", fit, mesh, BA = 1,
                 correction = "none", level = 420, diag_tresh = 1500,
                 midbin_tresh = 0),
        validate_ipm(new_ipm(
            IPM = list(`1` = ex), BA = 1, mesh = seq(90, 1500, length.out = 10),
            climatic = climate, clim_lab = "test",
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

