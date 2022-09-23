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

