test_that("utils works", {

    matrix <- matrix(0, ncol = 6, nrow = 6)
    N_int <- 2
    fill <- matrix(1, ncol = 6, nrow = N_int)
    fill[,2:4 ] <- 1:(3*N_int)

    exp <- structure(c(1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 4, 0,
                       0, 0, 0, 0, 5, 6, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1),
                     dim = c(6L, 6L))

    expect_identical(sub_diag(matrix, fill, dist = 0), exp)

    exp <- structure(c(0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3,
                       4, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                     dim = c(6L, 6L))

    expect_identical(sub_diag(matrix, fill, dist = 2), exp)


    expect_identical(sub_diag(matrix, fill, dist = 6), matrix)

    N_int <- 3
    fill <- matrix(1, ncol = 6, nrow = N_int)
    fill[,2:4 ] <- 1:(3*N_int)

    exp <- structure(c(1, 1, 1, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 4, 5, 6,
                       0, 0, 0, 0, 7, 8, 9, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1),
                     dim = c(6L, 6L))

    expect_identical(sub_diag(matrix, fill, dist = 0), exp)

    exp <- structure(c(0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 4,
                       5, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                     dim = c(6L, 6L))

    expect_identical(sub_diag(matrix, fill, dist = 2), exp)
})
