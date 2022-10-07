
# sp_name ####
test_that("sp_name works", {
    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    x <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(sp_name(x), "Yggdrasil")
    x <- old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(sp_name(x), "Yggdrasil")
})

test_that("climatic works", {
    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    x <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(climatic(x), "1")
    x <- old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1)
    expect_identical(climatic(x), "1")
})



# delay ####
test_that("delay dtCMatrix works", {
    x <- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L), p = c(0L, 3L,  5L, 6L),
             Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")

    expect_identical(
        delay(x, 2),
        new("dgCMatrix", i = c(1L, 2L,
                               2L, 3L, 4L, 3L, 4L, 4L),
            p = c(0L, 1L, 2L, 5L,  7L, 8L),
            Dim = c(5L, 5L), x = c(1, 1,
                                   1, 2, 3, 1, 2, 1))
    )

    expect_identical( delay(x, 0), x )
})


test_that("delay numeric works", {
    x <- as.numeric(1:10)

    expect_identical( delay(x, 2), c(0,0, x) )
    expect_identical( delay(x, 0), x )
})


test_that("delay ipm works", {

    rec <-  c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    mat<- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
              p = c(0L, 3L,  5L, 6L),
              Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec,
                 compress = FALSE)

    exp <- new_ipm(IPM = list( delay(mat, 2) ), BA = 1, mesh = c(0,0,1:3),
                   species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec,
                   compress = FALSE, delay = 2)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( delay(x, 2), exp )
    expect_identical( delay(x, 0), x )

    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec,
                 compress = TRUE)

    exp <- new_ipm(IPM = list( delay(mat * 1e-7, 2) ), BA = 1, mesh = c(0,0,1:3),
                   species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec,
                   compress = FALSE, delay = 2)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( delay(x, 2), exp )
})


test_that("delay species works", {

    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin", climatic = 1, clim_lab = "1",
    rec_params = c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1),
    compress = FALSE)
    validate_ipm(IPM)

    x <- new_species(IPM = IPM, init_pop = def_init,
                     harvest_fun = def_harv)

    exp <- new_species(delay(IPM, 2), init_pop = def_init,
                       harvest_fun = def_harv)
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
    species = "darwin", climatic = c(sgdd = 1), clim_lab = "1",
    rec_params = c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1),
    compress = FALSE)
    validate_ipm(IPM)

    sp <- new_species(IPM = IPM, init_pop = def_init,
                     harvest_fun = def_harv)

    x <- new_forest(list(darwin = sp))
    exp <- new_forest(list(darwin = delay(sp, 2)))
    # HACK ne pas nommer les elements ici fout la merde.

    # validate_forest(x)
    # validate_forest(exp)

    expect_identical( delay(x, 2), exp )
    expect_identical( delay(x, 0), x )

})


# correction ####
test_that("correction ipm works", {

    rec <-  c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    mat<- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
              p = c(0L, 3L,  5L, 6L),
              Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    new_mat<- new("dtCMatrix", i = c(0L, 1L, 1L),
              p = c(0L, 2L,  3L, 3L),
              Dim = c(3L, 3L), x = c(1, 2, 1), uplo = "L", diag = "N")
    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec, compress = FALSE)

    exp <- new_ipm(IPM = list( new_mat  ), BA = 1, mesh = 1:3,
                   species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec, compress = FALSE)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( correction(x, "cut"), exp )
    expect_identical( correction(x), x )

    x <- new_ipm(IPM = list(mat), BA = 1, mesh = 1:3,
                 species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec, compress = TRUE)

    exp <- new_ipm(IPM = list( new_mat * 1e-7 ), BA = 1, mesh = 1:3,
                   species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec, compress = FALSE)
    # validate_ipm(x)
    # validate_ipm(exp)

    expect_identical( correction(x, "cut"), exp )
})


test_that("correction species works", {

    rec <-  c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1)
    IPM <- new_ipm(IPM = list(
        new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L),
            p = c(0L, 3L,  5L, 6L),
            Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
    ), BA = 1, mesh = 1:3,
    species = "darwin", climatic = 1, clim_lab = "1", rec_params = rec, compress = FALSE)
    validate_ipm(IPM)


    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv)

    exp <- new_species(correction(IPM, "cut"), init_pop = def_init,
                       harvest_fun = def_harv)
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
    species = "darwin", climatic = 1, clim_lab = "1",
    rec_params = c(intercept = 2, BATOTSP = 0, BATOTNonSP = 1), compress = FALSE)
    validate_ipm(IPM)

    sp <- new_species(IPM = IPM, init_pop = def_init,
                      harvest_fun = def_harv)

    x <- new_forest(list(darwin = sp))
    exp <- new_forest(list(darwin = correction(sp, "cut")))
    # HACK ne pas nommer les elements ici fout la merde.

    expect_identical( correction(x, "cut"), exp )
    expect_identical( correction(x), x )

})
