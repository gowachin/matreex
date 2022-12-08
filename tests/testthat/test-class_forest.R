test_that("new_forest works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    Yggdrasil <- old_ipm2species("Yggdrasil", climatic = 1,
                                 path = path, replicat = 1)
    # browser()
    map_chr(list(Yggdrasil), climatic)
    expect_identical(new_forest(list(Yggdrasil)),
                     structure(list(
                         species = list(Yggdrasil = Yggdrasil),
                         harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                        freq = 1, alpha = 1),
                         info = list(species = c("Yggdrasil"),
                                     clim_lab = c(Yggdrasil = "1"))),
                         class = "forest"))

    expect_identical(
        new_forest(list(Yggdrasil, Yggdrasil)),
        structure(list(
            species = list(Yggdrasil = Yggdrasil,
                           Yggdrasil = Yggdrasil),
            harv_rules = c(Pmax = 0.25, dBAmin = 3,
                           freq = 1, alpha = 1),
            info = list(species = c("Yggdrasil", "Yggdrasil"),
                        clim_lab = c(Yggdrasil = "1", Yggdrasil = "1"))),
            class = "forest")
    )
})


test_that("validate_forest works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    Yggdrasil <- old_ipm2species("Yggdrasil", climatic = 1,
                                 path = path, replicat = 1)


    x <- forest(list(Yggdrasil, Yggdrasil))


    expect_identical(x, validate_forest(x))
    tmp <- x
    tmp$info$clim_lab <- c(Yggdrasil = 1, Yggdrasil = 2)
    expect_error(
        validate_forest(tmp),
        "Some ipm species are not defined with the same climatic name"
    )
})


test_that("old_ipm2species works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    Yggdrasil <- old_ipm2species("Yggdrasil", climatic = 1,
                                 path = path, replicat = 1)

    expect_identical(
        old_ipm2forest("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_forest(list(Yggdrasil))
    )

    expect_identical(
        old_ipm2forest(c("Yggdrasil", "Yggdrasil"), climatic = 1,
                       path = path, replicat = 1),
        new_forest(list(Yggdrasil, Yggdrasil))
    )

})

