test_that("new_forest works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    Yggdrasil <- old_ipm2species("Yggdrasil", climatic = 1,
                                 path = path, replicat = 1)

    expect_identical(new_forest(list(Yggdrasil)),
                     structure(list(
                         species = list(Yggdrasil = Yggdrasil),
                         info = list(species = c("Yggdrasil"),
                                     climatic = c(Yggdrasil = 1))),
                         class = "forest"))

    expect_identical(
        new_forest(list(Yggdrasil, Yggdrasil)),
        structure(list(
            species = list(Yggdrasil = Yggdrasil,
                           Yggdrasil = Yggdrasil),
            info = list(species = c("Yggdrasil", "Yggdrasil"),
                        climatic = c(Yggdrasil = 1, Yggdrasil = 1))),
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
    tmp$info$climatic <- c(Yggdrasil = 1, Yggdrasil = 2)
    expect_error(
        validate_forest(tmp),
        "All species are not defined for the same climatic."
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

