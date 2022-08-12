test_that("new_ipm works", {

    place <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                  "testthat", "testdata")
    file <- here(place, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- 1
    compress <- TRUE

    expect_identical(new_ipm(IPM, BA, mesh, species, climatic, compress),
                     structure(list(
                         IPM = IPM, BA = BA, mesh = mesh,
                         info = c(species = "Yggdrasil", climatic = "1",
                                  compress = "TRUE")), class = "ipm"))
})


test_that("validate_ipm works", {

    place <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                  "testthat", "testdata")
    file <- here(place, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- 1
    compress <- TRUE
    x <- new_ipm(IPM, BA, mesh, species, climatic, compress)


    expect_identical(x, validate_ipm(x))
    tmp <- x
    names(tmp) <- c("IPM", "BA", "meshpts", "info")
    expect_error(
        validate_ipm(tmp),
        "IPM class must be composed of elements IPM, BA, mesh and info"
    )
    tmp <- x
    names(tmp$info) <- c("species", "climatic", "compressed")
    expect_error(
        validate_ipm(tmp),
        "IPM class must have info of elements species, climatic and compress"
    )
})


test_that("old_ipm2ipm works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")
    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    x <- readRDS(file)

    IPM <- x[[1]]$LIPM
    BA <- 1:length(IPM)
    mesh <- x[[1]]$meshpts
    species <- "Yggdrasil"
    climatic <- 1
    compress <- TRUE

    expect_identical(
        old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_ipm(IPM, BA, mesh, species, climatic, compress)
    )

})
