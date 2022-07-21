test_that("Run_Sim_Deter single run works", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
                          # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                  climatic = 1, tlim = 5, correction = "cut")

    expect_equal(dim(x), c(210, 9))
    expect_equal(unique(x$name), c("state_eq", "BA", "N", "BAeq", "Neq"))
    expect_equal(unique(x$t), c(NA, 1:6))
    expect_equal(unique(x$m), c(1:30, NA))
    expect_equal(unique(x$delay), 0)
    expect_equal(unique(x$clim), 1)
    expect_equal(unique(x$harv), 0.006)
    expect_equal(unique(x$corr), c("cut"))
    expect_equal(unique(x$n), 1:5)
    expect_equal(range(x$value), c(0, 4.3653813))
})


test_that("Run_Sim_Deter single run talk", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    res<-evaluate_promise(
        x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                           climatic = 1, verbose = TRUE, tlim = 5,
                           correction = "cut")
    )

    expect_equal(res$messages[1], "proof done\n")

})


test_that("Run_Sim_Deter multi run works", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                       climatic = 1:2, tlim = 5, correction = "cut")

    expect_equal(dim(x), c(420, 9))
    expect_equal(unique(x$name), c("state_eq", "BA", "N", "BAeq", "Neq"))
    expect_equal(unique(x$t), c(NA, 1:6))
    expect_equal(unique(x$m), c(1:30, NA))
    expect_equal(unique(x$delay), 0)
    expect_equal(unique(x$clim), 1:2)
    expect_equal(unique(x$harv), 0.006)
    expect_equal(unique(x$corr), c("cut"))
    expect_equal(unique(x$n), 1:5)
    expect_equal(range(x$value), c(0, 4.4894449))
})


test_that("Run_Sim_Deter multi run talk", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    res<-evaluate_promise(
        x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                           climatic = 1:2, verbose = TRUE, tlim = 5,
                           correction = "cut")
    )

    expect_equal(res$messages[1], "proof done\n")
    expect_equal(res$messages[2], "Recursion happenning now !\n")
    expect_equal(res$messages[3], "proof done\n")
    expect_equal(res$messages[4], "proof done\n")

})


test_that("Run_Sim_Deter parallel run works", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                       climatic = 1:2, tlim = 5, correction = "cut",
                       parallel = TRUE)

    expect_equal(dim(x), c(420, 9))
    expect_equal(unique(x$name), c("state_eq", "BA", "N", "BAeq", "Neq"))
    expect_equal(unique(x$t), c(NA, 1:6))
    expect_equal(unique(x$m), c(1:30, NA))
    expect_equal(unique(x$delay), 0)
    expect_equal(unique(x$clim), 1:2)
    expect_equal(unique(x$harv), 0.006)
    expect_equal(unique(x$corr), c("cut"))
    expect_equal(unique(x$n), 1:5)
    # expect_equal(range(x$value), c(0, 4.981854))
    # TODO random seed in parallel makes it difficult to test... seed argument?
})


test_that("Run_Sim_Deter parallel run talk", {

    place <- here(ifelse(interactive(), "tests", ""),
                  "testthat", "testdata")

    init <- function(mesh, SurfEch = 0.03) {
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )

        return(res)
    }

    set.seed(666)
    res<-evaluate_promise(
        x <- run_sim_deter("Yggdrasil", init_pop = init, path = place,
                           climatic = 1:2, verbose = TRUE, tlim = 5,
                           parallel = TRUE,
                           correction = "cut")
    )

    expect_equal(res$messages[1], "proof done\n")
    expect_equal(res$messages[2], "Recursion happenning now !\n")
    expect_equal(res$messages[3], "going parallel on 2 cores.\n")

})

