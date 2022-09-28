test_that("sim_deter_forest simple", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    set.seed(666)
    res <- sim_deter_forest(Forest = model, tlim = 30, equil_time = 1e3,
                            correction = "cut")

    expect_equal(dim(res), c(63, 31))
    expect_equal(colnames(res), c(paste0("t", 1:30), "eqt260"))
})


test_that("sim_deter_species simple", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2species("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    set.seed(666)
    res <- sim_deter_forest(Forest = model, tlim = 30, equil_time = 1e3,
                            correction = "cut")

    expect_equal(dim(res), c(63, 31))
    expect_equal(colnames(res), c(paste0("t", 1:30), "eqt260"))
})


test_that("sim_deter_forest delay & cut", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2species("Yggdrasil", climatic = 1, path = path,
                            replicat = 1, delay = 1)

    res<-evaluate_promise({
        set.seed(42)
        new <- sim_deter_forest(Forest = model, tlim = 500, equil_time = 1e3,
                                correction = "cut",  verbose = TRUE)
    })

    expect_equal(res$messages[1], "apply a IPM cut correction\n")
    expect_equal(res$messages[2], "Starting while loop. Maximum t = 1000\n")
    expect_equal(res$messages[3], "time 500 | BA diff : 0.00\n")
    expect_equal(res$messages[4], "Simulation ended after time 500\n")
    expect_equal(res$messages[5], "BA stabilized at 2.38 with diff of 0.00 at time 500\n")

    expect_equal(dim(new), c(65, 501))
    expect_equal(colnames(new), c(paste0("t", 1:500), "eqt500"))
})


test_that("sim_deter_forest delay & cut", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2species("Yggdrasil", climatic = 1, path = path,
                            replicat = 1, delay = 5)

    res<-evaluate_promise({
        set.seed(666)
        new <- sim_deter_forest(Forest = model, tlim = 10, equil_time = 1e3,
                                correction = "cut",  verbose = TRUE)
    })

    expect_equal(res$messages[1], "apply a IPM cut correction\n")
    expect_equal(res$messages[2], "Starting while loop. Maximum t = 1000\n")
    expect_equal(res$messages[3], "Simulation ended after time 4\n")
    expect_equal(res$messages[4], "BA stabilized at 1.07 with diff of 0.08 at time 3\n")
    expect_equal(res$warnings[1], "Maximum Basal Area reached for this simulation.")

    expect_equal(dim(new), c(73, 11))
    expect_equal(colnames(new), c(paste0("t", 1:10), "eqt3"))
})


test_that("sim_deter_forest error", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    model$species$Yggdrasil$init_pop <- function(mesh, SurfEch = 0.03){
        ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
        ini <- exp(runif(1, -.005, .005) * mesh)
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        while(all(alea)){ # because god knows it's fucking possible.
            # and it will return NaN
            alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
        }
        ini[alea] <- 0
        res <- as.numeric(ini / sum(ct * ini) )
        res <- res * 300 # HACK to limit falling in floating point trap !

        return(res)
    }

    exp <- paste(
        "Maximum Basal Area reached for this simulation.",
        "This maximum is reached before iteration, check init_pop functions"
    )
    expect_error(
        sim_deter_forest(Forest = model, tlim = 500, equil_time = 1e3),
        exp)
})


test_that("sim_deter_forest error equil_time", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    err <- "equil_time must be higher or equal to tlim and equil_dist"
    expect_error( sim_deter_forest(Forest = model, tlim = 30, equil_time = 10,
                                   correction = "cut")
                  , err)
})


test_that("sim_deter format work", {

    res <- structure(
        c(0.26, 0.00, 0.00, 0.27, 0.27, 1.00, 2.57, 0, 0, 0, 0, 0, 0.00,
          0.28, 0.37, 0.05, 0.05, 0.18, 1.04, 2.83, 0, 0, 0, 0, 0, 0.01,
          0.28, 0.44, 0.26, 0.13, 0.09, 1.10, 3.12, 0, 0, 0, 0, 0, 0.02,
          0.28, 0.45, 0.33, 0.27, 0.17, 1.20, 3.45, 0, 0, 0, 0, 0, 0.02,
          0.27, 0.44, 0.35, 0.32, 0.26, 1.27, 3.73, 0, 0, 0, 0, 0, 0.02),
        dim = c(13L, 5L),
        dimnames = list( c(
            "Ygg.m1", "Ygg.m2", "Ygg.m3", "Ygg.m4", "Ygg.m5", "Ygg.BAsp", "Ygg.N",
            "Ygg.h1", "Ygg.h2", "Ygg.h3", "Ygg.h4", "Ygg.h5", "Ygg.H"
        ), c("t1", "t2", "t3", "t4", "t5")))
    res <- new_deter_sim(res)

    exp <- structure(list(
        species = c("Ygg", "Ygg", "Ygg", "Ygg", "Ygg", "Ygg"),
        var = c("m", "m", "m", "m", "m", "m"),
        time = c(1, 2, 3, 4, 5, 1),
        mesh = c(1, 1, 1, 1, 1, 2),
        size = rep(NA, 6),
        equil = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        value = c(0.26, 0.28, 0.28, 0.28, 0.27, 0)),
        row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"))


    done <-evaluate_promise({
        new <- head(tree_format(res))
    })

    expect_equal(done$warnings[1], "mesh attribute missing, size column will be composed of NA")
    expect_identical(new, exp)


    attributes(res)$mesh <- list(Ygg = 1:5/10)
    exp$size <- c(.1, .1, .1, .1, .1, .2)
    expect_identical(head(tree_format(res)), exp)


})

