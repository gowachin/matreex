test_that("sim_deter_forest simple", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    set.seed(666)
    res <- sim_deter_forest(Forest = model, tlim = 30, equil_time = 1e3,
                            correction = "cut")

    expect_equal(dim(res), c(63, 31))
    expect_equal(colnames(res), paste0("t", c(1:30, 260)))
})


test_that("sim_deter_forest simple", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2species("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    set.seed(666)
    res <- sim_deter_forest(Forest = model, tlim = 30, equil_time = 1e3,
                            correction = "cut")

    expect_equal(dim(res), c(63, 31))
    expect_equal(colnames(res), paste0("t", c(1:30, 260)))
})


test_that("sim_deter_forest delay & cut", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    res<-evaluate_promise({
        set.seed(42)
        new <- sim_deter_forest(Forest = model, tlim = 500, equil_time = 1e3,
                                correction = "cut",  verbose = TRUE, delay = 1)
    })

    expect_equal(res$messages[1], "apply a IPM delay of 1\n")
    expect_equal(res$messages[2], "apply a IPM cut correction\n")
    expect_equal(res$messages[3], "Starting while loop. Maximum t = 1000\n")
    expect_equal(res$messages[4], "time 500 | BA diff : 0.00\n")
    expect_equal(res$messages[5], "Simulation ended after time 500\n")
    expect_equal(res$messages[6], "BA stabilized at 2.38 with diff of 0.00 at time 500\n")

    expect_equal(dim(new), c(65, 501))
    expect_equal(colnames(new), paste0("t", c(1:500, 500)))
})


test_that("sim_deter_forest delay & cut", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    model <- old_ipm2forest("Yggdrasil", climatic = 1, path = path,
                            replicat = 1)

    res<-evaluate_promise({
        set.seed(666)
        new <- sim_deter_forest(Forest = model, tlim = 10, equil_time = 1e3,
                                correction = "cut",  verbose = TRUE, delay = 5)
    })

    expect_equal(res$messages[1], "apply a IPM delay of 5\n")
    expect_equal(res$messages[2], "apply a IPM cut correction\n")
    expect_equal(res$messages[3], "Starting while loop. Maximum t = 1000\n")
    expect_equal(res$messages[4], "Simulation ended after time 4\n")
    expect_equal(res$messages[5], "BA stabilized at 1.07 with diff of 0.08 at time 3\n")
    expect_equal(res$warnings[1], "Maximum Basal Area reached for this simulation.")

    expect_equal(dim(new), c(73, 11))
    expect_equal(colnames(new), paste0("t", c(1:10, 3)))
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
