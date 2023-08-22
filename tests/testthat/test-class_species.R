test_that("new_species works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                  "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file)
    raw_IPM <- raw_IPM[[1]]
    IPM <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)

    expect_identical(new_species(IPM = IPM, init_pop = def_init,
                                 harvest_fun = def_harv,
                                 disturb_fun = def_disturb
                                 ),
                     structure(list(
                         IPM = IPM, init_pop = def_init, harvest_fun = def_harv,
                         harv_lim = c(dth = 175, dha = 575, hmax = 1),
                         disturb_fun = def_disturb,
                         rdi_coef = NULL, disturb_coef = NULL,
                         recruit_fun = exp_recFun(params = IPM$fit$rec$params_m,
                                                  list_covs = IPM$climatic),
                         info = c(species = "Yggdrasil", clim_lab = "1", type = "Broadleaf")),
                         class = "species"))

    class(IPM) <- "mu_gr"

    expect_identical(
        new_species(IPM = IPM, init_pop = def_init, harvest_fun = def_harv,
                    disturb_fun = def_disturb),
        structure(list(
            IPM = IPM, init_pop = def_init, harvest_fun = def_harv,
            harv_lim = c(dth = 175, dha = 575, hmax = 1),
            disturb_fun = def_disturb,
            rdi_coef = NULL, disturb_coef = NULL, recruit_fun = "to define",
            info = c(species = "Yggdrasil", clim_lab = "1", type = "Broadleaf")),
            class = "species")
    )

    class(IPM) <- "mu_growth"

    expect_error(
        new_species(IPM = IPM, init_pop = def_init, harvest_fun = def_harv,
                    disturb_fun = def_disturb),
        "IPM must either be an ipm or mu_gr object."
    )
})


test_that("validate_species ipm works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                 "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file)
    raw_IPM <- raw_IPM[[1]]
    IPM <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)

    x <- new_species(IPM, def_init, def_harv, disturb_fun = def_disturb)


    expect_identical(x, validate_species(x))
    tmp <- x
    names(tmp) <- c("IPM", "init_pop", "harvest_func", "harv_lim", "rdi_coef",
                    "recruit_fun", "info")
    expect_error(
        validate_species(tmp),
        "species class must be composed of elements IPM, init_pop, harvest_fun, harv_lim, disturb_fun, rdi_coef, disturb_coef, recruit_fun and info"
    )
    tmp <- x
    names(tmp$info) <- c("sp", "clim_lab")
    expect_error(
        validate_species(tmp),
        "species class must have info of elements species, clim_lab and type"
    )

    class(x$IPM) <- "mu_growth"
    expect_error(
        validate_species(x),
        "IPM must either be an ipm or mu_gr object."
    )


})



test_that("validate_species mu works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                 "testthat", "testdata")

    mu <- make_mu_gr(species = "Picea_abies", fit = fit_Picea_abies,
                     mesh = c(m = 10, L = 90, U = 1200), stepMu = 1,
                     level = c(3, 10), midbin_tresh = 2)
    x <- species(mu, init_pop = def_init)


    expect_identical(x, validate_species(x))


})


test_that("old_ipm2species works", {

    path <- here(ifelse(testthat:::on_ci() | interactive() | covr::in_covr(),
                        "tests", "."),
                 "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file)
    raw_IPM <- raw_IPM[[1]]

    expect_identical(
        old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_species(
            old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
            def_init, def_harv, disturb_fun = def_disturb,
            rdi_coef = c(intercept = 13.99, slope = -2.18),
            disturb_coef = data.frame(
                disturbance = "biotic", species = "Yggdrasil",
                a0 = -5.81, a1 = 0, b = 2.9, c = 0.0052,
                dbh.intercept = -0.787, dbh.slope = 0.00793,
                logratio.intercept = 0.468, logratio.slope = 2.92, row.names = 64L)
            )
    )

    expect_identical(
        old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1, delay = 2),
        new_species(
                delay(old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
                      delay = 2
                      ),
            def_init, def_harv, disturb_fun = def_disturb,
            rdi_coef = c(intercept = 13.99, slope = -2.18),
            disturb_coef = data.frame(
                disturbance = "biotic", species = "Yggdrasil",
                a0 = -5.81, a1 = 0, b = 2.9, c = 0.0052,
                dbh.intercept = -0.787, dbh.slope = 0.00793,
                logratio.intercept = 0.468, logratio.slope = 2.92, row.names = 64L)
        )
    )

})


test_that("def_init works", {

    x <- seq(30, 40, by = 2)

    set.seed(666)
    expect_equal(
        def_init(x),
        c(10.56114746134755, 1e-04, 1e-04, 10.4850290591089, 10.4597783699702,
          1e-04)
    )

    set.seed(666)
    expect_equal(
        def_init_even(x),
        c(6.4829468383972, 6.51861894865491, 6.5544873430226, 6.59055310154436,
          6.62681731020718, 1e-15)
    )

    set.seed(666)
    expect_equal(
        def_init_even(delay(x, 2)),
        c(0, 0, 6.4829468383972, 6.51861894865491, 6.5544873430226, 6.59055310154436,
          6.62681731020718, 1e-15)
    )

    # Case with all(alea) require specific seed that I don't know.
    foo <- def_initBA(20)
    set.seed(666)
    expect_equal(
        foo(x),
        c(211.220949226951, 0, 0, 209.698581182179, 209.193567399403,
          0)
    )

    foo <- def_init_k(c(12, 2, 0, 40, 0, 5))
    expect_identical(foo(1:6, SurfEch = 1),c(12, 2, 0, 40, 0, 5))
    expect_error(
        foo(1:5),
        paste0("A species initiate with a define distribution with ",
               "different length that it's mesh. Check sp$init_pop ",
               "functions using def_init_k !"),
        fixed = TRUE
    )
 })


test_that("def_harv works", {

    x <- 1:10
    ct <- rep(1, 10)

    expect_identical(
        def_harv(x, ct = ct),
        x * 0.006
    )

    ct <- c(0, 0, rep(1, 8))

    expect_identical(
        def_harv(x, ct = ct),
        x *( 0.006 * (ct > 0))
    )

})


test_that("def_disturb works", {

    x <- 1:10

    expect_identical(
        def_disturb(x),
        x * 0
    )
    expect_warning(
        def_disturb(x, disturb = "highway to hell"),
        "default disturbance function does not impact populations. Please add your own disturbance function."
    )

})

