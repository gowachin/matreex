test_that("new_species works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                  "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file) # NOTE 10" to load...
    raw_IPM <- raw_IPM[[1]]
    IPM <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)

    expect_identical(new_species(IPM = IPM, init_pop = def_init,
                                 harvest_fun = def_harv,
                                 # recruit_fun = raw_IPM$RecFun
                                 ),
                     structure(list(
                         IPM = IPM, init_pop = def_init, harvest_fun = def_harv,
                         harv_lim = c(dth = 175, dha = 575, hmax = 1),
                         recruit_fun = exp_recFun(params = IPM$rec$params_m,
                                                  list_covs = IPM$climatic),
                         info = c(species = "Yggdrasil", clim_lab = "1")),
                         class = "species"))
})


test_that("validate_species works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file) # NOTE 10" to load...
    raw_IPM <- raw_IPM[[1]]
    IPM <- old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1)

    x <- new_species(IPM, def_init, def_harv,
                     # recruit_fun = exp_recFun(params = raw_IPM$rec$params_m,
                                # list_covs = raw_IPM$list_m)
    )


    expect_identical(x, validate_species(x))
    tmp <- x
    names(tmp) <- c("IPM", "init_pop", "harvest_func", "recruit_fun", "info")
    expect_error(
        validate_species(tmp),
        "species class must be composed of elements IPM, init_pop, harvest_fun, recruit_fun and info"
    )
    tmp <- x
    names(tmp$info) <- c("sp", "clim_lab")
    expect_error(
        validate_species(tmp),
        "species class must have info of elements species and clim_lab"
    )
})


test_that("old_ipm2species works", {

    path <- here(ifelse(interactive() | covr::in_covr(), "tests", ""),
                 "testthat", "testdata")

    file <- here(path, "output", "Yggdrasil", "IPM_Clim_1.Rds")
    raw_IPM <- readRDS(file) # NOTE 10" to load...
    raw_IPM <- raw_IPM[[1]]

    expect_identical(
        old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1),
        new_species(
            old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
            def_init, def_harv,
            # recruit_fun = exp_recFun(params = raw_IPM$rec$params_m,
                                           # list_covs = raw_IPM$list_m)
            )
    )

    expect_identical(
        old_ipm2species("Yggdrasil", climatic = 1, path = path, replicat = 1, delay = 2),
        new_species(
                delay(old_ipm2ipm("Yggdrasil", climatic = 1, path = path, replicat = 1),
                      delay = 2
                      ),
            def_init, def_harv,
            # recruit_fun = exp_recFun(params = raw_IPM$rec$params_m,
                                     # list_covs = raw_IPM$list_m)
        )
    )

})


test_that("def_init works", {

    x <- seq(30, 40, by = 2)

    set.seed(666)
    expect_identical(
        def_init(x),
        c(10.34167240524874, 1e-04, 1e-04, 10.513325511364572, 10.571174040583468,
          1e-04)
    )

    # NOTE Case with all(alea) require specific seed that I don't know.
    foo <- def_initBA(20)
    set.seed(666)
    expect_identical(
        foo(x),
        c(206.8314481049748, 0, 0, 210.26451022729145, 211.42148081166937,
          0)
    )

    foo <- def_init_k(c(12, 2, 0, 40, 0, 5))
    expect_identical(foo(1:6),c(12, 2, 0, 40, 0, 5))
    # FIXME can't test this !
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

