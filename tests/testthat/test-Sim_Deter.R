test_that("sim_deter_forest simple", {

    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    new <- sim_deter_forest(model, tlim = 20,
                            Xini = ini, SurfEch = SurfEch,
                            correction = "cut")

    expect_equal(new["BA", ncol(new)],  2.65286225)
    expect_equal(
        new["N", ],
        c(t1 = 9.5477662579656, t2 = 10.5242627967232, t3 = 11.7290585045797,
          t4 = 12.7890419188541, t5 = 13.7426339419704, t6 = 14.6558088273717,
          t7 = 15.5216561526095, t8 = 16.3363557372735, t9 = 17.1229282438257,
          t10 = 17.9101235706325, t11 = 18.7157691905626, t12 = 19.5384866494996,
          t13 = 20.3599432658585, t14 = 21.1554345998797, t15 = 21.9031537927613,
          t16 = 22.5867184312445, t17 = 23.1932248186191, t18 = 23.7110994623583,
          t19 = 24.1298017815526, t20 = 24.441350377244, t262 = 23.0754436306797
        )
    )
})


test_that("sim_deter_forest delay", {

    delay <- 1
    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))
    ct <- delay.numeric(ct, delay)

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- delay.numeric(ini, delay)
    ini <- as.numeric(ini / sum(ct * ini))

    Niter <- 20
    res<-evaluate_promise(
        new <- sim_deter_forest(model, tlim = Niter,
                         Xini = ini, SurfEch = SurfEch, delay = delay,
                         correction = "cut", verbose = TRUE)
    )

    expect_equal(new["BA", ncol(new)],  2.6696178)

    expect_equal(res$messages[1], "apply a IPM delay of 1\n")
    expect_equal(res$messages[4], "Simulation ended after time 264\n")
    expect_equal(res$messages[5], "BA stabilized at 2.67 with diff of 1.05\n")
})


test_that("sim_deter_forest error delay", {

    delay <- 1
    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    Niter <- 20
    err <- paste("Length of population \\(30\\) differs from IPM size \\(31\\),",
                 "please check dimensions. Is Xini delayed ?")
    expect_error(sim_deter_forest(model, tlim = Niter,
                                Xini = ini, SurfEch = SurfEch, delay = delay)
                 , err)
})


test_that("sim_deter_forest error equil_time", {

    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    Niter <- 20
    err <- "equil_time must be higher or equal to tlim and equil_dist"
    expect_error( sim_deter_forest(model, tlim = 400,
                                   Xini = ini, SurfEch = SurfEch,
                                   equil_time = 5)
                  , err)
})


test_that("sim_deter_forest message", {

    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    res<-evaluate_promise(
        sim_deter_forest(model, tlim = 20,
                         Xini = ini, SurfEch = SurfEch,
                         correction = "cut", verbose = TRUE)
    )

    expect_equal(res$messages[1], "apply a IPM cut correction\n")
    expect_equal(res$messages[2], "Starting while loop. Maximum t = 10000\n")
    expect_equal(res$messages[3], "Simulation ended after time 263\n")
    expect_equal(res$messages[4], "BA stabilized at 2.65 with diff of 1.03\n")
})


test_that("sim_deter_forest warning BA", {

    delay <- 5
    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))
    ct <- delay.numeric(ct, delay)

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- delay.numeric(ini, delay)
    ini <- as.numeric(ini / sum(ct * ini))

    Niter <- 20

    res<-evaluate_promise(
        sim_deter_forest(model, tlim = Niter,
                         Xini = ini, SurfEch = SurfEch, delay = delay,
                         correction = "cut", verbose = TRUE)
    )

    expect_equal(res$warnings[1], "Maximum Basal Area reached for this simulation.")
    expect_equal(res$messages[1], "apply a IPM delay of 5\n")
    expect_equal(res$messages[4], "Simulation ended after time 26\n")
    expect_equal(res$messages[5], "BA stabilized at 3.02 with diff of NA\n")
})


test_that("sim_deter_forest warning BA", {

    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    Niter <- 600

    res<-evaluate_promise(
        sim_deter_forest(model, tlim = Niter,
                         Xini = ini, SurfEch = SurfEch,
                         correction = "cut", verbose = TRUE)
    )

    expect_equal(res$messages[3], "time 500 | BA diff : 0.00\n")
    expect_equal(res$messages[4], "Simulation ended after time 601\n")
    expect_equal(res$messages[5], "BA stabilized at 2.65 with diff of 0.00\n")
})


test_that("sim_deter_forest summary", {

    SurfEch <- 0.1
    model <- treeforce:::try_mini

    ct <- drop(Buildct(model$meshpts, SurfEch=SurfEch))

    set.seed(666)
    ini <- state_init(model$meshpts)
    ini <- as.numeric(ini / sum(ct * ini))

    sim <- sim_deter_forest(model, tlim = 2,
                            Xini = ini, SurfEch = SurfEch,
                            correction = "cut")
    x <- summary(sim)

    expect_equal(x, structure(
        list(
            BA = c(t1 = 1.04860292734326, t2 = 1.10049667990408, t2 = 1.10049667990408),
            N = c(t1 = 9.5477662579656, t2 = 10.5242627967232, t2 = 10.5242627967232),
            state_eq = c(m1 = 1.00202124323669, m2 = 1.60161404109516,
                         m3 = 0.895476142635547, m4 = 0.371315666761291, m5 = 0.321499602934222,
                         m6 = 0.539933183152098, m7 = 0.534341341126861, m8 = 0.289016672591804,
                         m9 = 0.142082669585628, m10 = 0.197888154801501, m11 = 0.333818632219745,
                         m12 = 0.23370425338651, m13 = 0.111065089210931, m14 = 0.062983774507838,
                         m15 = 0.143920512491159, m16 = 0.321429094415444, m17 = 0.367386004967978,
                         m18 = 0.426699918133126, m19 = 0.329743840414589, m20 = 0.274267182821563,
                         m21 = 0.359101311242012, m22 = 0.300085502381768, m23 = 0.184901711783716,
                         m24 = 0.174937851178927, m25 = 0.302084096535827, m26 = 0.290040414992235,
                         m27 = 0.185681033611153, m28 = 0.105688091683625, m29 = 0.121535762824208,
                         m30 = 0),
            BA_eq = c(t2 = 1.10049667990408),
            time_eq = "2"),
        class = "summary_sim"))

    res<-evaluate_promise(print(sim))

    expect_equal(res$output, paste0(
        "Summary of deterministic simulation :\n",
        "Equilibrium reached after time 2 for an initial duration of 2 times\n",
        "BA stabilized at 1.100497 \n",
        "Mesh dimension : 30"
    ))

    res<-evaluate_promise(print(x))

    expect_equal(res$output, paste0(
        "Summary of deterministic simulation :\n",
        "Equilibrium reached after time 2 for an initial duration of 2 times\n",
        "BA stabilized at 1.100497 \n",
        "Mesh dimension : 30"
        ))
})


test_that("tree_format deter_sim", {
    library(tibble)
    x <- structure(matrix(
        c(1,3,4,2, 2, 4, 6, 3, 6, 8, 14, 7), nrow = 4,
        dimnames = list(c("m1", "m2", "N", "BA"), c("t1", "t2", "t7"))
        ), class = "deter_sim")

    res <- tibble(
        name = c("state_eq", "state_eq", "BA", "N", "BA", "N", "BAeq", "Neq"),
        t = rep(c(NA, 1, 2, 7), each = 2), m = c(1, 2, rep(NA, 6)),
        value = c(6, 8, 2, 4, 3, 6, 7, 14)
        )

    expect_equal(tree_format(x), res)
})

