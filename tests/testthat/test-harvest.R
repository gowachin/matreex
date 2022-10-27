
# Uneven ####
test_that("getPcutUneven works", {
  x <- c(0.2, 0.1, 2, 0.3)
  h <- c(0, 0.2, 0.3, 1)
  ct <- c(10, 20, 30, 40)
  hmax <- 1

  expect_equal(getPcutUneven(x, h, ct),
               c(0, 0, 0, 0))

  expect_equal(getPcutUneven(x, h, ct, targetBAcut = 20),
               c(0, 0.066155014, 0.131128762, 1))
})


test_that("getBAcutTarget works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)

    BA <- drop(x %*% ct)

    expect_equal(getBAcutTarget(BA), 19)
    expect_equal(getBAcutTarget(BA, targetBA = 80), 0)
})

test_that("uneven_harv works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)
    species <- list(IPM = list(mesh = ct),
                    harv_lim = c(dth = 15, dha = 35, hmax = 1))

    expect_equal(Uneven_harv(x, species, targetBAcut = 20, ct = ct),
                 c(0, 0.00000607, 0.26666331, 0.3))
    expect_equal(Uneven_harv(x, species, targetBAcut = 0,  ct = ct),
                 c(0, 0, 0, 0))
})

test_that("getBAstand works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)
    species <- list(IPM = list(mesh = ct),
                    harv_lim = c(dth = 15, dha = 35, hmax = 1))

    expect_equal(getBAstand(x, species, SurfEch = 0.03),
                 0.060737458)
})


# Even ####
test_that("RDI works", {
    x <- c(12, 31, 52, 60)
    tx <- sum(x)
    meshcm2 <- ((1:4)*10)^2
    RDI_int <- 15
    RDI_slo <- - 2

    expect_equal(RDI(x, RDI_int, RDI_slo, meshcm2, tx), 0.047843123)
})

test_that("getPcutEven works", {
    x <- c(1.2, 3.1, 5.2, 6.0)
    mesh <- 1:4
    RDIcoef <- c(intercept = 15, slope = - 2)
    targetRDI=0.65
    targetKg=0.9
    SurfEch = 0.03

    getPcutEven(x, mesh, RDIcoef)

    expect_equal(getPcutEven(x, mesh, RDIcoef),  rep(0, 4))

    x <- c(1120, 310, 520, 600)
    mesh <- 25:28
    expect_equal(getPcutEven(x, mesh, RDIcoef, targetRDI = 0.1),
                 c(0.41422628, 0.20711314, 0.13807543, 0.10355657))
})


test_that("even_harv works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)
    species <- list(IPM = list(mesh = ct),
                    harv_lim = c(dth = 15, dha = 35, hmax = 1),
                    rdi_coef = c(intercept = 15, slope = - 2) )
    expect_equal(Even_harv(x, species, targetRDI = 0.1, targetKg = 0.9,
                             ct = ct, SurfEch = 0.03),
                 rep(0, 4))

    x <- c(1120, 310, 520, 600)
    expect_equal(Even_harv(x, species, targetRDI = 0.1, targetKg = 0.9,
                           ct = ct, SurfEch = 0.03),
                 c(307.378325, 74.427492, 115.451172, 126.019981))
})

