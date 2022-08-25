
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

    expect_equal(Uneven_harv(x, species, targetBAcut = 20, ct),
                 c(0, 0.00000607, 0.26666331, 0.3))
    expect_equal(Uneven_harv(x, species, targetBAcut = 0,  ct),
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

