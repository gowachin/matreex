
test_that("getPcutUneven works", {
  x <- c(0.2, 0.1, 2, 0.3)
  h <- c(0, 0.2, 0.3, 1)
  ct <- c(10, 20, 30, 40)
  hmax <- 1

  expect_equal(getPcutUneven(x, h, ct),
               c(0, 0.055392952, 0.114819476, 1))

  expect_equal(getPcutUneven(x, h, ct, targetBA = 80),
               c(0, 0, 0, 0))
})


test_that("getBAcutTarget works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)

    expect_equal(getBAcutTarget(x, ct), 19)
    expect_equal(getBAcutTarget(x, ct, targetBA = 80), 0)
})

test_that("uneven_harv works", {
    x <- c(0.2, 0.1, 2, 0.3)
    ct <- c(10, 20, 30, 40)
    species <- list(IPM = list(mesh = ct),
                    harv_lim = c(dth = 15, dha = 35, hmax = 1))
    harv_rule <- c(Pmax = 1, dBAmin = 3, freq = 10)

    expect_equal(Uneven_harv(x, species, harv_rule, targetBA = 20,
                             ct, t = 1),
                 c(0, 0, 0, 0))
    expect_equal(Uneven_harv(x, species, harv_rule, targetBA = 20,
                             ct, t = 10),
                 c(0, 0.01723069, 1.38850837, 0.3))
})

