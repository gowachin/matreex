test_that("X2Pop works", {
  set.seed(42)
  expect_equal(X2Pop(c(1, 1, 2, 1.1, 0.5), c(0, 0, 90, 100, 110)), 
               c(-2, -2, -1, -1, -1, 90, 100, 100, 110)
              )
})

test_that(" new_indiv_sim works", {
  set.seed(42)
  expect_equal(new_indiv_sim(matrix(c(1, 1, 2, 2), ncol = 2), mesh = c(0, 0, 90, 100, 110)), 
               structure(c(1, 1, 2, 2), dim = c(2L, 2L), 
                         class = c("indiv_sim", "matrix"), 
                         mesh = c(0, 0, 90, 100, 110))
              )
})


