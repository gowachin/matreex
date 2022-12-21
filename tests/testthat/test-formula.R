test_that("qmd works", {


    expect_equal(
        QMD(c(1:5), c(3, 2, 1, 4, 2)),
        3.3416563
    )

})
