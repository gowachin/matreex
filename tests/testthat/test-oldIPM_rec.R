test_that("multi works", {
    df <- data.frame(var1 = "BATOTSP", var2 = NA_character_, params = -0.0179,
                     value.x = 1, value.y = 1, K = -0.0179)
    expect_identical(
        multi(x = "BATOTSP", df = df),
        call2("<-", expr(BATOTSP_in), call2("*", -0.0179, expr(BATOTSP)))
    )

    df <- data.frame(var1 = "BATOTSP2", var2 = NA_character_, params = -0.0179,
                     value.x = 1, value.y = 1, K = -0.0179)
    expect_identical(
        multi(x = "BATOTSP2", df = df),
        call2("<-", expr(BATOTSP2_in), call2("*", -0.0179, expr(BATOTSP ^ 2)))
    )

    df <- data.frame(var1 = "logBATOTSP", var2 = NA_character_, params = -0.0179,
                     value.x = 1, value.y = 1, K = -0.0179)
    expect_identical(
        multi(x = "logBATOTSP", df = df),
        call2("<-", expr(logBATOTSP_in), call2("*", -0.0179, expr(log(BATOTSP))))
    )
})


test_that("format_fit works", {
    params <- c(intercept = -0.864, BATOTSP = -0.018, sgddb = 286.813,
                wai = -0.057, wai2 = 0.288
    )

    list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)

    exp <- data.frame(
        var1 = c("intercept", "BATOTSP", "sgddb", "wai","wai2"),
        var2 = NA_character_,
        params = c(-0.864, -0.018, 286.813,             -0.057, 0.288),
        value.x = c(1, 1, 0, -0.187, 0.34),
        value.y = c(1, 1, 1, 1, 1),
        K = c(-0.864, -0.018, 0, 0.010659, 0.097920))

    expect_identical( format_fit(params, list_covs), exp)
})


test_that("format_fit works", {
    params <- c(intercept = -0.864, BATOTSP = -0.018, sgddb = 286.813,
                wai = -0.057, wai2 = 0.288
    )

    list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)

    empty <- function (BATOTSP, BATOTNonSP, mesh, SurfEch = 0.003) {
        intercept <- -0.755421
        res <- 0
        BATOTSP_in <- -0.018 * BATOTSP
        res <- res + intercept
        res <- res + BATOTSP_in
        mesh <- length(mesh)
        distrib <- c(rep(1/2, 2), numeric(mesh - 2))
        final <- exp(res) * SurfEch/0.03 * distrib
        return(final)
    }
    body(empty)[[2]][[3]] <- -0.755421
    body(empty)[[4]][[3]][[2]] <- -0.018

    expect_identical(formals(exp_recFun(params, list_covs)),  formals(empty))
    expect_identical(body(exp_recFun(params, list_covs)),  body(empty))
})
