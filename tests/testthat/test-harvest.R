
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
    meshcm2 <- ((1:4)*10)^2
    RDI_int <- 15
    RDI_slo <- - 2

    expect_equal(RDI(x, RDI_int, RDI_slo, meshcm2), 0.047843123)
})

test_that("getPcutEven works", {
    x <- list(sp1 = c(1.2, 3.1, 5.2, 6.0))
    mesh <- list(sp1 = 1:4)
    ct <- c(10, 20, 30, 40)
    targetRDI=0.65
    targetKg=0.9
    SurfEch = 0.03
    species <- list(sp1 = list(IPM = list(mesh = ct),
                               harv_lim = c(dth = 15, dha = 35, hmax = 1),
                               rdi_coef = c(intercept = 15, slope = - 2) ))

    getPcutEven(x, species, mesh)

    expect_equal(
        getPcutEven(x, species, mesh),
        list(sp1 = c(0.000516621828056597, 0.000368914027743261,
                     0.000302953390848195, .000263437494264849)))

    x <- list(sp1 = c(1120, 310, 520, 600))
    mesh <- list(sp1 = 25:28)
    expect_equal(getPcutEven(x, species, mesh, targetRDI = 0.1),
                 list(sp1 = c(0.742161821098047, 0.371080910549023,
                              0.247387273699349, 0.185540455274512)))
})

test_that("sim_rdikg works", {
    sim <- structure(list(
        species = c(
            "Picea_abies", "Picea_abies", "Picea_abies", "Picea_abies",
            "Abies_alba", "Abies_alba", "Abies_alba", "Abies_alba"
        ),
        var = c("n", "n", "h", "h", "n", "n", "h", "h"),
        time = c(90, 91, 90, 91, 90, 91, 90, 91),
        mesh = c(29, 29, 29, 29, 29, 29, 29, 29),
        size = c(113.086071428571, 113.086071428571, 113.086071428571, 113.086071428571,
                 91.0514285714286, 91.0514285714286, 91.0514285714286, 91.0514285714286),
        equil = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        value = c(1.05310592209055, 1.06290314862089, 0.319890215992106, 0,
                  3.50307832197179, 3.26533476705513, 1.06409094996663, 0)),
        row.names = c(NA, -8L), class = c("tbl_df", "tbl", "data.frame"
        ))

    expect_equal(
        sim_rdikg(sim),
        structure(list(
            species = c("All", "All",
                        "Abies_alba", "Abies_alba", "Picea_abies", "Picea_abies",
                        "All", "All", "All", "All", "All", "All"),
            time = c(90, 91, 90, 91, 90, 91, 90, 90, 90, 91, 91, 91),
            var = c("rdi", "rdi", "rdi", "rdi", "rdi", "rdi",
                    "Dg2", "Dgcut2", "Kg", "Dg2", "Dgcut2", "Kg"),
            value = c(0.000641478493937736, 0.000612878352467643,
                      0.000447943911740301, 0.000417543284579767, 0.000193534582197435, 0.000195335067887876,
                      93.3004239364435, 93.3004239364435, 1, 93.9497876295814, NaN, 0)),
            class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -12L)))
})
