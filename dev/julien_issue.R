sim.equil.in = readRDS("tryhard/sim_equil.rds")
climate.in = readRDS("tryhard/climate.rds")
species.in = readRDS("tryhard/species.rds")

library(dplyr)
devtools::load_all()

distrib.in <- filter(sim.equil.in, equil, var == "n") %>% pull(value) * 0.03

species.in$init_pop <- def_init_k(distrib.in*0.03)
species.in$info["type"] <- "Broadleaf"


range(species.in$IPM$mu_tab)
range( mu_growth)

species.in$IPM$info

forest.in = forest(species = list(mu = species.in), harv_rules = c(
    Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))

sim.dist.in = sim_deter_forest(
    forest.in, tlim = max(climate.in$t), climate = climate.in,
    equil_dist = max(climate.in$t),  equil_time = max(climate.in$t),
    verbose = TRUE, correction = "cut")

# range during simulation
# range( mu_growth)
# [1] -5.017209  1.743862

# fit <- fit_Abies_alba
# climate <- rename(climate.in, N = t)

getRangemu_dev <- function(climate,
                       fit,
                       BA = 0:200,
                       mesh = seq(90, 900, by = 1)) {

    assertIntegerish(BA, lower = 0, upper = 200)
    assertNumeric(mesh, lower = 0)

    climate_species <- climate
    N <- NULL # hack to bind global value
    n <- nrow(climate)

    fres <- data.frame(min = 1:n, max = 1:n)
    for (Nc in 1:n) { # TODO replace climate_species for climate, because it's confusing
        climate <- subset(climate_species, N == Nc, select = -N)
        climate <- drop(as.matrix(climate)) # we need it as a vector.
        list_covs <- c(climate, BATOTcomp = 0)
        res <- matrix(
            ncol = 2, nrow = length(BA),
            dimnames = list(NULL, c("min", "max"))
        )

        for (iBA in seq_along(BA)) {
            list_covs["BATOTcomp"] <- BA[iBA]

            grFun <- exp_sizeFun(fit$gr$params_m, list_covs)
            mu <- grFun(mesh)

            res[iBA, ] <- range(mu)
        }
        fres[Nc, ] <- c(min(res[, "min"]), max(res[, "max"]))
    }

    range <- c(min = min(fres$min), max = max(fres$max), sig = fit$gr$sigma)
    print(range)
    return(fres)
    return(range)
}

getRangemu_dev(climate, fit_Abies_alba, BA = 0:200, mesh = species.in$IPM$mesh)

climate.in
climate
plot(climate$sgdd, climate$wai, cex = climate$N/ 100)
climate <- rename(climate.in, N = t)
tmp <- cbind(rbind(
    climate[climate$sgdd == min(climate$sgdd), c("sgdd", "sgddb", "sgdd2")],
    climate[climate$sgdd == max(climate$sgdd), c("sgdd", "sgddb", "sgdd2")]),
    rbind(
        climate[climate$wai == min(climate$wai), c("wai", "waib", "wai2")],
        climate[climate$wai == max(climate$wai), c("wai", "waib", "wai2")]
    ), N = 1:2)
tmp

tmp2 <- cbind(rbind(
    climate[climate$sgdd == min(climate$sgdd), 1:6],
    climate[climate$sgdd == max(climate$sgdd), 1:6],
        climate[climate$wai == min(climate$wai), 1:6],
        climate[climate$wai == max(climate$wai), 1:6]
    ),
    N = 1:4)
tmp2

tmp3 <- cbind(rbind(
    climate[climate$sgdd == min(climate$sgdd), c("sgdd", "sgddb", "sgdd2")],
    climate[climate$sgdd == max(climate$sgdd), c("sgdd", "sgddb", "sgdd2")]),
    rbind(
        climate[climate$wai == max(climate$wai), c("wai", "waib", "wai2")],
        climate[climate$wai == min(climate$wai), c("wai", "waib", "wai2")]
    ), N = 1:2)
tmp3

plot(climate$sgdd, climate$wai, cex = climate$N/ 100)
points(tmp$sgdd, tmp$wai, col = "red")
points(tmp2$sgdd, tmp2$wai, col = "chartreuse")
points(tmp3$sgdd, tmp3$wai, col = "turquoise")
z <- getRangemu_dev(climate = tmp, fit_Abies_alba, BA = 0:200, mesh = species.in$IPM$mesh) # red points
# min        max        sig
# -7.8596011  2.3885288  0.5967262
z2 <- getRangemu_dev(climate = tmp2, fit_Abies_alba, BA = 0:200, mesh = species.in$IPM$mesh) # green
# min        max        sig
# -8.2948181  2.2820031  0.5967262
z4 <- getRangemu_dev(climate = tmp3, fit_Abies_alba, BA = 0:200, mesh = species.in$IPM$mesh) # blue
# min        max        sig
# -8.3224903  2.2870663  0.5967262
z4 <- getRangemu_dev(climate, fit_Abies_alba, BA = 0:200, mesh = species.in$IPM$mesh)
# min        max        sig
# -8.2948181  2.3177412  0.5967262
