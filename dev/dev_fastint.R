# Testing personnal function ####
library(ggplot2)
library(dplyr)
load_all()
rm(list = ls())
x <- make_mu_gr(species = "Picea_abies", fit_Picea_abies,
                 mesh = c(m = 700, L = 90,
                          U = as.numeric(fit_Picea_abies$info["max_dbh"]) * 1.1),
                 verbose = TRUE, stepMu = 0.001)
# fit_Picea_abies
lobstr::obj_size(x) # 2.87 MB

mu_Picea <- new_species(IPM = x, init_pop = def_initBA(20), harvest_fun = def_harv)
# mu_Picea$info["species"] <- "mPicea_abies"

climate <- subset(climate_species, sp == "Picea_abies" & N ==2, select = -c(N, sp))
climate <- drop(as.matrix(climate))
# BA <- 20 ; species <- mu_Picea ; mu_gr <- x ; verbose = TRUE
# load_all()
test <- get_step_IPM.mu_gr(x = x, BA = 20, climate = climate, sim_corr = "cut")

# IPM <- old_ipm2ipm(species = "Picea_abies", climatic = 2)
# max(abs((IPM$IPM[[20]] * 1e-7) - test))
# (IPM$IPM[[20]] * 1e-7)[1:5, 1:5]
# test[1:5, 1:5]


# working in sim deter
spsel <- "Picea_abies"
climate <- subset(climate_species, sp == spsel & N == 2, select = -sp)
climate <- drop(as.matrix(climate))

Picea_ipm <- make_IPM("Picea_abies", climate, "opt_clim", fit_Picea_abies,
                    mesh = c(m = 700, L = 90,
                             U = as.numeric(fit_Picea_abies$info["max_dbh"]) * 1.1),
                    BA = 0:200, verbose = TRUE
)

load_all()
microbenchmark::microbenchmark(
    o = get_step_IPM.mu_gr(x = x, BA = 20, climate = climate, sim_corr = "cut"),
    ipm = get_step_IPM.ipm(x = Picea_ipm, BA = 20, climate = climate, sim_corr = "cut")
)

forest <- new_forest(species = list(mu_Picea = mu_Picea))
time <- 1000
load_all()
set.seed(42)
profvis::profvis({
    sim_deter_forest.forest(forest, tlim = time, equil_dist = time, equil_time = time,
                     verbose = TRUE, correction = "none")
})


