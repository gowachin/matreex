#' Testing the effect of harvest parameters on final distribution ####
#' Script for simulations with different settings on dth and dha combinations.

## Library ####
# library(matreex)
devtools::load_all() # Using
library(dplyr)
library(ggplot2)


## IPM integration ####

# Load fitted model for a species
# fit_species # list of all species in dataset
data("fit_Picea_abies")

# Load associated climate
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)

Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    correction = "none",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:60, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)

## Equilibrium run ####

Picea_sp <- species(IPM = Picea_ipm, init_pop = def_initBA(30))

Picea_for <- forest(species = list(Picea = Picea_sp))
set.seed(42) # The seed is here for initial population random functions.
Picea_sim <- sim_deter_forest(
    Picea_for,
    tlim = 2000,
    equil_time = 3000, equil_dist = 250, equil_diff = 1,
    SurfEch = 0.03,
    verbose = TRUE
)

distrib_equil <- Picea_sim %>%
    dplyr::filter(grepl("^n$", var), equil) %>%
    dplyr::pull(value)
distrib_equil <- distrib_equil * 0.03

Picea_tmp <- Picea_sp
Picea_tmp$init_pop <- def_init_k(distrib_equil)
Picea_tmp$harvest_fun <- Uneven_harv

simul <- function(sp= Picea_tmp, dth = 175, dha = 575, hmax = 1){

    sp$harv_lim <- c(dth = dth, dha = dha, hmax = hmax)
    Picea_for_tmp <- forest(species = list(Picea = sp),
                            harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                           freq = 5, alpha = 1))
    res <- sim_deter_forest(
        Picea_for_tmp,
        tlim = 500,
        equil_time = 500, equil_dist = 10, equil_diff = 1,
        harvest = "Uneven", targetBA = 20, # We change the harvest and set targetBA.
        SurfEch = 0.03,
        verbose = TRUE
    )

    return(res)
}

# Modify harv_lim
Picea_sim_f20_10_57 <- simul(Picea_tmp, dth = 100, dha = 575, hmax = 1)
Picea_sim_f20_10_10 <- simul(Picea_tmp, dth = 100, dha = 100, hmax = 1)
Picea_sim_f20_30_57 <- simul(Picea_tmp, dth = 300, dha = 575, hmax = 1)
Picea_sim_f20_30_30 <- simul(Picea_tmp, dth = 300, dha = 300, hmax = 1)

sim <- bind_rows(th10_ha57 = Picea_sim_f20_10_57,
                 th30_ha57 = Picea_sim_f20_30_57,
                 th30_ha30 = Picea_sim_f20_30_30,
                 th10_ha10 = Picea_sim_f20_10_10,
                 .id = "simulation") %>%
    mutate(lastcut = time == max(time[value[var== "H"] > 0 & var== "H"]))

var <- sim %>%
    filter(var %in% c("BAsp", "BAstand", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = simulation)) +
    geom_line() +
    facet_wrap(. ~ var, scales = "free_y") +
    NULL

distrib <- sim %>%
    filter(var %in% c("n", "h"), lastcut, size > 0, value > 0) %>%
    # mutate(value = 1- value / max(value)) %>%
    mutate(var = factor(var, levels = c("n", "h"))) %>%
    ggplot(aes(x = size, y = value, color = simulation,
               linetype = var)) +
    geom_line() +
    geom_vline(xintercept = c(100, 300, 575), linetype = "dashed", alpha = .5) +
    NULL

ggarrange(
    distrib,  var, nrow = 2, common.legend = TRUE,
    labels = list("Last cut time size distribution", "Time evolution")
)

# Modify hmax
Picea_sim_f20_10_57 <- simul(Picea_tmp, dth = 100, dha = 575, hmax = .8)
Picea_sim_f20_10_10 <- simul(Picea_tmp, dth = 100, dha = 100, hmax = .8)
Picea_sim_f20_30_57 <- simul(Picea_tmp, dth = 300, dha = 575, hmax = .8)
Picea_sim_f20_30_30 <- simul(Picea_tmp, dth = 300, dha = 300, hmax = .8)

sim <- bind_rows(th10_ha57 = Picea_sim_f20_10_57,
                 th30_ha57 = Picea_sim_f20_30_57 ,
                 th30_ha30 = Picea_sim_f20_30_30,
                 th10_ha10 = Picea_sim_f20_10_10, .id = "simulation") %>%
    mutate(lastcut = time == max(time[value[var== "H"] > 0 & var== "H"]))


var <- sim %>%
    filter(var %in% c("BAsp", "BAstand", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = simulation)) +
    geom_line() +
    facet_wrap(. ~ var, scales = "free_y") +
    NULL

distrib <- sim %>%
    filter(var %in% c("n", "h"), lastcut, size > 0, value > 0) %>%
    # mutate(value = 1- value / max(value)) %>%
    mutate(var = factor(var, levels = c("n", "h"))) %>%
    ggplot(aes(x = size, y = value, color = simulation, linetype = var)) +
    geom_line() +
    geom_vline(xintercept = c(100, 300, 575), linetype = "dashed", alpha = .5) +
    NULL

ggarrange(
    distrib,  var, nrow = 2, common.legend = TRUE,
    labels = list("Last cut time size distribution", "Time evolution")
)

