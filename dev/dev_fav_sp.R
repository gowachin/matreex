
sp_name <- c("toto", "tutu")

lapply(sp_name, isTRUE)
sapply(sp_name, isTRUE, USE.NAMES = TRUE)
library(purrr)

library(tidyr)
# library(matreex)
devtools::load_all()
library(dplyr)
library(ggplot2)
options(W_matreex_edist = FALSE)
# Load fitted model for a species
# fit_species # list of all species in dataset
data("fit_Picea_abies")
data("fit_Fagus_sylvatica")

# Load associated climate
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
# see ?climate_species to understand the filtering of N.
climate

Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:70, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
Fagus_ipm <- make_IPM(
    species = "Fagus_sylvatica",
    climate = climate,
    fit = fit_Fagus_sylvatica,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
    BA = 0:70, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)

# <!--
# #### Favoured species
#
# In some case, we may want to favour some species compared to others.
# We note $Q$ the $q$ species we want to favour, and $P_Q=\sum_{i=1}^{q} p_i$
# We first compute the harvest rate $H_Q$ and $H_{1-Q}$ for respectively the favoured/other species.
#
# If $P_Q \geq 0.5$, we take $H_Q = H_{1-Q} = H$ (the species to be favoured are already dominant).
#
#
# If $P_Q < 0.5$, we compute $H_Q=f(P_Q) H$ and $H_{1-Q}=f(1-P_Q) H$ with $\alpha > 1$.
# By definition, we will get $H_Q \leq H$.
#
# We then apply for each species $i$ the harvest rate $H_i = H_Q$ or $H_i=H_{1-Q}$ depending on which group it belongs to.
# -->
devtools::load_all()
Picea_Uneven <- species(IPM = Picea_ipm, init_pop = def_initBA(30),
                        harvest_fun = Uneven_harv,
                        harv_lim = c(dth = 175, dha = 575, hmax = 1))
Fagus_Uneven <- species(IPM = Fagus_ipm, init_pop = def_initBA(30),
                        harvest_fun = Uneven_harv,
                        harv_lim = c(dth = 175, dha = 575, hmax = 1))

cases <- list(nofav = c(Picea_abies = FALSE, Fagus_sylvatica = FALSE),
             fav_fagus = c(Picea_abies = FALSE, Fagus_sylvatica = TRUE),
             fav_picea = c(Picea_abies = TRUE, Fagus_sylvatica = FALSE),
             fav_all = c(Picea_abies = TRUE, Fagus_sylvatica = TRUE)
)

sim_ls <- purrr::map(cases, function(vec) {
    PiFa_for_Uneven <- forest(species = list(Picea = Picea_Uneven,
                                             Fagus = Fagus_Uneven),
                              harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                             freq = 5, alpha = 2),
                              favoured_sp = vec
    )

    set.seed(42) # The seed is here for initial population random functions.
    Picea_sim_f20 <- sim_deter_forest(
        PiFa_for_Uneven,
        tlim = 260,
        equil_time = 260, equil_dist = 10, equil_diff = 1,
        harvest = "Favoured_Uneven", targetBA = 20, # We change the harvest and set targetBA.
        # harvest = "Uneven", targetBA = 20,
        SurfEch = 0.03,
        verbose = TRUE
    )

    return(Picea_sim_f20)
})


devtools::load_all()
vec <- c(Picea_abies = FALSE, Fagus_sylvatica = TRUE)
 PiFa_for_Uneven <- forest(species = list(Picea = Picea_Uneven,
                                                   Fagus = Fagus_Uneven),
                                    harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                                   freq = 5, alpha = 2),
                                    favoured_sp = vec
)

set.seed(42) # The seed is here for initial population random functions.
Picea_sim_f20 <- sim_deter_forest(
    PiFa_for_Uneven,
    tlim = 460,
    equil_time = 460, equil_dist = 10, equil_diff = 1,
    harvest = "Favoured_Uneven", targetBA = 20, # We change the harvest and set targetBA.
    # harvest = "Uneven", targetBA = 20,
    SurfEch = 0.03,
    verbose = TRUE
)

sim_ls <- list(fav_all = Picea_sim_f20)

sim_ls %>% dplyr::bind_rows(.id = "sims") %>%
    dplyr::filter(var %in% c("BAstand", "H"), ! equil, value > 0) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    facet_grid(sims ~ var, scales = "free_y") +
    geom_line(linewidth = .2) + geom_point(size = 0.4)
