load_all()
document()
data("fit_Picea_abies")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
climate
Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
Picea_sp <- species(IPM = Picea_ipm, init_pop = def_initBA(30))
Forest <- forest(species = list(Picea = Picea_sp))
sim1 <- sim_indiv_forest(Forest, verbose = TRUE, tlim = 5000)

library(ggplot2)
sim1  %>%
    ggplot(aes(x = time, y = value)) +
    geom_line(linewidth = .4) + ylab("BA") +
    facet_wrap(~ var, scales = "free_y")

