# Testing modification of

library(ggplot2)
devtools::load_all()
species <- "Pinus_pinea"
# You can choose any species in the dataset.
data("fit_species") # available species

# Data sets ####
# load data with fit models and place it in fit object.
data(list = paste0("fit_", species))

# load data et select optimum climate (N == 2, see ?climate_species)
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == species, select = -sp)
climate <- drop(as.matrix(climate)) # we need it as a vector.



dev_IPM <- dev_make_IPM(
    species = species, climate = climate, fit =  fit_Pinus_pinea,
    clim_lab = "optimum clim",
    mesh = c(m = 6, L = 100,  U = 1300), diag_tresh = 500, level = c(3, 4),
    BA = 0:30, verbose = TRUE
)

IPM <- make_IPM(
    species = species, climate = climate, fit =  fit_Pinus_pinea,
    clim_lab = "optimum clim",
    mesh = c(m = 6, L = 100,  U = 1300), diag_tresh = 500, level = 12,
    BA = 0:30, verbose = TRUE
)

# The algorithm change nothing so everything is fine.

dev_IPM <- dev_make_IPM(
    species = species, climate = climate, fit =  fit_Pinus_pinea,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90,
             U = as.numeric(fit_Pinus_pinea$info[["max_dbh"]]) * 1.1),
    level = c(1, 140),
    BA = 0:30, verbose = TRUE
)

IPM <- make_IPM(
    species = species, climate = climate, fit =  fit_Pinus_pinea,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90,
             U = as.numeric(fit_Pinus_pinea$info[["max_dbh"]]) * 1.1),
    level = 420,
    BA = 0:30, verbose = TRUE
)

microbenchmark::microbenchmark(
    dev1_IPM = dev_make_IPM(
        species = species, climate = climate, fit =  fit_Pinus_pinea,
        clim_lab = "optimum clim",
        mesh = c(m = 700, L = 90,
                 U = as.numeric(fit_Pinus_pinea$info[["max_dbh"]]) * 1.1),
        level = c(1, 140),
        BA = 0:200, verbose = TRUE
    ),
    dev_IPM = dev_make_IPM(
        species = species, climate = climate, fit =  fit_Pinus_pinea,
        clim_lab = "optimum clim",
        mesh = c(m = 700, L = 90,
                 U = as.numeric(fit_Pinus_pinea$info[["max_dbh"]]) * 1.1),
        level = c(3, 140),
        BA = 0:200, verbose = TRUE
    ),
    IPM = make_IPM(
        species = species, climate = climate, fit =  fit_Pinus_pinea,
        clim_lab = "optimum clim",
        mesh = c(m = 700, L = 90,
                 U = as.numeric(fit_Pinus_pinea$info[["max_dbh"]]) * 1.1),
        level = 420,
        BA = 0:200, verbose = TRUE
    ),
    times = 1
)

dev_IPM$int_log
IPM$int_log

n_ipm <- dev_IPM$IPM
o_ipm <- IPM$IPM

res <- purrr::map2(n_ipm, o_ipm, ~ .x - .y)
res <- purrr::map(res, ~ reshape2::melt(as.matrix(.x)))
res <- dplyr::bind_rows(res, .id = "BA")
res <- dplyr::mutate(res, BA = as.numeric(BA))

res %>%
    dplyr::filter(BA < 5) %>%
    ggplot(aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(na.value="transparent") +
    facet_wrap(~ BA, scales = "free_y") +
    theme_dark() +
    coord_cartesian(
        xlim = c(0, 50), ylim = c(0, 50)) +
    NULL

d <- 0
res %>%
    dplyr::mutate(diag = Var1 == Var2 +d, subdiag = Var1 == Var2 +d+1,
                  color = ifelse(diag, "diag", "subdiag")) %>%
    dplyr::filter(BA < 5,  diag | subdiag) %>%
    ggplot(aes(x = Var1, y = value, color = color)) +
    geom_line() +
    scale_fill_viridis_c(na.value="transparent") +
    facet_wrap(~ BA, scales = "free_y") +
    theme_dark() +
    NULL


old_time <- 1000

dev_species <- matreex::new_species(
    delay(dev_IPM, delay = 10), harvest_fun = def_harv, init_pop = def_init)
set.seed(42)
dev_sim <- sim_deter_forest(
        dev_species, tlim = old_time, equil_time = old_time,
        equil_dist = 2, verbose = TRUE, correction = "cut") %>%
    tree_format() %>% dplyr::filter(var == "BAsp")

o_species <- matreex::new_species(
    delay(IPM, delay = 10), harvest_fun = def_harv, init_pop = def_init)
set.seed(42)
sim <- sim_deter_forest(
    o_species, tlim = old_time, equil_time = old_time,
    equil_dist = 2, verbose = TRUE, correction = "cut") %>%
    tree_format() %>% dplyr::filter(var == "BAsp")

memor <- dplyr::bind_rows(list(int1 = dev_sim, int3 = sim), .id = "version")

memor %>%
    dplyr::filter(! equil, value != 0) %>%
    ggplot(aes(x = time, y = value, color = version)) +
    # facet_wrap(~ species, scales = "free_y") +
    geom_line(size = .4) +
    # coord_cartesian(xlim = c(400, 550), ylim = c(23, 25)) +
    NULL
