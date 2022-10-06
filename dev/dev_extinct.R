# Test if sim work on BA = 0

load_all()
fit  <- old_fit2fit("Yggdrasil", mean = TRUE)
# replace above line with data("fit_yggdrasil")


clim <- readRDS("output/Yggdrasil/plots_pred.Rds")
climate <- unlist(clim[2,])

# Testing 0 IPM ####
ipm <- make_IPM(species = "Yggdrasil", climate = climate,
                clim_lab = "test", fit = fit,
                mesh = c(m = 700, L = 90, U = 1600), BA = 0:6,
                verbose = TRUE) # 2.4min

Yggdrasil <- species(IPM = ipm, init_pop = def_initBA(5),
                     harvest_fun = def_harv)

Yggdrasil$recruit_fun <- function(BATOTSP, BATOTNonSP, mesh, SurfEch = 0.03){
    res <- rep(0, length(mesh))
    return(res)
}

load_all()
set.seed(666)
x <- sim_deter_forest(Yggdrasil, tlim = 300, verbose = TRUE)
y <- tree_format(x)
head(y)
library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
y %>%
    mutate(across(where(is_character),as.factor)) %>%
    filter(var == "BAsp") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    NULL

y %>%
    mutate(across(where(is_character),as.factor)) %>%
    filter(var == "BAsp") %>% select(value) %>% range()


# Testing highet BA increment ####
ipm <- make_IPM(species = "Yggdrasil", climate = climate,
                clim_lab = "test", fit = fit,
                mesh = c(m = 700, L = 90, U = 1600), BA = seq(0, 80, by = 4),
                verbose = TRUE) # 2.4min

Yggdrasil <- species(IPM = ipm, init_pop = def_initBA(10),
                     harvest_fun = def_harv)

load_all()
set.seed(666)
x <- sim_deter_forest(Yggdrasil, tlim = 1000, verbose = TRUE)
y <- tree_format(x)
head(y)
library(ggplot2, quietly = TRUE)
library(dplyr, quietly = TRUE)
y %>%
    mutate(across(where(is_character),as.factor)) %>%
    filter(var == "BAsp") %>%
    ggplot(aes(x = time, y = value)) +
    geom_line() +
    NULL
