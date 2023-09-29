# Mosaic idea ####

#' Idea is to simulate multiple forest in parallel and put recrutment in common
#' or from a regional pool of species to a forest.
#'
#' I think the most difficult aspect of this is to define a correct initiation state
#' for the plots and/or the regional pool.
#'
#' I should not duplicate the ipm when creating the object class


#' Notes from Georges
#' For the mosaic with  N patches the recruitment rate of species $s$ is given
#' by $1/N \times \sum_{k=1}^N \alpha_s *BA_{sk}^{\beta_s}$


if(FALSE){
## Common values to test on it ####

# library(matreex)
devtools::load_all()
library(dplyr)
library(ggplot2)
data("fit_Abies_alba")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Abies_alba", select = -sp)
climate

Abies_ipm <- make_IPM(
    species = "Abies_alba",
    climate = climate,
    fit = fit_Abies_alba,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Abies_alba) * 1.1),
    BA = 0:120,
    verbose = TRUE
)

Abies_ipm$fit$rec

options(W_matreex_edist = FALSE)
Abies_sp <- species(IPM = Abies_ipm, init_pop = def_initBA(30))

Abies_sp$recruit_fun

Abies_for <- forest(species = list(Abies = Abies_sp))
set.seed(42) # The seed is here for initial population random functions.
Abies_sim <- sim_deter_forest(
    Abies_for,
    tlim = 2000,
    equil_time = 2000, equil_dist = 1, equil_diff = 1,
    SurfEch = 0.03,
    verbose = TRUE
)

data("fit_Fagus_sylvatica")
Fagus_ipm <- make_IPM(
    species = "Fagus_sylvatica",
    climate = climate,
    fit = fit_Fagus_sylvatica,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
    BA = 0:120, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
Fagus_sp <- species(IPM = Fagus_ipm, init_pop = def_initBA(30))
AbFa_for <- forest(species = list(Abies = Abies_sp,  Fagus = Fagus_sp))
set.seed(42) # The seed is here for initial population random functions.
AbFa_sim <- sim_deter_forest(
    AbFa_for,
    tlim = 2000,
    equil_time = 2000, equil_dist = 1, equil_diff = 1,
    SurfEch = 0.03,
    verbose = TRUE
)

## Regional pool ####

### regional forest class ####
#' species are initiated with their own functions and only the regional abundance
#' is set in the forest object as well as the migration rate

#' Constructor for forest class
#'
#' Only used in the matreex package
#'
#' @param species List of species created with matreex package.
#' @param harv_rules Vector for harvest rules at the scale of the forest.
#' \describe{
#'   \item{Pmax}{maximum proportion of BAcut / BA}
#'   \item{dBAmin}{the minimum BA to perform cut}
#'   \item{freq}{Frequence at which the harvest will be executed.}
#' }
#' @param regional_abundance list of vector with size distribution for each species.
#' This is not a direct species relative abundance, but I don't know how to implement this...help ,
#' @param migration_rate numeric vector with a migration rate in percentage between 1 et and 0.
#'
#' @importFrom purrr map_chr
#'
#' @keywords internal
#' @export
new_forest <- function(species = list(),
                       harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                       regional_abundance = NULL,
                       migration_rate = NULL
){

    sp <- map_chr(species, sp_name)
    names(species) <- sp
    if(!is_null(regional_abundance)){
        names(regional_abundance) <- sp
        names(migration_rate) <- sp
    }
    forest <- list(
        species = species, harv_rules = harv_rules,
        info = list(species = sp,
                    clim_lab = map_chr(species, climatic)),
        regional_abundance = regional_abundance,
        migration_rate = migration_rate
    )

    if(!is_null(regional_abundance)){
        class(forest) <- c("forest", "reg_forest")
    } else {
        class(forest) <- "forest"
    }

    return(forest)
}

#' validator for forest class.
#'
#' @param x forest class object
#'
#' @import checkmate
#'
#' @noRd
validate_forest <- function(x){

    regional <- inherits(x, "reg_forest")
    values <- unclass(x)
    names <- attr(x, "names")

    map(values$species, validate_species)
    # TODO check forest harv rules

    clim_lab <- values$info$clim_lab
    if(length(unique(clim_lab)) > 1){
        clim_ipm <- clim_lab[clim_lab != "mu_gr"]
        if(length(clim_ipm) > 1){ # D & F
            stop(paste0("Some ipm species are not defined with the same climatic name.",
                        "Check it with : map_chr(species, climatic)"))
        }
    }

    # check the regional pool settings
    if(regional){

        assertNumeric(values$migration_rate, len = length(values$species),
                      lower = 0, upper = 1)
        if(all(values$migration_rate == 0)){
            warning("All migration rate are 0, the regional pool of this forest is deleted")
            x$regional_abundance <- NULL
            x$migration_rate <- NULL
            class(x) <- "forest"

            return(invisible(x))
        }

        assertSubset(names(values$migration_rate), names(values$species))
        # length_meshs <- map_dbl(values$species, ~ length(.x$IPM$mesh))

        # assertList(values$regional_abundance, types = "numeric",
                   # len = length(values$species))
        # if(any(lengths(values$regional_abundance) != length_meshs)){
            # stop("regional abundance numeric vector should be the length of the species mesh.")
        # }

        assertSubset(names(values$regional_abundance), names(values$species))
    }

    invisible(x)
}

#' Create a new forest for simulation
#'
#' A forest is a group of one of multiple species to silumate along time using
#' the IPM defined for each species and harvest rules.
#'
#' @inheritParams new_forest
#'
#' @export
forest <- function(species = list(),
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                   regional_abundance = NULL,
                   migration_rate = NULL
){

    res <- validate_forest(new_forest(
        species = species,
        harv_rules = harv_rules,
        regional_abundance = regional_abundance,
        migration_rate = migration_rate
    ))

    return(res)
}

### Trying simulation ####

equil_dist <- dplyr::filter(Abies_sim, equil, var == "n") %>% dplyr::pull(value)
equil_BA <- dplyr::filter(Abies_sim, equil, var == "BAsp") %>% dplyr::pull(value)

Abies_spE <- Abies_sp
Abies_spE$init_pop <- def_init_k(equil_dist)

regional_Abies <- forest(
    species = list(Abies = Abies_spE),
    regional_abundance = equil_BA,
    migration_rate = c(Abies = 0.05)
)


e2_dist <- dplyr::filter(AbFa_sim, equil, var == "n") %>%
    group_by(species) %>%
    dplyr::group_split() %>% map(pull, value)
e2_BA <- dplyr::filter(AbFa_sim, equil, var == "BAsp") %>% dplyr::pull(value)


# Testing invasive species
Fagus_spE0 <- Fagus_sp
Fagus_spE0$init_pop <- def_init_k(e2_dist[[2]] * 0)


invasive_AbFa <- forest(
    species = list(Abies = Abies_spE, Fagus = Fagus_spE0),
    regional_abundance = e2_BA,
    migration_rate = c(Abies = 0.1, Fagus = 0.1)
)


rm(AbFa_sim, Abies_ipm, Abies_sim, Abies_sp, Abies_spE, climate_species,
   e2_BA, e2_dist, Fagus_sp, Fagus_spE0, fit_Abies_alba, fit_Fagus_sylvatica,
   equil_BA, equil_dist)


save.image("dev/dev_mosa.RData")

}


library(dplyr)
library(purrr)
library(ggplot2)
devtools::load_all()
load("dev/dev_mosa.RData")
source("dev/dev_env.R")
## Mosaic forest ####
#' The most important thing is to prevent RAM overload
forests <- list(A = AbFa_for,
                B = invasive_AbFa)
Mosaic <- mosaic(forests)
names(Mosaic$forests)

disturb <- list(A = data.frame(type = "storm", intensity = 0.2,
                                          IsSurv = FALSE, t = 100),
                B = data.frame(type = "storm", intensity = 0.2,
                               IsSurv = FALSE, t = 900))


# Simulations ####
# profvis::profvis({
    res <- sim_deter_mosaic(Mosaic, tlim = 1000, harvest = "default",
                            disturbance = NULL,
                            verbose = TRUE)
# })


disturb_fun <- function(x, species, disturb = NULL, ...){

    dots <- list(...)
    qmd <- dots$qmd
    size <- species$IPM$mesh
    coef <- species$disturb_coef
    if(any(disturb$type %in% coef$disturbance)){
        coef <- subset(coef, disturbance == disturb$type)
    } else {
        stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                     sp_name(species), disturb$type))
    }

    # edits for delay
    size[size == 0] <- min(size[size !=0])

    logratio <-  log(size / qmd)
    dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
    logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
    Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                        coef$b * disturb$intensity ^(coef$c * dbh.scaled))

    return(x* Pkill) # always return the mortality distribution
}

Mosaic$forests$A$species$Abies_alba$disturb_fun <- disturb_fun
Mosaic$forests$A$species$Fagus_sylvatica$disturb_fun <- disturb_fun

Mosaic$forests$A$species$Abies_alba$disturb_coef <- filter(
    matreex::disturb_coef, species == "Abies_alba")

Mosaic$forests$A$species$Fagus_sylvatica$disturb_coef <- filter(
    matreex::disturb_coef, species == "Fagus_sylvatica")

res <- sim_deter_mosaic(Mosaic, tlim = 1000, harvest = "default",
                        disturbance = disturb,
                        verbose = TRUE)


res %>%
    filter(var %in% c("N", "BAsp")) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line() +
    facet_wrap(plot ~ var, scales = "free_y") +
    NULL

filter(res, is.na(value))
res$value


# testing timing ####
sim_deter_forest(AbFa_for, tlim = 1000, harvest = "default", verbose = TRUE)
# Starting while loop. Maximum t = 10000
# time 500 | BA diff : 5.01
# time 1000 | BA diff : 1.16
# Simulation ended after time 1160
# BA stabilized at 56.42 with diff of 1.00 at time 1160
# Time difference of 14 secs
res <- sim_deter_mosaic(mosaic(list(A = AbFa_for)),
                        tlim = 1000, harvest = "default", verbose = TRUE)
# Starting while loop. Maximum t = 1000
# time 500
# Simulation ended after time 999
# Time difference of 14 secs
res <- sim_deter_mosaic(mosaic(list(A = AbFa_for, B = AbFa_for)),
                        tlim = 1000, harvest = "default", verbose = TRUE)
# Starting while loop. Maximum t = 1000
# time 500
# Simulation ended after time 999
# Time difference of 28.2 secs
res <- sim_deter_mosaic(mosaic(list(A = AbFa_for, B = AbFa_for, C = AbFa_for)),
                        tlim = 1000, harvest = "default", verbose = TRUE)
# Starting while loop. Maximum t = 1000
# time 500
# Simulation ended after time 999
# Time difference of 42.5 secs
