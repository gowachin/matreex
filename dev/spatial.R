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
#'
#' For the immigration from a regional pool with constant species abundance
#' $R_s$ the recruitment would be
#' $\alpha_s *BA_{sk}^{\beta_s} + m \times R_s$


## Common values to test on it ####

# library(matreex)
devtools::load_all()
data("fit_Picea_abies")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
# N here is a climate defined in Kunstler et al 2021.
# N == 2 is the optimum climate for the species.
# see ?climate_species for more info.
climate

Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:60, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
options(W_matreex_edist = FALSE)
Picea_sp <- species(IPM = Picea_ipm, init_pop = def_initBA(30))
Picea_for <- forest(species = list(Picea = Picea_sp))
set.seed(42) # The seed is here for initial population random functions.
Picea_sim <- sim_deter_forest(
    Picea_for,
    tlim = 200,
    equil_time = 300, equil_dist = 50, equil_diff = 1,
    SurfEch = 0.03,
    verbose = TRUE
)

## Regional pool ####

### regional forest class ####
#' species are initiated with their own functions and only the regional abundance

regional_forest <- list(
    species = list(),
    regional_abondance = list(), # list of abundance. names are species and this should have the length of the
    migration_rate = c(), # species named migration rate ?
    harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), # Keeping it ?
    info = list(species = sp, clim_lab = map_chr(species, climatic))
)

class(regional_forest) <- c("forest", "reg_forest")



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
#' @importFrom purrr map_chr
#'
#' @keywords internal
#' @export
new_forest <- function(species = list(),
                       harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                      freq = 1, alpha = 1)
                       regional_abondance = NULL,
                       migration_rate = NULL
){

    sp <- names(species) <- map_chr(species, sp_name)
    forest <- list(
        species = species, harv_rules = harv_rules,
        info = list(species = sp,
                    clim_lab = map_chr(species, climatic))
        regional_abondance = regional_abondance,
        migration_rate = migration_rate
    )

    if(!is_null(regional_abondance)){
        class(forest) <- c("forest", "reg_forest")
    } else {
        class(forest) <- "forest"
    }

    return(forest)
}

#' validator for IPM class.
#'
#' @param x IPM class object
#'
#' @import checkmate
#'
#' @noRd
validate_forest <- function(x){

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
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1)
){

    res <- validate_forest(new_forest(
        species = species,
        harv_rules = harv_rules
    ))

    return(res)
}




#' I will use sim_deter_forest.forest function and just edit it to accept the regional migration I guess ?

### Trying simulation ####

# the idea is to do this ^^
a <- 1
class(a) <- "cat"
b <- 2
class(b) <- c("cat", "anis")
foo <- function(x){UseMethod("foo")}
foo.cat <- function(x){
    print(x)

    if(inherits(x, "anis")){
        print("I love this cat")
    }
}
foo(a)
foo(b)

## Mosaic forest ####



## Simulations ####
sim_deter_mosaic <- function(x,
                             tlim = 3e3,
                             equil_dist = 250,
                             equil_diff = 1,
                             equil_time = 1e4,
                             harvest = c("default", "Uneven", "Even"),
                             targetBA = 20,
                             targetRDI = 0.9,
                             targetKg = 0.9,
                             final_harv = 100,
                             climate = NULL,
                             # require a disturb table that indicates which plot to disturb
                             disturbance = NULL,
                             correction = "none",
                             SurfEch = 0.03,
                             verbose = FALSE){
    res <- NULL
    return(res)
}
