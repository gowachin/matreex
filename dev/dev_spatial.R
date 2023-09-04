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
load("dev/dev_mosa.RData")
## Mosaic forest ####

#' The most important thing is to prevent RAM overload


forests <- list(A = AbFa_for,
                B = invasive_AbFa)


# each forest require the migration rate I guess
mosaic <- function(forests = list()){

    sp <- map(forests, ~ .x$info$species) %>% flatten_chr() %>% unique()
    sp_ipms <- vector("list", length(sp))
    names(sp_ipms) <- sp

    # TODO write a test where all species are the same in all forests !

    for(i in seq_along(forests)){
        spi <- forests[[i]]$info$species
        tmp_sub <- names(sp_ipms[spi])

        for(sp_i in tmp_sub){
            if(is.null(sp_ipms[[sp_i]])){
                sp_ipms[[sp_i]] <- forests[[i]]$species[[sp_i]]$IPM
            }

            forests[[i]]$species[[sp_i]]$IPM$IPM <- NULL
            forests[[i]]$species[[sp_i]]$IPM$fit$fec <- forests[[i]]$species[[sp_i]]$IPM$fit$rec
            forests[[i]]$species[[sp_i]]$IPM$fit$fec$params_m[c("BATOTSP", "BATOTNonSP")] <- 0
        }
    }

    mosaic <- list(ipms = sp_ipms,
                   forests = forests)

    class(mosaic) <- "mosaic"

    return(mosaic)

}

Mosaic <- mosaic(forests)

names(Mosaic$forests)

source("dev/dev_env.R")

## Simulations ####
# sim_deter_mosaic <- function(Mosaic,
#                              tlim = 3e3,
#                              harvest = c("default", "Uneven", "Even"),
#                              targetBA = 20,
#                              targetRDI = 0.9,
#                              targetKg = 0.9,
#                              final_harv = 100,
#                              climate = NULL,
#                              # require a disturb table that indicates which plot to disturb
#                              disturbance = NULL,
#                              correction = "none",
#                              SurfEch = 0.03,
#                              verbose = FALSE){


    # browser()
    tlim = 100
    harvest = "default"
    climate = NULL
    disturbance = NULL
    correction = "none"
    SurfEch = 0.03
    verbose = TRUE

    # Idiot Proof ####
    # validate_forest(Forest) # TEMP dev
    assertCount(tlim)
    harvest <- match.arg(harvest, c("default", "Uneven", "Even"))
    # assertNumber(targetBA, lower = 0)
    # assertNumber(targetRDI, lower = 0, upper = 1) # FIXME single or species target ?
    # assertNumber(targetKg, lower = 0, upper = 1)
    IPM_cl <- map_chr(Mosaic$ipms, class)
    if(all(IPM_cl == "ipm") && !is.null(climate)) {
        # no climate needed
        warning(paste0("Because all species are fully integrated on a climate, ",
                       "providing one now is unnecessary"))
        clim_i <- which(IPM_cl == "ipm")[[1]]
        climate <- t(Mosaic$ipms[[clim_i]]$climatic)
        climate <-  as.matrix(bind_cols(climate, t = 1:tlim))
    } else
        {
        if(any(IPM_cl == "ipm")){
            if(!is.null(climate)){
                warning(
                    paste0("At least one species is fully integrated on a ",
                           "climate, so this climate will be used for simulation"))
            }
            clim_i <- which(IPM_cl == "ipm")[[1]]
            climate <- t(Mosaic$ipms[[clim_i]]$climatic)
        }
        if(inherits(climate, "data.frame")){
            climate <- as.matrix(climate)
        } else if(inherits(climate, "numeric")){
            climate <- t(climate)
        }
        assertMatrix(climate)
        if(nrow(climate) != 1 & nrow(climate) < tlim){
            stop(paste0("climate matrix is not defined for each time until",
                        " tlim. This matrix require a row per time or ",
                        "single one."))
        }
        if(nrow(climate) == 1){
            climate <-  as.matrix(bind_cols(climate, t = 1:tlim))
        }
        }


    run_disturb <- !is.null(disturbance)
    if(run_disturb){
        # TODO idiot proof disturbance
        t_disturb <- logical(equil_time)
        t_disturb[disturbance$t] <- TRUE
        if(any(disturbance$intensity <= 0 | disturbance$intensity > 1)){
            warning("Disturbances with intensity outside ]0;1] have been removed")
            disturbance <- disturbance[disturbance$intensity > 0 &
                                           disturbance$intensity < 1,]
            if(nrow(disturbance) == 0){
                warning("There is no disturbances left with correct intensity.")
                disturbance <- NULL
                run_disturb <- FALSE
            }
        }
    }

    correction <- match.arg(correction, c("cut", "none"))
    assertNumber(SurfEch, lower = 0)
    assertLogical(verbose, any.missing = FALSE, len = 1)

    start <- Sys.time()

    # Initialisation ####
    init_sim <- function(nsp, tlim, mesh){ # TODO : set function outside of here
        res <- vector("list", nsp)
        res <- map2(lengths(mesh), names(mesh), function(x, y) {
            tmp <- matrix(
                data = NA_real_, ncol = tlim ,
                nrow = x + 3 + x + 1
            )
            colnames(tmp) <- c(paste0("t", 1:tlim)) #, "sp")
            rownames(tmp) <- c(paste0(y, ".n", 1:x),
                               paste0(y, c(".BAsp", ".BAstand", ".N")),
                               paste0(y, ".h", 1:x), paste0(y,".H"))
            return(tmp)
        })
        names(res) <- names(mesh)
        return(res)
    }
    nplot <- length(Mosaic$forests)
    nms_plot <- names(Mosaic$forests)
    disturb_surv <- TRUE

    ## Modify IPM ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
    }
    # correct also decompress integer to double with x * 1e-7 app
    Mosaic$ipms <- map(Mosaic$ipms, correction.ipm, correction = correction)

    ## Initiate variables and populations ####
    landscape <- map(nms_plot,
                     ~ init_forest_env(Mosaic, index = .x,
                                       tlim = tlim, SurfEch = SurfEch)
    )
    # save first pop
    landscape <- map(landscape, ~ save_step_env(.x, t = 1))

    # Create sim IPM ####
    sim_clim <- climate[1, , drop = TRUE] # why this line ?
    landscape <- map(landscape, ~ get_step_env(.x, Mosaic, t= 1,
                              climate, correction))
    if (verbose) {
        message("Starting while loop. Maximum t = ", tlim)
    }

    # While tlim & eq ####
    t <- 2
    the <- NA_real_ # real time of simulation ending in case of outbound BA

    while (t < tlim ) {

        sim_clim <- climate[t, , drop = TRUE]
        # Growth, Disturbance, Harvesting
        landscape <- map(landscape, ~ growth_mortal_env(
            .x, t = t, harvest = harvest, run_disturb
            ))
        ### Recruitment ####
        landscape <- recrut_env(landscape, t)
        ## Save BA ####
        landscape <- map(landscape, ~ save_step_env(.x, t = t))
        ## Get sim IPM ####
        landscape <- get_step_env(landscape, Mosaic, t= t, climate, correction)
        ## Loop Verbose ####
        if (t %% 500 == 0 && verbose) {
            message(sprintf(
                "time %i",
                t#, diff(range(sim_BA[max(1, t - equil_dist):t]))
            ))
        }
        t <- t + 1
    }

#     # Format output ####
#     landsim <- map(landscape, function(plot){
#
#         tmp <- new_deter_sim(plot$sim_X, mesh = plot$meshs)
#         return(tree_format(tmp))
#     })
#     names(landsim) <- nms_plot
#
#     if (verbose) {
#         message("Simulation ended after time ", ifelse(is.na(the), t-1, the))
#         tmp <- Sys.time() - start
#         message("Time difference of ", format(unclass(tmp), digits = 3),
#                 " ", attr(tmp, "units"))
#     }
#
#     # Return ####
#     final <- dplyr::bind_rows(landsim, .id = "plot")
#
#     return(final)
# }