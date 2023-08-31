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


## Mosaic forest ####

#' The most important thing is to prevent RAM overload


forests <- list(A = regional_Abies,
                B = invasive_AbFa)


# each forest require the migration rate I guess
mosaic <- function(forests = list()){

    sp <- map(forests, ~ .x$info$species) %>% flatten_chr() %>% unique()
    sp_ipms <- vector("list", length(sp))
    names(sp_ipms) <- sp

    for(i in seq_along(forests)){
        spi <- forests[[i]]$info$species
        tmp_sub <- names(sp_ipms[spi])

        for(sp_i in tmp_sub){
            if(is.null(sp_ipms[[sp_i]])){
                sp_ipms[[sp_i]] <- forests[[i]]$species[[sp_i]]$IPM
            }

            forests[[i]]$species[[sp_i]]$IPM <- NULL
        }
    }

    mosaic <- list(ipms = sp_ipms,
                   forests = forests)

    class(mosaic) <- "mosaic"

    return(mosaic)

}

Mosaic <- mosaic(forests)

names(Mosaic$forests)

## Simulations ####
sim_deter_mosaic <- function(Mosaic,
                             tlim = 3e3,
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


    # browser()
    tlim = 500
    harvest = "default"
    climate = NULL
    disturbance = NULL
    correction = "none"
    SurfEch = 0.03
    verbose = TRUE

    # Idiot Proof ####
    # validate_forest(Forest) # TEMP dev
    assertCount(tlim)
    assertCount(equil_dist)
    assertNumber(equil_diff)
    assertNumber(equil_time)
    if (equil_time < tlim || equil_time < equil_dist) {
        stop("equil_time must be higher or equal to tlim and equil_dist")
    }

    if(getOption("W_matreex_edist") && equil_dist > 1 && equil_dist < 1000){
        warning(paste0(
            "The equil_dist value is low and could lead to inappropriate ",
            "equilibrium states. A recommended value is 1000. ",
            "\nThis warning  will be printed once by session and is ",
            "desactivable with options(W_matreex_edist = FALSE)"
        ))
        options(W_matreex_edist = FALSE)
    }


    harvest <- match.arg(harvest)
    # assertNumber(targetBA, lower = 0)
    # assertNumber(targetRDI, lower = 0, upper = 1) # FIXME single or species target ?
    # assertNumber(targetKg, lower = 0, upper = 1)
    IPM_cl <- map_chr(Mosaic$sp_ipms, class)
    if(all(IPM_cl == "ipm") && !is.null(climate)) {
        # no climate needed
        warning(paste0("Because all species are fully integrated on a climate, ",
                       "providing one now is unnecessary"))
        clim_i <- which(IPM_cl == "ipm")[[1]]
        climate <- t(Forest$species[[clim_i]]$IPM$climatic)
        climate <-  as.matrix(bind_cols(climate, t = 1:equil_time))
    } else
        {
        if(any(IPM_cl == "ipm")){
            if(!is.null(climate)){
                warning(
                    paste0("At least one species is fully integrated on a ",
                           "climate, so this climate will be used for simulation"))
            }
            clim_i <- which(IPM_cl == "ipm")[[1]]
            climate <- t(Forest$species[[clim_i]]$IPM$climatic)
        }
        if(inherits(climate, "data.frame")){
            climate <- as.matrix(climate)
        } else if(inherits(climate, "numeric")){
            climate <- t(climate)
        }
        assertMatrix(climate)
        if(nrow(climate) != 1 & nrow(climate) < equil_time){
            stop(paste0("climate matrix is not defined for each time until",
                        " equil_time. This matrix require a row per time or ",
                        "single one."))
        }
        if(nrow(climate) == 1){
            climate <-  as.matrix(bind_cols(climate, t = 1:equil_time))
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
                data = NA_real_, ncol = tlim + 1,
                nrow = x + 3 + x + 1
            )
            colnames(tmp) <- c(paste0("t", 1:(tlim+1))) #, "sp")
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
    Mosaic$sp_ipms <- map(Mosaic$sp_ipms, correction.ipm, correction = correction)

    ## Initiate variables and populations ####
    landscape <- map(nms_plot,
                     ~ init_forest_env(Mosaic, index = .x,
                                       tlim = tlim,  equil_time = equil_time,
                                       SurfEch = SurfEch)
    )
    # save first pop
    landscape <- map(landscape, ~ save_step_env(.x, t = 1))

    # Create sim IPM ####
    sim_clim <- climate[1, , drop = TRUE] # why this line ?
    landscape <- get_step_env(landscape, Mosaic, t= 1,
                              climate, correction, disturb_surv)
    if (verbose) {
        message("Starting while loop. Maximum t = ", equil_time)
    }

    # While tlim & eq ####
    t <- 2
    the <- NA_real_ # real time of simulation ending in case of outbound BA

    while (t < tlim ) {

        sim_clim <- climate[t, , drop = TRUE]
        # TODO continue ici
        landscape <- map(landscape, ~ growth_mortal_env(
            .x, t = 1, harvest = harvest, run_disturb
            ))
        # actual in dev

        #***********************************************************************
        ### Recruitment ####
        #' separer la repro de la compet. regrouper la repro et la diviser par
        #' n patch pour apres appliquer la compet locale.

        #' compute recruitment per plot, mean by species for nplot
        rec <- map(Forest$species, sp_rec.species, sim_clim)

        recrues <- imap(
            rec,
            function(x, .y, basp, banonsp, mesh, SurfEch){
                exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch)
            }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
            mesh = meshs, SurfEch = SurfEch)

        # X <- map2(X, recrues, `+`) # gain time
        X <- sapply(names(X), function(n, x, y) x[[n]] + y[[n]],
                    X, recrues, simplify = FALSE)
        #***********************************************************************

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

    # Format output ####
    # browser()
    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(x / SurfEch, ba[[.y]], bast[[.y]], sum(x) / SurfEch,
          harv[[.y]] / SurfEch, sum(harv[[.y]]) / SurfEch)
    },
    ba = sim_BAsp[t-1,,drop = FALSE],
    bast = sim_BAstand[t-1,,drop = FALSE],
    harv = Harv)

    tmp <- do.call("c", tmp)
    sim_X[, tlim +1] <- tmp

    colnames(sim_X)[tlim + 1] <- paste0("eqt", t-1)


    if (verbose) {
        message("Simulation ended after time ", ifelse(is.na(the), t-1, the))
        message(sprintf(
            "BA stabilized at %.2f with diff of %.2f at time %i",
            sim_BA[t - 1],
            diff(range(sim_BA[max(1, t - equil_dist - 1):(t - 1)])),
            t -1
        ))
        tmp <- Sys.time() - start
        message("Time difference of ", format(unclass(tmp), digits = 3),
                " ", attr(tmp, "units"))
    }
    sim_X <- new_deter_sim(sim_X, mesh = meshs)

    sim_X <- tree_format(sim_X)

    # Return ####
    # return(sim_X)

    res <- NULL
    return(res)
}
