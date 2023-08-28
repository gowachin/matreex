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
library(dplyr)
library(ggplot2)
data("fit_Abies_alba")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Abies_alba", select = -sp)
# N here is a climate defined in Kunstler et al 2021.
# N == 2 is the optimum climate for the species.
# see ?climate_species for more info.
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
options(W_matreex_edist = FALSE)
Abies_sp <- species(IPM = Abies_ipm, init_pop = def_initBA(30))
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
comp_time <- 500

regional_Abies <- forest(
    species = list(Abies = Abies_spE),
    regional_abundance = equil_BA,
    migration_rate = c(Abies = 0.05)
)

load_all()
regional_sim <- sim_deter_forest(
    regional_Abies,
    tlim = comp_time,
    equil_time = comp_time, equil_dist = 10,
    SurfEch = 0.03,
    verbose = TRUE
)

isol_sim <- sim_deter_forest(
    forest(
        species = list(Abies = Abies_spE),
        regional_abundance = list(Abies = equil_dist),
        migration_rate = c(Abies = 0)
    ),
    tlim = comp_time,
    equil_time = comp_time, equil_dist = 10,
    SurfEch = 0.03,
    verbose = TRUE
)


regional_pool <- dplyr::bind_rows(regional002 = regional_sim, isolation = isol_sim,
                                  .id = "migration")

regional_pool %>%
    dplyr::filter(var %in% c("BAsp", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = migration)) +
    geom_line(linewidth = .4) +
    facet_wrap(~ var, scales =  "free_y") +
    NULL

#' I guess it works, but questions :
#' - What basal area use for regional pool
#' - Migration rate values ?
#' - SurfEch for the regional pool ?

# Testing 2 species simulations
AbFa_sim %>%
    dplyr::filter(var %in% c("BAsp", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linewidth = .4) +
    facet_wrap(~ var, scales =  "free_y") +
    NULL

e2_dist <- dplyr::filter(AbFa_sim, equil, var == "n") %>%
    group_by(species) %>%
    dplyr::group_split() %>% map(pull, value)
e2_BA <- dplyr::filter(AbFa_sim, equil, var == "BAsp") %>% dplyr::pull(value)


Abies_spEn <- Abies_sp
Abies_spEn$init_pop <- def_init_k(e2_dist[[1]])
Fagus_spEn <- Fagus_sp
Fagus_spEn$init_pop <- def_init_k(e2_dist[[2]])


regional_AbFa <- forest(
    species = list(Abies = Abies_spEn, Fagus = Fagus_spEn),
    regional_abundance = e2_BA,
    migration_rate = c(Abies = 0.1, Fagus = 0.1)
)

sp2_sim <- sim_deter_forest(
    regional_AbFa,
    tlim = comp_time,
    equil_time = comp_time, equil_dist = 10,
    SurfEch = 0.03,
    verbose = TRUE
)

sp2_sim %>%
    dplyr::filter(var %in% c("BAsp", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linewidth = .4) +
    facet_wrap(~ var, scales =  "free_y") +
    NULL


# Testing invasive species
Fagus_spE0 <- Fagus_sp
Fagus_spE0$init_pop <- def_init_k(e2_dist[[2]] * 10e-7)


invasive_AbFa <- forest(
    species = list(Abies = Abies_spE, Fagus = Fagus_spE0),
    regional_abundance = e2_BA,
    migration_rate = c(Abies = 0.1, Fagus = 0.1)
)

invas_sim <- sim_deter_forest(
    invasive_AbFa,
    tlim = comp_time,
    equil_time = comp_time, equil_dist = 10,
    SurfEch = 0.03,
    verbose = TRUE
)

invas_sim %>%
    dplyr::filter(var %in% c("BAsp", "N"), !equil) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linewidth = .4) +
    facet_wrap(~ var, scales =  "free_y") +
    NULL

## Mosaic forest ####

#' The most important thing is to prevent RAM overload


forests <- list(A = regional_Abies,
                B = invasive_AbFa)

map(forests, ~ .x$info$species) %>% flatten_chr() %>% unique()

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

## Simulations ####
sim_deter_mosaic <- function(Mosaic,
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


    # browser()
    tlim = 500
    equil_dist = 1
    equil_diff = 1
    equil_time = 500
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


    # TODO continue ici
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
    disturb_surv <- TRUE

    ## Modify IPM ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
    }
    # correct also decompress integer to double with x * 1e-7 app
    Mosaic$sp_ipms <- map(Mosaic$sp_ipms, correction.ipm, correction = correction)

    meshs <- map(Forest$sp_ipms, ~ .x$IPM$mesh)
    stand_above_dth <- map2(meshs, Forest$species, ~ .x > .y$harv_lim["dth"])
    # delay <- map(Forest$species, ~ as.numeric(.x$IPM$info["delay"]))

    ## Create output ####
    sim_X <- init_sim(nsp, tlim, meshs)
    sim_X <- do.call("rbind", sim_X)
    sim_BAstand <- sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, names(Forest$species))
    ))
    sim_BA <- rep(NA_real_, equil_time)
    sim_BAnonSp <- rep(NA_real_, nsp)

    ## Initiate pop ####
    X <- map2(map(Forest$species, `[[`, "init_pop"),
              meshs,
              exec, SurfEch = SurfEch)
    Harv <- map(lengths(meshs), ~ rep(0, .x))
    ct <- map(meshs, Buildct, SurfEch = SurfEch)

    BAsp <- map(Forest$species, ~ .x$IPM$BA)
    # save first pop
    sim_BAsp[1, ] <- map2_dbl(X, ct, ~ .x %*% .y )
    standX <- map2(X, stand_above_dth, `*`)
    sim_BAstand[1, ] <- map2_dbl(standX, ct, `%*%`)
    sim_BA[1] <- sum(sim_BAsp[1,])
    sim_BAnonSp <- map2_dbl( - sim_BAsp[1, ,drop = FALSE], sim_BA[1],  `+`)

    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(x / SurfEch , ba[[.y]], bast[[.y]], sum(x) / SurfEch,
          harv[[.y]] / SurfEch, sum(harv[[.y]]) / SurfEch )
    },
    ba = sim_BAsp[1,, drop = FALSE],
    bast = sim_BAstand[1,,drop = FALSE],
    harv = Harv )

    tmp <- do.call("c", tmp)
    sim_X[, 1] <- tmp

    if (any(map2_lgl(sim_BA[1], BAsp, ~ ! between(.x, min(.y), max(.y -1))))) {
        stop(paste(
            "Border Basal Area reached for this simulation.",
            "This maximum is reached before iteration, check init_pop functions"
        ))
    }

    if(regional){
        if(verbose){
            message("Simulation with regional pool")
        }
        reg_ba <- Forest$regional_abundance
        reg_banonsp <- sum(reg_ba) - reg_ba
        migrate <- Forest$migration_rate
    } else {
        migrate <- map(X, ~ 0)
    }

    # Create sim IPM ####
    start_clim <- climate[1, , drop = TRUE]

    sim_ipm <- map(Forest$species, ~ get_step_IPM(
        x = .x$IPM, BA = sim_BA[1], climate = start_clim, sim_corr = correction,
        IsSurv = disturb_surv
    ))

    if (verbose) {
        message("Starting while loop. Maximum t = ", equil_time)
    }

    # While tlim & eq ####
    t <- 2
    the <- NA_real_ # real time of simulation ending in case of outbound BA
    # Harv cst
    alpha <- Forest$harv_rule["alpha"]
    Pmax <- Forest$harv_rule["Pmax"]
    dBAmin <- Forest$harv_rule["dBAmin"]
    disturb <- FALSE

    while (t < tlim || (t <= equil_time && (t <= tlim || diff(
        range(sim_BA[max(1, t - 1 - equil_dist):max(1, t - 1)])
    ) > equil_diff))) {


        ## t size distrib ####
        X <- map2(X, sim_ipm, ~ drop( .y %*% .x ) )# Growth

        ## Disturbance ####
        if(run_disturb && t_disturb[t]){
            disturb <- TRUE

            if (verbose) {
                message(sprintf(
                    "time %i | Disturbance : %s I = %.2f",
                    t, disturbance[disturbance$t == t, "type"],
                    disturbance[disturbance$t == t, "intensity"]
                )
                )
            }

            qmd <- QMD(size = unlist(meshs), n = unlist(X))
            # TODO remove unborn size from X before computations

            # TODO compute percentage of coniferous (relative share in number of stems)
            total_stem <- purrr::reduce(X, sum, .init = 0)
            sp_stem <- map_dbl(X, ~ sum(.x) / total_stem)
            perc_coni <- sum(sp_stem[names(types[types == "Coniferous"])])
            # browser()

            Disturb <- imap(
                map(Forest$species, `[[`, "disturb_fun"),
                function(f, .y, X, sp, disturb, ...){
                    exec(f, X[[.y]], sp[[.y]], disturb, ...)
                }, X = X, sp = Forest$species,
                disturb = disturbance[disturbance$t == t, ],
                qmd = qmd, perc_coni = perc_coni
            )

            X <- map2(X, Disturb, `-`)

        }

        ## Harvest ####
        if(!disturb && t %% Forest$harv_rule["freq"] == 0 &&
           harvest %in% c("Uneven", "Favoured_Uneven")){
            ### Uneven ####
            BAstandsp <- map2_dbl(X, Forest$species, getBAstand, SurfEch)
            BAstand <- sum(BAstandsp)
            BAcut <- getBAcutTarget(BAstand, targetBA, Pmax, dBAmin )

            sfav <- sum(Forest$favoured_sp)
            if( harvest == "Favoured_Uneven" && (sfav == 0 || sfav == length(Forest$favoured_sp))){
                print('!!!!!!!!!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!!!!!!')
                # warning("No species are favoured in the forest object, harvest mode 'Favoured_Uneven' is replaced with 'Uneven'")
                harvest <- "Uneven"
            }

            if(harvest == "Uneven"){
                pi <- BAstandsp / BAstand
                Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
                targetBAcut <- Hi * BAstandsp
            } else { # Favoured_Uneven
                p_fav <- sum(BAstandsp[Forest$favoured_sp])/BAstand
                # cat(p_fav)
                if(p_fav > 0.5){
                    # ici qu'il faut modifier en fait !
                    # cat(" - let's fav \n")
                    Hi <- BAcut / BAstand
                }  else {
                    # cat(" \n")
                    pi <- ifelse(Forest$favoured_sp, p_fav, 1-p_fav)
                    Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
                }
                targetBAcut <- Hi * BAstandsp
            }
            # browser()

            Harv <- imap(
                map(Forest$species, `[[`, "harvest_fun"),
                function(f, .y, X, sp, bacut, ct, ...){
                    exec(f, X[[.y]], sp[[.y]],
                         targetBAcut = bacut[[.y]],
                         ct = ct[[.y]], ...)
                }, X = X, sp = Forest$species, bacut = targetBAcut,
                ct = ct, t = t
            )

            X <- map2(X, Harv, `-`)
        } else if(!disturb && harvest == "Even"){
            ### Even ####
            if(t %% final_harv == 0){
                Harv <- X
                X <- map2(map(Forest$species, `[[`, "init_pop"),
                          meshs,
                          exec, SurfEch = SurfEch)
            } else if(t %% Forest$harv_rule["freq"] == 0){
                Harv <- imap(
                    map(Forest$species, `[[`, "harvest_fun"),
                    function(f, .y, X, sp, tRDI, tKg, ct, ...){
                        exec(f, X[[.y]], sp[[.y]],
                             targetRDI = tRDI[[.y]],
                             targetKg = tKg[[.y]],
                             ct = ct[[.y]],
                             ...)
                    }, X = X, sp = Forest$species, tRDI = targetRDI,
                    tKg = targetKg, ct = ct, t = t, SurfEch = SurfEch
                )
                X <- map2(X, Harv, `-`)

            } else {
                Harv <- map(meshs, ~ rep(0, length(.x)))
            }
        } else if (!disturb && t %% Forest$harv_rule["freq"] == 0 && harvest == "default") {
            ### Nothing ####
            Harv <- imap(
                map(Forest$species, `[[`, "harvest_fun"),
                function(f, .y, X, sp, ct, ...){
                    exec(f, X[[.y]], sp[[.y]], ct = ct[[.y]], ...)
                }, X = X, sp = Forest$species, ct = ct, t = t, SurfEch = SurfEch
            )

            X <- map2(X, Harv, `-`)
        } else if(disturb){
            Harv <- Disturb
            disturb <- FALSE
        } else {
            Harv <- map(meshs, ~ rep(0, length(.x)))
        }

        ### Recruitment ####
        sim_clim <- climate[t, , drop = TRUE]
        rec <- map(Forest$species, sp_rec.species, sim_clim)

        recrues <- imap(
            rec,
            function(x, .y, basp, banonsp, mesh, SurfEch, mig){
                exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch) * (1 - mig[[.y]])
            }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
            mesh = meshs, SurfEch = SurfEch, mig = migrate )

        if(regional){
            reg_recrues <- imap(
                rec,
                function(x, .y, basp, banonsp, mesh, SurfEch, mig){
                    exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch) * mig[[.y]]
                }, basp = reg_ba, banonsp = reg_banonsp,
                mesh = meshs, SurfEch = SurfEch, migrate )
        } else {
            reg_recrues <- map(recrues, ~ .x * 0)
        }

        # X <- map2(X, recrues, `+`) # gain time
        X <- sapply(names(X), function(n, x, y, z) x[[n]] + y[[n]] + z[[n]],
                    X, recrues, reg_recrues, simplify = FALSE)

        ## Save BA ####
        # compute new BA for selecting the right IPM and save values
        sim_BAsp[t, ] <- map2_dbl(X, ct, `%*%`)
        standX <- map2(X, stand_above_dth, `*`)
        sim_BAstand[t, ] <- map2_dbl(standX, ct, `%*%`)
        sim_BA[t] <- sum(sim_BAsp[t,])
        sim_BAnonSp <- map2_dbl( - sim_BAsp[t, ,drop = FALSE], sim_BA[t],  `+`)

        # Update X and extract values per ha
        if (t <= tlim) {
            tmp <- imap(X, function(x, .y, ba, bast, harv){
                c(x / SurfEch, ba[[.y]], bast[[.y]], sum(x) / SurfEch,
                  harv[[.y]] / SurfEch, sum(harv[[.y]]) / SurfEch)
            },
            ba = sim_BAsp[t,,drop = FALSE],
            bast = sim_BAstand[t,,drop = FALSE],
            harv = Harv)

            tmp <- do.call("c", tmp)
            sim_X[, t] <- tmp
        }


        ## Stop loop if BA larger than LIPM largest BA ####
        if (any(map2_lgl(sim_BA[t], BAsp, ~ ! between(.x, min(.y), max(.y))))) {
            warning("Maximum Basal Area reached for this simulation.")
            # TODO  say which species reached BA limit !
            the <- t
            break()
        }

        ## Get sim IPM ####
        # Is there a disturbance ?
        if(run_disturb && t < equil_time){ # IDEA rewrite this ?
            if(t_disturb[t+1]){
                disturb_surv <- disturbance[disturbance$t == t+1, "IsSurv"]
            } else {
                disturb_surv <- TRUE
            }
        }

        sim_ipm <- map(Forest$species, ~ get_step_IPM(
            x = .x$IPM, BA = sim_BA[t], climate = sim_clim, sim_corr = correction,
            IsSurv = disturb_surv
        ))

        ## Loop Verbose ####
        if (t %% 500 == 0 && verbose) {
            message(sprintf(
                "time %i | BA diff : %.2f",
                t, diff(range(sim_BA[max(1, t - equil_dist):t]))
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
