#' Tree size distrib to BA of a plot
#'
#' Build the vector to pass from a tree size distribution for a sampled plot size
#' to a basal area per ha
#'
#' @param mesh Vector of all values between L and U of the IPM
#' (L being the smallest size and U the largest.). Length of this vector is
#' the number of size class in the IPM. num.
#' @param SurfEch Value of plot size surface in ha
#'
#' @return
#' numeric vector
#'
#' @examples
#' Buildct(seq(10,50, 10), 5e-4)
#'
#' @import checkmate
#'
#' @noRd
Buildct <- function(mesh, SurfEch= 0.03){
    # Idiot proof
    assertNumeric(mesh, lower = 0)
    assertNumber(SurfEch, lower = 0)

    ct <- pi*(mesh/2*1e-3)^2 # in m^2
    ct <- ct/SurfEch

    return(ct)
}


#' Deterministic Tree population simulation
#'
#' Simulate a population size state during \eqn{[1, t_{lim}]} times.
#' The state at time \eqn{t+1} is dependent on state at time \eqn{t} and
#' the projection matrix (IPM).
#' The simulation will run until \eqn{t_{lim}} and if the equilibrium is not reached,
#' it will continue. Only simulation in \eqn{[1, t_{lim}]}, and equilibrium
#' state are returned.
#'
#' At each iteration, the basal area is evaluated to select the corresponding
#' IPM matrix.
#'
#' @param Forest Group of species that each contains IPM for deterministic
#' transition for \eqn{Z_{t}} state in a population to \eqn{Z_{t+1}} state.
#' A species is also defined with recruitment and harvest functions, please
#' see \code{\link[matreex]{species}} for more information.
#' @param tlim Number of simulation iterations (years) in the future. single int.
#' @param equil_dist Number of last n time for which the range difference
#' should not exceed \code{equil_diff} during the equilibrium research.
#' single int.
#' @param equil_diff Difference threshold of the basal area under which
#' equilibrium is assumed. single real.
#' @param equil_time Total maximum time simulation allowed in equilibrium
#' research. Must be higher or equal to tlim and equil_dist. single int.
#' @param harvest Choice of harvest rules between default, Uneven and Even.
#' This indicate what settings will be used. See Details.
#' @param targetBA BA value per ha that is targeted when using uneven harvesting.
#' Single numeric in \eqn{m^2}.
#' @param targetRDI RDI value that is targeted when using even harvesting.
#' RDI is the ratio between the number of trees and the maximum number of trees
#' given the self-thinning boundary for the corresponding mean diameter and
#' species.
#' @param targetKg Kg value that is targeted when using even harvesting.
#' Kg is the ratio between mean quadratic diameter of killed trees and mean
#' quadratic diameter of trees before harvesting.
#' @param final_harv Final harvest time used when \code{harvest} is set to
#' "Even". This parameter drives the final cut time for even stands.
#' @param climate Optional, climate matrix if climate variation along time is
#' needed. Climate variation rely on species created with mu_gr class objects.
#' This matrix require as many rows as time steps until equil_time.
#' If the climate does not variate, a single row can given and will be reused.
#' @param disturbance `r lifecycle::badge("experimental")` parameter.
#' @param correction Choice of correction of the IPM between \code{"none"}
#' (default) and \code{"cut"}. The second option set the last column to 0 in the
#' IPM so that no individual can grow outside of the defines classes.
#' @param SurfEch Value of plot size surface in ha
#'
#' @param verbose Print message. FALSE by default
#'
#' @details
#' Basic simulations input are illustrated in the main vignette.
#' The harvesting scenario and theory is explained in the harvesting vignette.
#'
#'
#' @return
#' Data.frame with long tidyverse format : a row for each observation and a
#'  column per variable. Columns are listed below, some may contains NA values,
#'  as for example species when there is a non-specific variable (BA).
#'  \describe{
#'   \item{species}{ Name of the species.}
#'   \item{var}{Variable of interest}
#'   \item{time}{Time step of the simulation. If the equilibrium is the last
#'   time in \code{tlim} input, this time will occur twice in the table.}
#'   \item{mesh}{Mesh class number, from 1 to n class.}
#'   \item{size}{Size corresponding to the mesh class.}
#'   \item{equil}{Logical if this time step is the equilibrium or last step of
#'   simulation}
#'   \item{value}{Numeric values of the variables.}
#'  }
#'
#' The variables are :
#' \describe{
#'  \item{n}{Distribution of density by mesh along time per ha.}
#'  \item{N}{Sum of density per ha. (colSums for n)}
#'  \item{BAsp}{Basal area of the population per ha and species}
#'  \item{BAstand}{Basal area of the population per ha and species when
#'  excluding size class below dth. See Harvesting vignette.}
#'  \item{h}{Distribution of harvest density by mesh along time per ha.}
#'  \item{H}{Sum of harvested density per ha. (colSums for h)}
#' }
#'
#' @import Matrix
#' @import checkmate
#' @import purrr
#' @importFrom dplyr between bind_cols
#'
#' @name sim_deter_forest
#' @export
sim_deter_forest  <- function(Forest,
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
                             disturbance = NULL,
                             correction = "none",
                             SurfEch = 0.03,
                             verbose = FALSE) {
    UseMethod("sim_deter_forest")
}

#' @method sim_deter_forest species
#' @export
sim_deter_forest.species  <- function(Forest,
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
                              disturbance = NULL,
                              correction = "none",
                              SurfEch = 0.03,
                              verbose = FALSE) {
    sim_deter_forest(
        Forest = forest(species = list(Forest)),
        tlim = tlim,

        equil_dist = equil_dist,
        equil_diff = equil_diff,
        equil_time = equil_time,

        harvest = harvest,
        targetBA = targetBA,
        targetRDI = targetRDI,
        targetKg = targetKg,
        final_harv = final_harv,
        climate = climate,
        disturbance = disturbance,
        correction = correction,
        SurfEch = SurfEch,
        verbose = verbose
    )
}

#' @method sim_deter_forest forest
#' @export
sim_deter_forest.forest  <- function(Forest,
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
                                     disturbance = NULL,
                                     correction = "none",
                                     SurfEch = 0.03,
                                     verbose = FALSE) {
    # browser()
    # tlim = 500
    # equil_dist = 500
    # equil_diff = 500
    # equil_time = 500
    # harvest = "default"
    # targetBA = 20
    # targetRDI = 0.9
    # targetKg = 0.9
    # final_harv = 100
    # correction = "none"
    # SurfEch = 0.03
    # verbose = FALSE

    # TEMP dev
    targetRDI <- map_dbl(Forest$species, ~ targetRDI)
    targetKg <- map_dbl(Forest$species, ~ targetKg)
    # TEMP dev

    # Idiot Proof ####
    # validate_forest(Forest) # TEMP dev
    assertCount(tlim)
    assertCount(equil_dist)
    assertNumber(equil_diff)
    assertNumber(equil_time)
    if (equil_time < tlim || equil_time < equil_dist) {
        stop("equil_time must be higher or equal to tlim and equil_dist")
    }

    harvest <- match.arg(harvest)
    assertNumber(targetBA, lower = 0)
    # assertNumber(targetRDI, lower = 0, upper = 1) # FIXME single or species target ?
    # assertNumber(targetKg, lower = 0, upper = 1)
    IPM_cl <- map_chr(Forest$species, ~ class(.x$IPM))
    if(all(IPM_cl == "ipm") && !is.null(climate)) {
        # no climate needed
        warning(paste0("Because all species are fully integrated on a climate, ",
                       "providing one now is unnecessary"))
        clim_i <- which(IPM_cl == "ipm")[[1]]
        climate <- t(Forest$species[[clim_i]]$IPM$climatic)
        climate <-  as.matrix(bind_cols(climate, t = 1:equil_time))
    } else {
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
    nsp <- length(Forest$species)
    disturb_surv <- TRUE

    ## Modify IPM ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
    }
    # correct also decompress integer to double with x * 1e-7 app
    Forest <- correction(Forest, correction = correction)
    meshs <- map(Forest$species, ~ .x$IPM$mesh)
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

            Disturb <- imap(
                map(Forest$species, `[[`, "disturb_fun"),
                function(f, .y, X, sp, disturb, ...){
                    exec(f, X[[.y]], sp[[.y]], disturb, ...)
                }, X = X, sp = Forest$species,
                disturb = disturbance[disturbance$t == t, ],
                qmd = qmd
            )

            X <- map2(X, Disturb, `-`)

        }

        ## Harvest ####
        if(!disturb && t %% Forest$harv_rule["freq"] == 0 && harvest == "Uneven"){
            ### Uneven ####
            BAstandsp <- map2_dbl(X, Forest$species, getBAstand, SurfEch)
            BAstand <- sum(BAstandsp)
            BAcut <- getBAcutTarget(BAstand, targetBA, Pmax, dBAmin )
            pi <- BAstandsp / BAstand
            Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
            targetBAcut <- Hi * BAstandsp

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
                function(f, .y, X, sp, ct){
                    exec(f, X[[.y]], sp[[.y]], ct = ct[[.y]])
                }, X = X, sp = Forest$species, ct = ct
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
            function(x, .y, basp, banonsp, mesh, SurfEch){
                exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch)
            }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
            mesh = meshs, SurfEch = SurfEch )

        # X <- map2(X, recrues, `+`) # gain time
        X <- sapply(names(X), function(n, x, y) x[[n]] + y[[n]], X, recrues,
               simplify = FALSE)

        # browser()
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
    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(x / SurfEch, ba[[.y]], bast[[.y]], sum(x) / SurfEch,
          harv[[.y]] / SurfEch, sum(harv[[.y]]) / SurfEch)
    },
    ba = sim_BAsp[t,,drop = FALSE],
    bast = sim_BAstand[t,,drop = FALSE],
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
    return(sim_X)
}


#' Class of deterministic simulation
#'
#' @param x a matrix.
#' @param mesh mesh size values to be set as attributes.
#'
#' @details Format is specified in \code{\link[matreex]{sim_deter_forest}}
#'
#' @noRd
new_deter_sim <- function(x = matrix(), mesh = NULL){
    assertMatrix(x)
    structure(x, class = c("deter_sim", "matrix"),
              mesh = mesh)

}

#' tree_format generic
#'
#' Format simulation output from sim_deter_forest function to a more tidyverse
#' format (long format) for ggplot2 and filtering.
#'
#' @param x Simulations created with sim_deter_forest
#'
#' @name tree_format
#' @keywords internal
#' @export
tree_format <- function(x){
    UseMethod("tree_format")
}


#' @rdname tree_format
#' @keywords internal
#' @export
tree_format.default <- function(x){
    warning(paste0("tree_format is now deprecated. It is already integrated to",
                   " sim_deter_forest function to simplify the simulation pipeline"))
    x
}


#' @rdname tree_format
#' @importFrom dplyr relocate
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map flatten_dbl
#' @keywords internal
#' @export
tree_format.deter_sim <- function(x){

    mesh <- attributes(x)$mesh %>%
        purrr::map( ~ c(.x, NA, NA, NA, .x, NA)) %>%
        purrr::flatten_dbl()

    if(length(mesh) == 0){
        warning("mesh attribute missing, size column will be composed of NA")
        mesh <- rep(NA_real_, nrow(x))
    }

    # cut rownames in different variables.
    var <- rownames(x)
    pattern <- "([[:alpha:]]+_?[[:alpha:]]*)[.]([[:alpha:]]+)([[:digit:]]*)"
    rnms <- data.frame(
        species = sub(pattern, "\\1", var, perl = TRUE),
        var = sub(pattern, "\\2", var, perl = TRUE),
        mesh = as.numeric(sub(pattern, "\\3", var, perl = TRUE)))

    eq_lgl <- function(x){ # function used in pivot_longer
        x == "eq"
    }

    # pivot_longer all of it.
    res <- cbind(as.data.frame(x), rnms, size = mesh)
    res <- tidyr::pivot_longer(
        res,
        -c("var", "species", "mesh", "size"),
        names_to = c("equil", "time"), names_pattern = "(eq)?t(.*)",
        names_transform = list(time = as.numeric, equil = eq_lgl)
    )

    res <- dplyr::relocate(res, "species", "var", "time", "mesh",
                           "size", "equil", "value")
    return(res)
}
