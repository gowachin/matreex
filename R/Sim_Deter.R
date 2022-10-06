#' Tree size distrib to BA of a plot
#'
#' Build the vector to pass from a tree size distribution to
#' a basal area computed for a given plot size
#'
#' @param mesh Vector of all values between L and U of the IPM
#' (L being the smallest size and U the largest.). Length of this vector is
#' the number of size class in the IPM. num.
#' @param SurfEch Value of plot size surface in ha
#'
#' @return
#' numeric vector of basal area computed depending on plot size
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
#' A species is also defined with recruitment and harvest functions.
#' test # TODO link to class species
#' @param tlim Number of simulation iterations (years) in the future. single int.
#' @param equil_dist Number of last n time for which the range difference
#' should not exceed \code{equil_diff} during the equilibrium research.
#' single int.
#' @param equil_diff Difference threshold of the basal area under which
#' equilibrium is assumed. single real.
#' @param equil_time Total maximum time simulation allowed in equilibrium
#' research. Must be higher or equal to tlim and equil_dist. single int.
#' @param targetBA BA value that is targetted by the Harvest module.
#' Single numeric in \eqn{m^2}.
#' @param correction Choice of correction of the IPM between \code{"none"}
#' (default) and \code{"cut"}. The second option set the last column to 0 in the
#' IPM so that no individual can grow outside of the defines classes.
#' @param SurfEch Value of plot size surface in ha
#'
#' @param verbose Print message. FALSE by default
#'
#'
#' @return
#' Matrix of the population states for time in \eqn{[1, t_{lim}]},
#' plus a last column that is the state at equilibrium.
#' In row are the states of the population, plus two last rows :
#' N, the total number of individual in the population at
#' time \eqn{t} and BA, the basal area of the population.
#' \code{colnames} are set as tn with n the time of simulation.
#' \code{rownames} are labelled mi with i the different state of the population.
#'
#' @import Matrix
#' @import checkmate
#' @import purrr
#' @importFrom dplyr between
#'
#' @name sim_deter_forest
#' @export
sim_deter_forest  <- function(Forest,
                             tlim = 3e3,
                             equil_dist = 250,
                             equil_diff = 1,
                             equil_time = 1e4,
                             targetBA = 20,
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
                              targetBA = 20,
                              correction = "none",
                              SurfEch = 0.03,
                              verbose = FALSE) {
    sim_deter_forest(
        Forest = forest(species = list(Forest)),
        tlim = tlim,

        equil_dist = equil_dist,
        equil_diff = equil_diff,
        equil_time = equil_time,

        targetBA = targetBA,
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
                                     targetBA = 20,
                                     correction = "none",
                                     SurfEch = 0.03,
                                     verbose = FALSE) {


    # Idiot Proof ####
    validate_forest(Forest)
    assertCount(tlim)
    assertCount(equil_dist)
    assertNumber(equil_diff)
    assertNumber(equil_time)
    if (equil_time < tlim || equil_time < equil_dist) {
        stop("equil_time must be higher or equal to tlim and equil_dist")
    }

    assertNumber(targetBA, lower = 0)
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
                nrow = x + 2 + x + 1
            )
            colnames(tmp) <- c(paste0("t", 1:(tlim+1))) #, "sp")
            rownames(tmp) <- c(paste0(y, ".m", 1:x), paste0(y, c(".BAsp", ".N")),
                               paste0(y, ".h", 1:x), paste0(y,".H"))
            return(tmp)
        })
        names(res) <- names(mesh)
        return(res)
    }

    nsp <- length(Forest$species)

    ## Modify IPM ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
    }
    # correct also decompress integer to double with x * 1e-7 app
    Forest <- correction(Forest, correction = correction)
    meshs <- map(Forest$species, ~ .x$IPM$mesh)
    delay <- map(Forest$species, ~ as.numeric(.x$IPM$info["delay"]))

    ## Create output ####
    sim_X <- init_sim(nsp, tlim, meshs)
    sim_X <- do.call("rbind", sim_X)
    sim_BAsp <- as.data.frame(matrix(
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
    sim_BA[1] <- sum(sim_BAsp[1,])
    sim_BAnonSp <- map2_dbl( - sim_BAsp[1, ,drop = FALSE], sim_BA[1],  `+`)

    tmp <- imap(X, function(x, .y, ba, harv){
        c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]) )
    }, ba = sim_BAsp[1,, drop = FALSE], harv = Harv )

    tmp <- do.call("c", tmp)
    sim_X[, 1] <- tmp

    if (any(map2_lgl(sim_BA[1], BAsp, ~ ! between(.x, min(.y), max(.y))))) {
        stop(paste(
            "Maximum Basal Area reached for this simulation.",
            "This maximum is reached before iteration, check init_pop functions"
        ))
    }

    # Create sim IPM ####

    low_id <- map_dbl(BAsp, ~ which(.x == max(.x[.x <= sim_BA[1]])) )
    high_id <- map_dbl(BAsp, ~ which(.x == min(.x[.x > sim_BA[1]])) )
    lower_ba <- map2_dbl(BAsp, low_id, ~ .x[.y] )
    higher_ba <- map2_dbl(BAsp, high_id, ~ .x[.y] )

    low_ba <- map2(Forest$species, low_id, ~ .x$IPM$IPM[[.y]])
    high_ba <- map2(Forest$species, high_id, ~ .x$IPM$IPM[[.y]])

    sim_ipm <- lapply(
        seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
            low_ba[[i]] * (1 - (floor(ba) - nipm[i])) +
                high_ba[[i]] * ( floor(ba)  - nipm[i] )
        }, low_ba, high_ba, sim_BA[1], lower_ba
    )

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

    while (t < tlim || (t <= equil_time && (t <= tlim || diff(
        range(sim_BA[max(1, t - 1 - equil_dist):max(1, t - 1)])
    ) > equil_diff))) {

        ## t size distrib ####
        X <- map2(X, sim_ipm, ~ drop( .y %*% .x ) )# Growth

        ### Harvest ####
        if(t %% Forest$harv_rule["freq"] == 0){
            BAstandsp <- map2_dbl(X, Forest$species, getBAstand, SurfEch)
            BAstand <- sum(BAstandsp)
            BAcut <- getBAcutTarget(BAstand, targetBA, Pmax, dBAmin )
            pi <- BAstandsp / BAstand
            Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
            targetBAcut <- Hi * BAstandsp

            Harv <- imap(
                map(Forest$species, `[[`, "harvest_fun"),
                function(f, .y, X, sp, bacut, ct){
                    exec(f, X[[.y]], sp[[.y]], bacut[[.y]], ct[[.y]])
                }, X = X, sp = Forest$species, bacut = targetBAcut, ct = ct
            )

            X <- map2(X, Harv, `-`)
        } else {
            Harv <- map(meshs, ~ rep(0, length(.x)))
        }

        ### Recruitment ####
        recrues <- imap(
            map(Forest$species, `[[`, "recruit_fun"),
            function(x, .y, basp, banonsp, mesh, SurfEch){
                exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch)
            }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
            mesh = meshs, SurfEch = SurfEch )

        X <- map2(X, recrues, `+`)

        ## Save BA ####
        # compute new BA for selecting the right IPM and save values
        sim_BAsp[t, ] <- map2_dbl(X, ct, `%*%`)
        sim_BA[t] <- sum(sim_BAsp[t,])
        sim_BAnonSp <- map2_dbl( - sim_BAsp[t, ,drop = FALSE], sim_BA[t],  `+`)

        # Update X
        if (t <= tlim) {
            tmp <- imap(X, function(x, .y, ba, harv){
                c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]))
            }, ba = sim_BAsp[t,,drop = FALSE], harv = Harv )

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
        # IDEA make a function for this
        # input : BAsp, sim_BA, delay, Forest, t
        # output : sim_ipm
        low_id <- map_dbl(BAsp, ~ which(.x == max(.x[.x <= sim_BA[t]])) )
        high_id <- map_dbl(BAsp, ~ which(.x == min(.x[.x > sim_BA[t]])) )
        lower_ba <- map2_dbl(BAsp, low_id, ~ .x[.y] )
        higher_ba <- map2_dbl(BAsp, high_id, ~ .x[.y] )

        low_ba <- map2(Forest$species, low_id, ~ .x$IPM$IPM[[.y]])
        high_ba <- map2(Forest$species, high_id, ~ .x$IPM$IPM[[.y]])

        sim_ipm <- lapply(
            # NOTE : bottleneck of the function because of sum of sparse matrix !
            # But using as.matrix is even longer so long live the sparse matrix !
            seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
                low_ba[[i]] * (1 - (floor(ba) - nipm[i])) +
                    high_ba[[i]] * ( floor(ba)  - nipm[i] )
            }, low_ba, high_ba, sim_BA[t], lower_ba
        )
        # eof idea

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
    tmp <- imap(X, function(x, .y, ba, harv){
        c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]))
    }, ba = sim_BAsp[t-1,,drop = FALSE], harv = Harv )
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

    # Return ####
    return(sim_X)
}


#' Class of deterministic simulation
#'
#' @param x a matrix.
#' @param mesh mesh size values to be set as attributes.
#'
#' @details Format is specified in \code{\link{treeforce}{sim_deter_forest}}
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
#' @export
tree_format <- function(x){
    UseMethod("tree_format")
}


#' @rdname tree_format
#' @importFrom dplyr mutate relocate
#' @importFrom tidyr pivot_longer separate
#' @importFrom rlang .data
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map
#' @keywords internal
#' @export
tree_format.deter_sim <- function(x){

    mesh <- attributes(x)$mesh %>%
        map( ~ c(.x, NA, NA, .x, NA)) %>%
        purrr::flatten_dbl()

    if(length(mesh) == 0){
        warning("mesh attribute missing, size column will be composed of NA")
        mesh <- rep(NA_real_, nrow(x))
    }

    res <- x %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "var") %>%
        mutate(size = mesh) %>%
        tidyr::pivot_longer( -c(.data$var, .data$size), names_to = "time") %>%
        mutate( species = sub( "\\..*$", "", .data$var, perl = TRUE)) %>%
        mutate( var = sub("^.*\\.", "", .data$var, perl = TRUE)) %>%
        dplyr::mutate(
            equil = grepl("eq", .data$time),
            time = as.numeric(
                sub("([[:alpha:]]+)([[:digit:]]+)", "\\2", .data$time,
                    perl = TRUE)
            )) %>%
        dplyr::mutate( mesh = as.numeric(
            sub("([[:alpha:]]+)([[:digit:]]*)", "\\2", .data$var, perl = TRUE)
        )) %>%
        dplyr::mutate( var = sub("([[:alpha:]]+)([[:digit:]]*)",
                                 "\\1", .data$var, perl = TRUE)) %>%
        dplyr::relocate(.data$species, .data$var, .data$time, .data$mesh,
                        .data$size, .data$equil, .data$value)
    return(res)
}
