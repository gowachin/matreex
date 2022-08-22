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
#' @keywords internal
#'
#' @export
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
#' @param correction Choice of correction of the IPM between \code{"none"}
#' (default) and \code{"cut"}. The second option set the last column to 0 in the
#' IPM so that no individual can grow outside of the defines classes.
#' @param SurfEch Value of plot size surface in ha
#' @param delay Number of year delay between the recruitment of an individual
#' and it's inclusion in the IPM. This will enlarge the IPM and add sub diagonal
#' values of 1. # TODO see code{link{treeforce}{delay.ipm}}.
#'
#' @param verbose Print message, used for debugs only. FALSE by default
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
#' @export
sim_deter_forest <- function(Forest,
                             tlim = 3e3,
                             equil_dist = 250,
                             equil_diff = 1,
                             equil_time = 1e4,
                             # Harv = 0.006, # TEMP dev
                             # BAsup = 200, # TEMP dev
                             correction = "none",
                             SurfEch = 0.03,
                             delay = 0,
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
    correction <- match.arg(correction, c("cut", "none"))
    assertNumber(SurfEch, lower = 0)
    assertCount(delay)
    assertLogical(verbose, any.missing = FALSE, len = 1)

    # Initialisation ####
    get_mesh <- function(x){ # TODO : set function outside of here
        return(x$IPM$mesh)
    }
    get_ipm <- function(x, n){
        return(x$IPM$IPM[[n]])
    }
    init_sim <- function(nsp, tlim, mesh){ # TODO : set function outside of here
        res <- vector("list", nsp)
        res <- map2(lengths(mesh), names(mesh), function(x, y) {
            # tmp <- as.data.frame(
            tmp <- matrix(
                data = NA_real_, ncol = tlim + 1, nrow = x + 2 + x + 1
            )
            # )
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
    # QUESTION group these two modif as one ?
    if (delay > 0) { ### Apply Delay ####
        if (verbose) {
            message("apply a IPM delay of ", delay)
        }
        Forest <- delay(Forest, delay = delay)
    }
    ### Apply correction ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
    }
    # correct also decompress integer to double with x * 1e-7 app
    Forest <- correction(Forest, correction = correction)

    ## Create output ####
    sim_X <- init_sim(nsp, tlim, map(Forest$species, get_mesh))
    sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, names(Forest$species))
    ))
    sim_BA <- rep(NA_real_, equil_time)
    sim_BAnonSp <- rep(NA_real_, nsp)

    ## Initiate pop ####
    X <- map2(map(Forest$species, `[[`, "init_pop"),
              map(Forest$species, get_mesh),
              exec, SurfEch = SurfEch)
    Harv <- map(lengths(map(Forest$species, get_mesh)), ~ rep(0, .x))
    ct <- map(map(Forest$species, get_mesh),
              Buildct, SurfEch = SurfEch)
    BAsp <- map(Forest$species, ~ .x$IPM$BA)
    # save first pop
    sim_BAsp[1, ] <- map2_dbl(X, ct, ~ .x %*% .y )
    sim_BA[1] <- sum(sim_BAsp[1,])
    sim_BAnonSp <- map2_dbl( - sim_BAsp[1, ,drop = FALSE], sim_BA[1],  `+`)

    # tmp <- map2(X, sim_BAsp[1, ,drop = TRUE], ~ c(.x, .y, sum(.x)))
    tmp <- imap(X, function(x, .y, ba, harv){
        c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]) )
    }, ba = sim_BAsp[1,, drop = FALSE], harv = Harv )

    sim_X <- map2(sim_X, tmp, ~ `[<-`(.x, 1:nrow(.x), 1, value = .y))

    if (any(map2_lgl(sim_BA[1], BAsp, ~ ! between(.x, min(.y), max(.y))))) {
        stop(paste(
            "Maximum Basal Area reached for this simulation.",
            "This maximum is reached before iteration, check init_pop functions"
        ))
    }

    # Create sim IPM ####
    lower_ba <- map_dbl(BAsp, ~ .x[which(.x == max(.x[.x <= sim_BA[1]]))] )
    higher_ba <- map_dbl(BAsp, ~ .x[which(.x == min(.x[.x > sim_BA[1]]))] )
    low_ba <- map2(Forest$species, lower_ba, get_ipm)
    high_ba <- map2(Forest$species, higher_ba, get_ipm)

    sim_ipm <- lapply(
        seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
            low_ba[[i]] * (1 - (floor(ba) - nipm[i])) + # TODO check if this is correct formula !
                high_ba[[i]] * ( floor(ba)  - nipm[i] )
        }, low_ba, high_ba, sim_BA[1], lower_ba
    )
    if(delay > 0){
        # add sub diag
        sim_ipm <- map( sim_ipm, function(x) {
            for (i in 1:delay) {
                x[i + 1, i] <- 1
            }
            x} )
    }

    if (verbose) {
        message("Starting while loop. Maximum t = ", equil_time)
    }

    # While tlim & eq ####
    t <- 2
    the <- NA_real_ # real time of simulation ending in case of outbound BA
    while (t <= equil_time && (t <= tlim || diff(
        range(sim_BA[max(1, t - 1 - equil_dist):max(1, t - 1)])
    ) > equil_diff)) {

        ## t size distrib ####
        old_X <- X # TEMP dev
        X <- map2(sim_ipm, X, ~ drop( .x %*% .y ) ) # Growth
        Harv <- map2(map(Forest$species, `[[`, "harvest_fun"), X, exec) # Harvest
        X <- map2(X, Harv, `-`)

        recrues <- imap(
            map(Forest$species, `[[`, "recruit_fun"),
            function(x, .y, basp, banonsp, mesh, SurfEch){
                exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch)
            }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
            mesh = map(Forest$species, get_mesh), SurfEch = SurfEch )

        X <- map2(X, recrues, `+`) # Recruitment
        # compute new BA for selecting the right IPM and save values
        sim_BAsp[t, ] <- map2_dbl(X, ct, `%*%`)
        sim_BA[t] <- sum(sim_BAsp[t,])
        sim_BAnonSp <- map2_dbl( - sim_BAsp[t, ,drop = FALSE], sim_BA[t],  `+`)

        # Update X
        if (t <= tlim) {
            # tmp <- map2(X, sim_BAsp[t, ,drop = TRUE], ~ c(.x, .y, sum(.x)))
            tmp <- imap(X, function(x, .y, ba, harv){
                c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]))
            }, ba = sim_BAsp[t,,drop = FALSE], harv = Harv )
            sim_X <- map2(sim_X, tmp, ~ `[<-`(.x, 1:nrow(.x), t, value = .y))
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
        # input : BAsp, sim_BA, delay, Forest, ( t ? )
        # output : sim_ipm
        lower_ba <- map_dbl(BAsp, ~ .x[which(.x == max(.x[.x <= sim_BA[t]]))] )
        higher_ba <- map_dbl(BAsp, ~ .x[which(.x == min(.x[.x > sim_BA[t]]))] )
        low_ba <- map2(Forest$species, lower_ba, get_ipm)
        high_ba <- map2(Forest$species, higher_ba, get_ipm)

        sim_ipm <- lapply(
            seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
                low_ba[[i]] * (1 - (floor(ba) - nipm[i])) + # TODO check if this is correct formula !
                    high_ba[[i]] * ( floor(ba)  - nipm[i] )
            }, low_ba, high_ba, sim_BA[t], lower_ba
        )
        if(delay > 0){
            # add sub diag
            sim_ipm <- map( sim_ipm, function(x) {
                for (i in 1:delay) {
                    x[i + 1, i] <- 1
                }
                x} )
        }
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
    # tmp <- map2(X, sim_BAsp[t-1, ,drop = TRUE], ~ c(.x, .y, sum(.x)))
    tmp <- imap(X, function(x, .y, ba, harv){
        c(x, ba[[.y]], sum(x), harv[[.y]], sum(harv[[.y]]))
    }, ba = sim_BAsp[t-1,,drop = FALSE], harv = Harv )
    sim_X <- map2(sim_X, tmp, ~ `[<-`(.x, 1:nrow(.x), tlim+1, value = .y))
    sim_X <- map(sim_X, ~ {
        `<-`(`[`(colnames(.x), tlim+1), paste0("t", t-1))
        .x
    })
    sim_X <- do.call(rbind, sim_X)


    if (verbose) {
        message("Simulation ended after time ", ifelse(is.null(the), t-1, the))
        message(sprintf(
            "BA stabilized at %.2f with diff of %.2f at time %i",
            sim_BA[t - 1], diff(range(sim_BA[max(1, t - equil_dist - 1):(t - 1)])),
            t -1
        ))
    }

    # Return ####
    return(sim_X)
}

#' Class of deterministic simulation
#'
#' @param x a matrix.
#'
#' @details Format is specified in \code{\link{treeforce}{sim_deter_forest}}
#'
#' @noRd
new_deter_sim <- function(x = matrix()){
    assertMatrix(x)
    structure(x, class = "deter_sim")
}

#' Summary for treeforce package objects.
#'
#' @param object deter_sim class
#' @param ... Ignored
#'
#' @keywords internal
#' @method summary deter_sim
#' @export
summary.deter_sim <- function(object, ...){

    res <- list(
        BA = object["BA",],
        N = object["N",],
        state_eq = head(object[, ncol(object)], nrow(object) - 2)
    )

    res$BA_eq <- tail(res$BA, 1)
    res$time_eq <- sub("t", "", names(res$BA_eq))

    class(res) <- "summary_sim"
    return(res)
}

#' Summary for treeforce package objects.
#'
#' @param object summary_sim class
#' @param ... Ignored
#'
#' @keywords internal
#' @method summary summary_sim
#' @export
summary.summary_sim <- function(object, ...){
    return(object)
}

#' Generics for treeforce classes.
#'
#' @rdname print_treeforce
#'
#' @param x deter_sim class
#' @param ... Ignored
#'
#' @keywords internal
#' @method print deter_sim
#' @export
print.deter_sim <- function(x, ...){
    print(summary(x))

    return(x)
}

#' @rdname print_treeforce
#'
#' @param x summary_sim class
#' @param ... Ignored
#'
#' @keywords internal
#' @method print summary_sim
#' @export
print.summary_sim <- function(x, ...){

    cat("Summary of deterministic simulation :\n")
    cat("Equilibrium reached after time", x$time_eq, "for an initial duration ")
    cat("of", length(x$N) -1, "times\n")
    cat("BA stabilized at", x$BA_eq, "\n")
    cat("Mesh dimension :", length(x$state_eq))

    return(invisible(x))
}

#' tree_format generic
#'
#' @param x any treeforce package object available.
#'
#' @name tree_format
#' @export
tree_format <- function(x){
    UseMethod("tree_format")
}


#' @rdname tree_format
#' @importFrom dplyr mutate relocate
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data
#' @importFrom tibble rownames_to_column
#' @importFrom purrr %>%
#' @export
tree_format.deter_sim <- function(x) {
  summ <- summary(x)
  time <- summ$time_eq

  equi <- as.data.frame(summ["state_eq"]) %>%
    rownames_to_column("m") %>%
    mutate(m = as.numeric(sub("m", "", .data$m))) %>%
    pivot_longer(cols = -.data$m) %>%
    mutate(t = NA_real_) %>%
    relocate(.data$name, .data$t, .data$m, .data$value)

  along <- as.data.frame(summ[c("BA", "N")]) %>%
    rownames_to_column("t") %>%
    mutate(t = as.numeric(sub("t", "", .data$t))) %>%
    pivot_longer(cols = -.data$t) %>%
    mutate(m = NA_real_) %>%
    relocate(.data$name, .data$t, .data$m, .data$value)

  res <- rbind(equi, along)

  res[which(res$t == time), "name"] <- paste0(
    res[which(res$t == time), "name", drop = TRUE], "eq"
  )

  return(res)
}

#' @rdname tree_format
#' @export
tree_format.summary_sim <- tree_format.deter_sim
