#' Tree size distrib to BA of a plot
#'
#' Build the vector to pass from a tree size distribution to
#' a basal area computed for a given plot size
#'
#' @param mesh Vector of all values between L and U of the IPM
#' (L being the smallest size and U the largest.). Length of this vector is
#' the number of size class in the IPM. num.
#' @param SurfEch Value of plot size surface in m^2
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

    ct <- t(pi*(mesh/2*1e-3)^2) # in m^2
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
#' @param Model IPM for deterministic transition for \eqn{Z_{t}} state
#' in a population to \eqn{Z_{t+1}} state. ipm object class # TODO link to class ipm
#' @param Xini Population state at T0. Must have the same mesh as Model.
#'  pop_state object class # TODO : link to class pop_state
#' @param tlim Number of simulation iterations (years) in the future. single int.
#' @param equil_dist Number of last n time for which the range difference
#' should not exceed \code{equil_diff} during the equilibrium research.
#' single int.
#' @param equil_diff Difference threshold of the basal area under which
#' equilibrium is assumed. single real.
#' @param equil_time Total maximum time simulation allowed in equilibrium
#' research. Must be higher or equal to tlim and equil_dist. single int.
#' @param Harv Percentage of harvest in the population per year.
#' @param BAsup Maximum basal area per hectare. This should be inferior or equal
#' to the BAsup of the Model. single int.
#' @param correction Choice of correction of the IPM between \code{"none"}
#' (default) and \code{"cut"}. The second option set the last column to 0 in the
#' IPM so that no individual can grow outside of the defines classes.
#' @param SurfEch Value of plot size surface in \eqn{m^2}
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
#'
#' @export
sim_deter_forest <- function(Model, # require the mesh, LIPM and RecFun only
                             Xini,
                             tlim = 3e3,
                             equil_dist = 250,
                             equil_diff = 1,
                             equil_time = 1e4,
                             Harv = 0.006,
                             BAsup = 200,
                             correction = "none",
                             SurfEch = 0.03,
                             delay = 0,
                             verbose = FALSE) {

    # Idiot Proof ####
    # assertClass(Model, "ipm") # TEMP write this function and class, this is an IPM with mesh and
    # assertClass(Xini, "pop_state") # TEMP write this function and class, this is Xt !
    # TODO check if mesh are the same for pop_state and ipm !
    assertCount(tlim)
    assertCount(equil_dist)
    assertNumber(equil_diff)
    assertNumber(equil_time)
    if (equil_time < tlim || equil_time < equil_dist) {
        stop("equil_time must be higher or equal to tlim and equil_dist")
    }
    assertNumber(Harv, lower = 0, upper = 1)
    assertCount(BAsup)
    match.arg(correction, c("cut", "none"))
    assertNumber(SurfEch, lower = 0)
    assertCount(delay)
    assertLogical(verbose, any.missing = FALSE, len = 1)

    # Initialisation ####
    # Apply delay
    if (delay > 0) {
        if (verbose) {
            message("apply a IPM delay of ", delay)
        }
        class(Model) <- "ipm" # TEMP this class will be created one day !
        Model <- delay(Model, delay = delay)
    }

    # link all Model elements to local env
    # list2env(Model, environment())
    # HACK this was a good idea but no visible binding note
    meshpts <- Model$meshpts
    LIPM <- Model$LIPM
    RecFun <- Model$RecFun

    if(length(Xini) != length(meshpts)){
        stop(paste0("Length of population (", length(Xini), ") ",
             "differs from IPM size (", length(meshpts), "), ",
             "please check dimensions. Is Xini delayed ?"))
    }

    # HACK : this is because data is stored as integer...
    # so yeah I guess I'll need to add an attribute compressed in ipm class
    if (max(LIPM[[1]]) > 1e6){
        if(verbose){
            message("Decompress the IPM matrices.")
        }
        for (i in 1:length(LIPM)){
            LIPM[[i]] <- LIPM[[i]] * 1e-7
        }
    }

    ## Modify IPM for eviction limite ####
    if (correction == "cut") {
        if (verbose) {
            message("apply a IPM cut correction")
        }
        for (i in 1:length(LIPM)) {
            LIPM[[i]][, length(meshpts)] <- 0
            LIPM[[i]][length(meshpts), ] <- 0
        }
    }
    # browser()
    ct <- Buildct(meshpts, SurfEch = SurfEch)
    b <- c(rep(1 / 2, 2), numeric(length(meshpts) - 2))

    # TODO : initialize X outside of sim_deter_forest
    # if (missing(Xini) && IsServ == 0){
    #     x <- exp(runif(1,-.005,.005)*meshpts)
    #     x[which(rbinom(length(meshpts),1,runif(1,.6,.9)) == 1)] <- 0
    #     if (delay>0){x <- c(rep(0, delay), x)}
    #     x <- x / (as.numeric(ct %*% x)) * RNG
    # }else{
    x <- Xini
    # }

    ## Find starting IPM with BA ####
    ba <- floor(as.numeric(ct %*% x)) # IDEA : replace with sum(ct * x)
    NIPM <- max(1, min((BAsup - 1), as.integer(ba))) # IDEA : why integer ?

    # FIXME use BA list from the ipm class for selecting an IPM !
    # This could be nice to limit the number of IPM we use !
    # This could also limit the usage of BAsup value.
    P <- (LIPM[[NIPM]] * (1 - (ba - NIPM)) +
              LIPM[[NIPM + 1]] * (ba - NIPM)) * (1 - Harv)
    if (delay > 0) {
        # add sub diag
        for (i in 1:delay) {
            P[i + 1, i] <- 1
        }
    }

    BA <- rep(NA_real_, equil_time)
    XT <- matrix(NA_real_,
                 ncol = tlim + 1,
                 nrow = length(meshpts) + 2)

    # NOTE : maybe save initial state x in the table ?

    if (verbose) {
        message("Starting while loop. Maximum t = ", equil_time)
    }

    # While tlim & eq ####
    t <- 1
    the <- NULL # real time of simulation ending in case of outbound BA
    while (t <= equil_time && (t <= tlim || diff(
        range(BA[max(1, t - 1 - equil_dist):max(1, t - 1)])
    ) > equil_diff)) {

        # compute new size repartition
                 # QUESTION is this BA below ??
        x1 <- P %*% x + exp(RecFun(ct %*% x)) * SurfEch / (300 * 1e-4) * b
        x <- x1 # IDEA why not drop this variable x1 ?
        # compute new BA for selecting the right IPM
        ba <- as.numeric(ct %*% x)
        NIPM <- max(1, as.integer(floor(ba)))
        ## Stop loop if BA larger than LIPM largest BA ####
        if (is.na(NIPM) || NIPM >= (length(LIPM) - 1)) {
            warning("Maximum Basal Area reached for this simulation.")
            the <- t
            t <- equil_time + 1 # jump to the end
        }

        P <- (LIPM[[NIPM]] * (1 - (ba - NIPM)) +
                  LIPM[[NIPM + 1]] * (ba - NIPM)) * (1 - Harv)
        if (delay > 0) {
            for (i in 1:delay) {
                P[i + 1, i] <- 1
            }
        }

        BA[t] <- ba
        if (t <= tlim) {
            XT[, t] <- c(as.numeric(x), sum(x), ba)
        }

        ## Verbose ####
        if (t %% 500 == 0 && verbose) {
            message(sprintf(
                "time %i | BA diff : %.2f",
                t, diff(range(BA[max(1, t - equil_dist):t]))
            ))
        }
        t <- t + 1
    }

    XT[, ncol(XT)] <- c(as.numeric(x), sum(x), ba)
    colnames(XT) <- paste0("t", c(1:tlim, t - 1))
    rownames(XT) <- c(paste0("m", 1:length(meshpts)), "N", "BA")
    XT <- new_deter_sim(XT)

    if (verbose) {
        message("Simulation ended after time ", ifelse(is.null(the), t, the))
        message(sprintf(
            "BA stabilized at %.2f with diff of %.2f",
            BA[t - 1], diff(range(BA[max(1, t - equil_dist - 1):t - 1]))
        ))
    }
    # Return ####
    return(XT)
}

#' Simulate multiple scenarii for a species
#'
#' Wrapper function of the \code{sim_deter_forest} that loop on different
#' climatic values by loading the corresponding IPM and running the replicated
#' simulations.
#'
#' @param species Name of the species to run simulation on. Single char.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most usefull when creating random
#' population.
#' @param climatic Vector of climatic situations to run on. IPM must exist for
#' each one or else this climatic value will be skipped. int.
#' @inheritParams sim_deter_forest
#'
#' @param path Place to save the resulting file. Single Char.
#' @param parallel When multiple climatic values are given, it's possible to
#' run computation on multiple cores. FALSE by default.
#' @param save Save the simulation in an \code{.Rds} file located at
#' @param verbose Print message, used for debugs only. FALSE by default
#'
#' @details
#' When using parallel computations, the number of cores used will be the
#' minimum between the number of climate and the number of core of the computer
#' - 1.
#'
#' @return
#' Dataframe (tibble class) formated in long format. Column follow this order :
#' \itemize{
#'  \item{name - name of the variable. Contains state_eq, BA, N, BAeq and Neq}
#'  \item{t - time for the specified value. BA and N follow t from 1 to tlim,
#'  whereas equilibrium value are all a the time of equilibrium}
#'  \item{m - mesh element for the state_eq variable.}
#'  \item{delay - Delay setup used in the simulation}
#'  \item{clim - Climatic condition for IPM selection.}
#'  \item{harv - Harvest setup used in the simulation}
#'  \item{corr - Correction setup used in the simulation.}
#'  \item{n - Number of the simulation, if multiple IPM are provided per
#'  climate}
#'  \item{value - final value for the variable}
#' }
#'
#' @import here
#' @import checkmate
#' @importFrom dplyr relocate mutate
#' @importFrom parallel detectCores mclapply
#'
#' @export
run_sim_deter <- function(species,
                          init_pop, # IDEA  add init_pop as default function ? mais quelle pop par defaut ?
                          climatic = 1,
                          tlim = 3e3,
                          Harv = 0.006,
                          BAsup = 200,
                          correction = "none",
                          SurfEch = 0.03,
                          delay = 0,
                          path = here(),
                          parallel = FALSE,
                          save = FALSE,
                          verbose = FALSE) {

    # Idio Proof ####
    assertCharacter(species, len = 1)
    assertFunction(init_pop,args = c("mesh", "SurfEch"))
    assertNumeric(climatic)
    assertCharacter(path, len = 1)
    assertLogical(save, any.missing = FALSE, len = 1)
    if (save) {
        save_file <- here(path, "outputSim", species, "deter_sim.Rds")
        assertPathForOutput(save_file) # check we can save before we simulate !
    }
    assertLogical(parallel, any.missing = FALSE, len = 1)
    assertLogical(verbose, any.missing = FALSE, len = 1)
    if(verbose){
        message("proof done")
    }

    # This replace save_SimNonDem function
    nclim <- length(climatic)
    # Recursion ####
    if (nclim > 1) {
        if (verbose) {
            message("Recursion happenning now !")
        }
        if (parallel) {
            tmp <- c(nclim, detectCores() - 1)
            n.core <- tmp[which.min(tmp)]
            if (verbose) {
                message("going parallel on ", n.core, " cores.")
                system(paste("echo ''"))
            }
        } else {
            n.core <- 1
        }
        res <- mclapply(
            climatic, run_sim_deter,
            species = species, init_pop = init_pop, tlim = tlim,
            Harv = Harv, BAsup = BAsup, correction = correction,
            SurfEch = SurfEch, delay = delay, path = path,
            verbose = verbose, parallel = n.core > 1,
            mc.cores = n.core
        )

        res <- do.call(rbind, res)
        ## Early return ####
        if (save) {
            saveRDS(res, file = save_file)
        } else {
            return(res)
        }
    }
    # ..................................................................... ####
    # Single sim proof ####
    assertCount(climatic) # past recursive check
    if (verbose && parallel && interactive()) {
        # HACK message lost in parallel running
        system(paste("echo 'climatic", climatic, "start'"))
    }

    # Read IPM ####
    fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
    IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
    nIPM <- length(IPM)
    res <- vector(mode = "list", length = nIPM)

    # Run sims ####
    for (i in seq_along(IPM)) {
        if (verbose && i %% 20 == 0) {
            message(" ") ; message("run ", i, " on ", nIPM)
        }
        ## Init pop ####
        # ct <- drop(Buildct(IPM[[i]]$meshpts, SurfEch = SurfEch))
        # ini <- state_init(IPM[[i]]$meshpts)
        # ini <- as.numeric(ini / sum(ct * ini) * 112)
        ini <- do.call(init_pop, list(mesh = IPM[[i]]$meshpts,
                                      SurfEch = SurfEch))
        assertNumeric(ini, any.missing = FALSE) # TEMP move this into sim_deter as idiot proof !
        ## single sim ####
        sim_rest <- sim_deter_forest(
            Model = IPM[[i]], # delayed inside fct.
            Xini = delay(ini, delay), tlim = tlim, Harv = Harv, BAsup = BAsup,
            correction = correction, SurfEch = SurfEch, delay = delay,
            verbose = (verbose && i %% 20 == 0)
        )
        res[[i]] <- mutate(tree_format(sim_rest), n = i)
    }
    # IDEA txtplot here to resume computations ?
    # txtcurve(sin(pi*x),from=0,to=2)

    # formating ####
    res <- do.call(rbind, res)
    res <- mutate(
        res, clim = climatic, harv = Harv, corr = correction, delay = delay
    )
    res <- relocate(
        res, .data$name, .data$t, .data$m, .data$delay,
        .data$clim, .data$harv, .data$corr, .data$n, .data$value
    )

    # res %>% mutate_if(is.character,as.factor) %>% summary() # TEMP dev
    # Return ####
    if (verbose && parallel && interactive()) {
        # HACK message lost in parallel running
        system(paste("echo '", climatic, "done'"))
    }
    if (save) {
        saveRDS(res, file = save_file)
    } else {
        return(res)
    }
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
#' @import magrittr
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
