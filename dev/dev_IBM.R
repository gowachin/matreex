exp_allFun <- function(params, list_covs,
                       mod = c("rec", "gr", "sv")){

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(df2$K[! (df2$var1 %in% invar | df2$var2 %in% invar)])


    exp_invar <- map(invar, multi, df2)
    invar <- unlist(map(exp_invar, ~ attributes(.x)$var))
    invar <- unique(invar[!is.na(invar)])
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    # SurfEch <- NULL # hack rm the note in devtools::check() about unbinded
    final_res <- list(
        # expr(mesh <- length(mesh)),
        # expr(distrib <- c(rep(1/2, 2), numeric(mesh - 2)) ),
        # expr(final <- exp(res)), # * SurfEch / 0.03 * distrib ),
        # expr(return(final))
        expr(return(res))
    )
    calls <- c(exp_invar, add_invar, final_res)

    # if(regional){
    #     empty <- function(BATOTSP, BAFecSP, BATOTNonSP, mesh, SurfEch = 0.03){}
    # } else {
    #     empty <- function(BATOTSP, BATOTNonSP, mesh, SurfEch = 0.03){}
    # }

    # IDEA for latter
    empty <- function(){}
    # arguments <- vector("pairlist", length(invar) + 2)
    # names(arguments) <- c(invar, "mesh", "SurfEch")
    # arguments <- lapply(arguments, function(y) substitute(x, alist(x=)))
    # c("BATOTSP", if(TRUE){"BAFecSP"}else{ NULL}, "BATOTNonSP")
    # c("BATOTSP", if(FALSE){"BAFecSP"}else{ NULL}, "BATOTNonSP")
    arguments <- setNames(rep(alist(x=), length(invar) + 1), c(invar, "..."))
    # arguments$SurfEch <- 0.03
    formals(empty) <- arguments

    body(empty)[[2]] <- call2("<-", expr(intercept), inter)
    body(empty)[[3]] <- expr(res <- 0)
    for(i in seq_along(calls)){
        body(empty)[[i + 3]] <- calls[[i]]
    }

    # this is messy and is to remove binded env.
    env_unbind(env = environment(empty), c("i", "calls", "final_res",
                                           "add_invar", "exp_invar",
                                           "inter", "invar", "df2"),
               inherit = FALSE)
    # empty
    return(empty)
}

load_all()
species <- "Abies_alba"
data(list = paste0("fit_", species))
fit <- eval(parse(text=paste0("fit_", species)))
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == species,
                  select = -c(sp,N, PC1, PC2, SDM))

# source("dev/dev_IBM.R")
g <- exp_allFun(params = fit$gr$params_m, list_covs = climate)
r <- exp_allFun(params = fit$rec$params_m, list_covs = climate)
s <- exp_allFun(params = fit$sv$params_m, list_covs = climate)



# size <- c(300, 900)
# Harv <- 0.006
# l <- length(size) # length of pop
# # trouver une manière d'écrire ça mieux
# args <- c(list(BATOTcomp = 30, BATOTNonSP = 0, BATOTSP = 30, size = size),
#           as.list(climate))
#
# Grmean <- do.call(g, args = args)
# size + rlnorm(l, meanlog=Grmean, sdlog = fit$gr$sigma)
#
# Survmean <- do.call(s, args = args)
# P_sv <- (1 - fit$sv$family$linkinv(Survmean)) * (1 - Harv)
# Surv <- rbinom(l, 1, P_sv)
#
# Surv[size > as.numeric(fit$info["max_dbh"])] <- 0
#
# Recmean <- do.call(r, args = args)
# # Nrec <- rnbinom(1, mu=exp(Recmean), size=fit$rec$sigma)
# Nrec <- rnbinom(1, mu=exp(Recmean), size=0.8)


#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||#

load_all()
data("fit_Picea_abies")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
climate
Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
Picea_sp <- species(IPM = Picea_ipm, init_pop = def_initBA(30))
Forest <- forest(species = list(Picea = Picea_sp))

X2Pop <- function(X, mesh){

    # TEMP dev
    # X <- X[[1]]
    # mesh <- meshs[[1]]
    # TEMP dev

    lag <- sum(mesh == 0)
    mesh[1:lag] <- -c(lag:1)
    Pop <- rpois(length(X), X)
    res <- rep(mesh, times = Pop)
    # rep is here super simple and efficient, maybe take size in unif in mesh when multiple indiv
    # BA <- sum(pi*(res[res>0]/2*1e-3)^2 / 0.03)
    return(res)
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
#' @param climate Optional, climate matrix if climate variation along time is
#' needed. Climate variation rely on species created with mu_gr class objects.
#' This matrix require as many rows as time steps until equil_time.
#' If the climate does not variate, a single row can given and will be reused.
#' @param disturbance `r lifecycle::badge("experimental")` parameter.
#' @param SurfEch Value of plot size surface in ha
#'
#' @param verbose Print message. FALSE by default
#'
#' @details
#' Basic simulations input are illustrated in the main vignette.
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
#'   \item{mesh}{Mesh class number, from 1 to n class. NA for variable for individual}
#'   \item{size}{Size corresponding to the mesh class. NA for variable for individual}
#'   \item{equil}{Logical if this time step is the equilibrium or last step of
#'   simulation}
#'   \item{value}{Numeric values of the variables.}
#'  }
#'
#' The variables are :
#' \describe{
#'  \item{N}{Sum of density for sampled area (SurfEch).}
#'  \item{BAsp}{Basal area of the population per sampled area (SurfEch) and species}
#'  \item{BAstand}{Basal area of the population per sampled area (SurfEch)
#'  and species when excluding size class below dth. See Harvesting vignette.}
#'  \item{H}{Sum of harvested density per sampled area (SurfEch).}
#' }
#'
#' @import Matrix
#' @import checkmate
#' @import purrr
#' @importFrom dplyr between bind_cols
#'
#' @name sim_deter_forest
#' @export
sim_indiv_forest  <- function(Forest,
                              tlim = 3e3,
                              climate = NULL,
                              disturbance = NULL,
                              harvest = 0.006,
                              SurfEch = 0.03,
                              verbose = FALSE) {
    UseMethod("sim_indiv_forest")
}

#' @method sim_deter_forest forest
#' @export
sim_indiv_forest.forest  <- function(Forest,
                                     tlim = 3e3,
                                     climate = NULL,
                                     disturbance = NULL,
                                     harvest = 0.006,
                                     SurfEch = 0.03,
                                     verbose = FALSE) {


    # TEMP dev
    # tlim = 500
    # SurfEch = 0.03
    # # climate = NULL
    # disturbance = NULL
    # verbose = FALSE
    # harvest = 0.006
    # TEMP dev

    # Idiot Proof ####
    # validate_forest(Forest) # TEMP dev
    assertCount(tlim)
    IPM_cl <- map_chr(Forest$species, ~ class(.x$IPM))
    if(all(IPM_cl == "ipm") && !is.null(climate)) {
        # no climate needed
        warning(paste0("Because all species are fully integrated on a climate, ",
                       "providing one now is unnecessary"))
        clim_i <- which(IPM_cl == "ipm")[[1]]
        climate <- t(Forest$species[[clim_i]]$IPM$climatic)
        climate <-  as.matrix(bind_cols(climate, t = 1:tlim))
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
        t_disturb <- logical(tlim)
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
    assertNumber(SurfEch, lower = 0)
    assertLogical(verbose, any.missing = FALSE, len = 1)

    start <- Sys.time()

    # Initialisation ####
    init_sim <- function(nsp, tlim, mesh){ # TODO : set function outside of here
        res <- vector("list", nsp)
        res <- map(names(mesh), function(x) {
            tmp <- matrix(
                data = NA_real_, ncol = tlim + 1,
                nrow = 4
            )
            colnames(tmp) <- c(paste0("t", 1:(tlim+1))) #, "sp")
            rownames(tmp) <- paste0(x, c(".BAsp", ".BAstand", ".N", ".H"))
            return(tmp)
        })
        names(res) <- names(mesh)
        return(res)
    }
    nsp <- length(Forest$species)
    disturb_surv <- TRUE

    ## Modify IPM ####
    regional <- inherits(Forest, "reg_forest")

    meshs <- map(Forest$species, ~ .x$IPM$mesh)
    stand_above_dth <- map_dbl(Forest$species, ~ .x$harv_lim["dth"])
    lag_rec <- map(Forest$species, ~ as.numeric(.x$IPM$info["delay"]))
    gr_sig <- map_dbl(Forest$species, ~ .x$IPM$fit$gr$sigma)
    sv_sig <- map_dbl(Forest$species, ~ .x$IPM$fit$gr$sigma) # TODO change for rec
    maxdbh <- map_dbl(Forest$species, ~ as.numeric(.x$IPM$fit$info["max_dbh"]))
    linkinv <- map(Forest$species, ~.x$IPM$fit$sv$family$linkinv)

    ## Create output ####
    sim_X <- init_sim(nsp, tlim, meshs)
    sim_X <- do.call("rbind", sim_X)
    sim_BAstand <- sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, names(Forest$species))
    ))
    sim_BA <- rep(NA_real_, tlim)

    ## Initiate pop ####
    X <- map2(map(Forest$species, `[[`, "init_pop"),
              meshs,
              exec, SurfEch = SurfEch)
    X <- map2(X, meshs, X2Pop)
    Harv <- map_dbl(meshs, ~ 0)
    start_clim <- climate[1, , drop = TRUE]
    # ct <- map(meshs, Buildct, SurfEch = SurfEch)

    bas <- c(list(BATOTcomp = NA, BATOTNonSP = NA, BATOTSP = NA),
             as.list(start_clim))

    # save first pop
    sim_BAsp[1, ] <- bas$BATOTSP<- map_dbl(X, ~ sum(pi*(.x[.x>0]/2*1e-3)^2 / SurfEch))
    sim_BAstand[1, ] <- map2_dbl(X, stand_above_dth,
                                 ~ sum(pi*(.x[.x>.y]/2*1e-3)^2 / SurfEch))
    sim_BA[1] <- bas$BATOTcomp <- sum(sim_BAsp[1,])
    bas$BATOTNonSP <- map2_dbl( - sim_BAsp[1, ,drop = FALSE], sim_BA[1],  `+`)

    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(ba[[.y]], bast[[.y]], length(x), harv)
    },
    ba = sim_BAsp[1,, drop = FALSE],
    bast = sim_BAstand[1,,drop = FALSE],
    harv = Harv )

    tmp <- do.call("c", tmp)
    sim_X[, 1] <- tmp

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

    # While tlim & eq ####
    t <- 2
    the <- NA_real_ # real time of simulation ending in case of outbound BA
    # Harv cst
    disturb <- FALSE

    # vital functions
    g_fun <- map(Forest$species, ~ exp_allFun(params =.x$IPM$fit$gr$params_m,
                                              list_covs = start_clim))
    r_fun <- map(Forest$species, ~ exp_allFun(params =.x$IPM$fit$rec$params_m,
                                              list_covs = start_clim))
    s_fun <- map(Forest$species, ~ exp_allFun(params =.x$IPM$fit$sv$params_m,
                                              list_covs = start_clim))

    if (verbose) {
        message("Starting while loop. Maximum t = ", tlim)
    }
    while (t < tlim ) {

        #' x is population, one size per individual,
        #' n is size in mm,
        #' 0 code for dead,
        #' -n code for lag until real size

        # # population that work well..
        # set.seed(974)
        # X <- map2(map(Forest$species, `[[`, "init_pop"),
        #           meshs,
        #           exec, SurfEch = SurfEch)
        # X <- map2(X, meshs, X2Pop)
        # X

        # print(X)
        # Growth
        X <- imap(
            g_fun,
            function(.g, .y, x, sigma, bas, ...){
                # Debug growth ||||||||||||||||||||||||||||||#
                # .g <- g_fun[[1]]
                # x <- X[[1]]
                # sigma = Forest$species[[1]]$IPM$fit$gr$sigma
                # bas = bas
                #||||||||||||||||||||||||||||||||||||||||||||#
                #
                # grow trees
                x <- x[[.y]]
                sigma <- sigma[[.y]]


                Grmean <- do.call(.g, args = c(list(size = x[x > 0]),
                                               as.list(bas)))
                x[x > 0] <- x[x > 0] + rlnorm(sum(x>0), meanlog=Grmean,
                                              sdlog = sigma)
                # grow lag
                x[x == -1] <- 90 # TODO == 0 or == -1 ?? -1 allow clean survival function
                x[x < 0] <- x[x < 0] + 1
                return(x)
            },
            x = X,
            sigma = gr_sig,
            bas = bas
        )

        # Survival
        X <- imap(
            s_fun,
            function(.s, .y, x, link, maxdbh, bas, harv, ...){
                # Debug growth ||||||||||||||||||||||||||||||#
                # .y <- "Picea_abies"
                # .s <- s_fun[[.y]]
                # x <- X
                # bas = bas
                # link = linkinv
                # maxdbh = maxdbh
                # harv = 0.006
                #||||||||||||||||||||||||||||||||||||||||||||#
                #
                # grow trees
                x <- x[[.y]]
                link <- link[[.y]]

                Survmean <- do.call(.s, args = c(list(size = x[x > 0]),
                                                 as.list(bas)))
                P_sv <- (1 - link(Survmean)) * (1 - harv)
                # print(P_sv)
                Surv <- rbinom(sum(x>0), 1, P_sv)
                # remove individual above maxdbh
                Surv[x[x > 0] > maxdbh[[.y]]] <- 0
                x[x > 0] <- x[x > 0] * Surv
                Harv[[.y]] <<- sum(x == 0) # Is this really a good idea ?
                x <- x[x != 0] # 0 code for dead

                return(x)
            },
            x = X,
            link = linkinv,
            maxdbh = maxdbh,
            harv = harvest,
            bas = bas
        )

        # Recruitment
        X <- imap(
            r_fun,
            function(.r, .y, x, sigma, bas, lag, ...){
                # Debug growth ||||||||||||||||||||||||||||||#
                # .r <- r_fun[[1]]
                # x <- X[[1]]
                # # sigma = Forest$species[[1]]$IPM$fit$gr$sigma
                # sigma = 0.8
                # bas = bas
                # lag = as.numeric(Forest$species[[1]]$IPM$info["delay"])
                #||||||||||||||||||||||||||||||||||||||||||||#
                x <- x[[.y]]
                bas$BATOTSP <- bas$BATOTSP[[.y]]
                bas$BATOTNonSP <- bas$BATOTNonSP[[.y]]

                # grow trees
                Recmean <- do.call(.r, args = c(list(size = x[x > 0]),
                                               as.list(bas)))
                Nrec <- rnbinom(1, mu=exp(Recmean), size=sigma)
                # add lag
                x <- c(rep(-lag[[.y]], times = Nrec), x)
                return(x)
            },
            x = X,
            # sigma = Forest$species[[.y]]$IPM$fit$gr$sigma,
            sigma = 0.8, # TODO edit the fit dataset from Julien !
            bas = bas,
            lag = lag_rec
        )

        ## Save BA ####
        # compute new BA for selecting the right IPM and save values
        sim_BAsp[t, ] <- bas$BATOTSP <- map_dbl(X, ~ sum(pi*(.x[.x>0]/2*1e-3)^2 / SurfEch))
        sim_BAstand[t, ] <- map2_dbl(X, stand_above_dth,
                                     ~ sum(pi*(.x[.x>.y]/2*1e-3)^2 / SurfEch))
        sim_BA[t] <- bas$BATOTcomp <- sum(sim_BAsp[t,])
        bas$BATOTNonSP <- map2_dbl( - sim_BAsp[t, ,drop = FALSE], sim_BA[t],  `+`)

        # update climate
        bas[colnames(climate)] <- as.list(climate[t,])

        # Update X and extract values per ha
        tmp <- imap(X, function(x, .y, ba, bast, harv){
            c(ba[[.y]], bast[[.y]], length(x), harv)
        },
        ba = sim_BAsp[t,, drop = FALSE],
        bast = sim_BAstand[t,,drop = FALSE],
        harv = Harv )

        tmp <- do.call("c", tmp)
        sim_X[, t] <- tmp

        # ## Stop loop if no populations left ####
        if (all(lengths(X) == 0)) {
            the <- t
            break()
        }

        ## Loop Verbose ####
        if (t %% 500 == 0 && verbose) {
            message(sprintf(
                "time %i | BA diff : %.2f",
                t, diff(range(sim_BA[max(1, t - tlim):t]))
            ))
        }
        t <- t + 1
    }

    # Format output ####
    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(ba[[.y]], bast[[.y]], length(x), harv)
    },
    ba = sim_BAsp[t-1,, drop = FALSE],
    bast = sim_BAstand[t-1,,drop = FALSE],
    harv = Harv )

    tmp <- do.call("c", tmp)
    sim_X[, tlim +1] <- tmp

    colnames(sim_X)[tlim + 1] <- paste0("eqt", t-1)


    if (verbose) {
        message("Simulation ended after time ", ifelse(is.na(the), t-1, the))
        tmp <- Sys.time() - start
        message("Time difference of ", format(unclass(tmp), digits = 3),
                " ", attr(tmp, "units"))
    }
    sim_X <- new_indiv_sim(sim_X, mesh = meshs)
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
new_indiv_sim <- function(x = matrix(), mesh = NULL){
    assertMatrix(x)
    structure(x, class = c("indiv_sim", "matrix"),
              mesh = mesh)

}

#' @rdname tree_format
#' @importFrom dplyr relocate
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map flatten_dbl
#' @keywords internal
#' @export
tree_format.indiv_sim <- function(x){

    mesh <- rep(NA_real_, nrow(x))

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

sim1 <- sim_indiv_forest(Forest, verbose = TRUE, tlim = 5000)

library(ggplot2)
sim1  %>%
    # dplyr::filter(var == "BAsp", ! equil) %>%
    ggplot(aes(x = time, y = value)) +
    geom_line(linewidth = .4) + ylab("BA") +
    facet_wrap(~ var, scales = "free_y")


sim1 %>%
    dplyr::filter(!is.na(value))
