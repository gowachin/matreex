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



size <- c(300, 900)
Harv <- 0.006
l <- length(size) # length of pop
# trouver une manière d'écrire ça mieux
args <- c(list(BATOTcomp = 30, BATOTNonSP = 0, BATOTSP = 30, size = size),
          as.list(climate))

Grmean <- do.call(g, args = args)
size + rlnorm(l, meanlog=Grmean, sdlog = fit$gr$sigma)

Survmean <- do.call(s, args = args)
P_sv <- (1 - fit$sv$family$linkinv(Survmean)) * (1 - Harv)
Surv <- rbinom(l, 1, P_sv)

Surv[size > as.numeric(fit$info["max_dbh"])] <- 0

Recmean <- do.call(r, args = args)
# Nrec <- rnbinom(1, mu=exp(Recmean), size=fit$rec$sigma)
Nrec <- rnbinom(1, mu=exp(Recmean), size=0.8)




ind <- list(
    sp = "Abies_alba",
    size = 90,
    tinit = 0,
    final = "alive",
    hist = NULL
)

ind
class(ind) <- paste0("ind_", ind$sp)
ind

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

sim_indiv_forest.forest  <- function(Forest,
                                     tlim = 3e3,
                                     climate = NULL,
                                     disturbance = NULL,
                                     SurfEch = 0.03,
                                     verbose = FALSE) {


    # TEMP dev
    tlim = 500
    SurfEch = 0.03
    # climate = NULL
    disturbance = NULL
    verbose = FALSE
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
    types <- map_chr(Forest$species, ~ .x$info["type"])
    stand_above_dth <- map_dbl(Forest$species, ~ .x$harv_lim["dth"])
    delay <- map(Forest$species, ~ as.numeric(.x$IPM$info["delay"]))

    ## Create output ####
    sim_X <- init_sim(nsp, tlim, meshs)
    sim_X <- do.call("rbind", sim_X)
    sim_BAstand <- sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, names(Forest$species))
    ))
    sim_BA <- rep(NA_real_, tlim)
    sim_BAnonSp <- rep(NA_real_, nsp)

    ## Initiate pop ####
    X <- map2(map(Forest$species, `[[`, "init_pop"),
              meshs,
              exec, SurfEch = SurfEch)
    X <- map2(X, meshs, X2Pop)
    Harv <- map_dbl(meshs, ~ 0)
    # ct <- map(meshs, Buildct, SurfEch = SurfEch)

    # save first pop
    sim_BAsp[1, ] <- map_dbl(X, ~ sum(pi*(.x[.x>0]/2*1e-3)^2 / SurfEch))
    sim_BAstand[1, ] <- map2_dbl(X, stand_above_dth,
                                 ~ sum(pi*(.x[.x>.y]/2*1e-3)^2 / SurfEch))
    sim_BA[1] <- sum(sim_BAsp[1,])
    sim_BAnonSp <- map2_dbl( - sim_BAsp[1, ,drop = FALSE], sim_BA[1],  `+`)

    tmp <- imap(X, function(x, .y, ba, bast, harv){
        c(ba[[.y]], bast[[.y]], length(x), harv)
    },
    ba = sim_BAsp[1,, drop = FALSE],
    bast = sim_BAstand[1,,drop = FALSE],
    harv = Harv )

    tmp <- do.call("c", tmp)
    sim_X[, 1] <- tmp

    # if (any(map2_lgl(sim_BA[1], BAsp, ~ ! between(.x, min(.y), max(.y -1))))) {
    #     stop(paste(
    #         "Border Basal Area reached for this simulation.",
    #         "This maximum is reached before iteration, check init_pop functions"
    #     ))
    # }

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
    alpha <- Forest$harv_rule["alpha"]
    Pmax <- Forest$harv_rule["Pmax"]
    dBAmin <- Forest$harv_rule["dBAmin"]
    disturb <- FALSE

    # vital functions
    start_clim <- climate[1, , drop = TRUE]
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

        # Growth
        g_fun
        # Survival
        s_fun
        # Recruitment
        r_fun


        size <- c(300, 900)
        Harv <- 0.006
        l <- length(size) # length of pop
        # trouver une manière d'écrire ça mieux
        args <- c(list(BATOTcomp = 30, BATOTNonSP = 0, BATOTSP = 30, size = size),
                  as.list(climate))

        Grmean <- do.call(g, args = args)
        size + rlnorm(l, meanlog=Grmean, sdlog = fit$gr$sigma)

        Survmean <- do.call(s, args = args)
        P_sv <- (1 - fit$sv$family$linkinv(Survmean)) * (1 - Harv)
        Surv <- rbinom(l, 1, P_sv)

        Surv[size > as.numeric(fit$info["max_dbh"])] <- 0

        Recmean <- do.call(r, args = args)
        # Nrec <- rnbinom(1, mu=exp(Recmean), size=fit$rec$sigma)
        Nrec <- rnbinom(1, mu=exp(Recmean), size=0.8)


        # ## t size distrib ####
        # X <- map2(X, sim_ipm, ~ drop( .y %*% .x ) )# Growth
        #
        # ## Disturbance ####
        # if(run_disturb && t_disturb[t]){
        #     disturb <- TRUE
        #
        #     if (verbose) {
        #         message(sprintf(
        #             "time %i | Disturbance : %s I = %.2f",
        #             t, disturbance[disturbance$t == t, "type"],
        #             disturbance[disturbance$t == t, "intensity"]
        #         )
        #         )
        #     }
        #
        #     qmd <- QMD(size = unlist(meshs), n = unlist(X))
        #     # TODO remove unborn size from X before computations
        #
        #     total_stem <- purrr::reduce(X, sum, .init = 0)
        #     sp_stem <- map_dbl(X, ~ sum(.x) / total_stem)
        #     perc_coni <- sum(sp_stem[names(types[types == "Coniferous"])])
        #
        #     Disturb <- imap(
        #         map(Forest$species, `[[`, "disturb_fun"),
        #         function(f, .y, X, sp, disturb, ...){
        #             exec(f, X[[.y]], sp[[.y]], disturb, ...)
        #         }, X = X, sp = Forest$species,
        #         disturb = disturbance[disturbance$t == t, ],
        #         qmd = qmd, perc_coni = perc_coni
        #     )
        #
        #     X <- map2(X, Disturb, `-`)
        #
        # }
        #
        # ## Harvest ####
        # if(!disturb && t %% Forest$harv_rule["freq"] == 0 &&
        #    harvest %in% c("Uneven", "Favoured_Uneven")){
        #     ### Uneven ####
        #     BAstandsp <- map2_dbl(X, Forest$species, getBAstand, SurfEch)
        #     BAstand <- sum(BAstandsp)
        #     BAcut <- getBAcutTarget(BAstand, targetBA, Pmax, dBAmin )
        #
        #     sfav <- sum(Forest$favoured_sp)
        #     if( harvest == "Favoured_Uneven" && (sfav == 0 || sfav == length(Forest$favoured_sp))){
        #         print('!!!!!!!!!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!!!!!!')
        #         # warning("No species are favoured in the forest object, harvest mode 'Favoured_Uneven' is replaced with 'Uneven'")
        #         harvest <- "Uneven"
        #     }
        #
        #     if(harvest == "Uneven"){
        #         pi <- BAstandsp / BAstand
        #         Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
        #         targetBAcut <- Hi * BAstandsp
        #     } else { # Favoured_Uneven
        #         p_fav <- sum(BAstandsp[Forest$favoured_sp])/BAstand
        #
        #         if(p_fav > 0.5){
        #             Hi <- BAcut / BAstand
        #         }  else {
        #             pi <- ifelse(Forest$favoured_sp, p_fav, 1-p_fav)
        #             Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
        #         }
        #         targetBAcut <- Hi * BAstandsp
        #     }
        #
        #     Harv <- imap(
        #         map(Forest$species, `[[`, "harvest_fun"),
        #         function(f, .y, X, sp, bacut, ct, ...){
        #             exec(f, X[[.y]], sp[[.y]],
        #                  targetBAcut = bacut[[.y]],
        #                  ct = ct[[.y]], ...)
        #         }, X = X, sp = Forest$species, bacut = targetBAcut,
        #         ct = ct, t = t
        #     )
        #
        #     X <- map2(X, Harv, `-`)
        # } else if(!disturb && harvest == "Even"){
        #     ### Even ####
        #     if(t %% final_harv == 0){
        #         Harv <- X
        #         X <- map2(map(Forest$species, `[[`, "init_pop"),
        #                   meshs,
        #                   exec, SurfEch = SurfEch)
        #     } else if(t %% Forest$harv_rule["freq"] == 0){
        #
        #         rdi_sp <- map2_dbl(X, Forest$species, RDI_sp,
        #                            SurfEch = SurfEch)
        #         rdi <- sum(rdi_sp)
        #
        #         if (rdi < targetRDI) {
        #             Harv <- map(meshs, ~ rep(0, length(.x)))
        #         } else {
        #             Pcut <- getPcutEven(x = X, sp = Forest$species,
        #                                 meshs = meshs,
        #                                 targetRDI = targetRDI,
        #                                 targetKg = targetKg,
        #                                 SurfEch = SurfEch
        #             )
        #
        #             Harv <- map2(X, Pcut, ~ .x * .y)
        #         }
        #         X <- map2(X, Harv, `-`)
        #
        #     } else {
        #         Harv <- map(meshs, ~ rep(0, length(.x)))
        #     }
        # } else if (!disturb && t %% Forest$harv_rule["freq"] == 0 && harvest == "default") {
        #     ### Nothing ####
        #     Harv <- imap(
        #         map(Forest$species, `[[`, "harvest_fun"),
        #         function(f, .y, X, sp, ct, ...){
        #             exec(f, X[[.y]], sp[[.y]], ct = ct[[.y]], ...)
        #         }, X = X, sp = Forest$species, ct = ct, t = t, SurfEch = SurfEch
        #     )
        #
        #     X <- map2(X, Harv, `-`)
        # } else if(disturb){
        #     Harv <- Disturb
        #     disturb <- FALSE
        # } else {
        #     Harv <- map(meshs, ~ rep(0, length(.x)))
        # }
        #
        # ### Recruitment ####
        # sim_clim <- climate[t, , drop = TRUE]
        # rec <- map(Forest$species, sp_rec.species, sim_clim)
        #
        # recrues <- imap(
        #     rec,
        #     function(x, .y, basp, banonsp, mesh, SurfEch, mig){
        #         if(basp[[.y]] == 0){ # if species is absent, no recruitment
        #             return(mesh[[.y]] * 0)
        #         }
        #         exec(x, basp[[.y]], banonsp[.y], mesh[[.y]], SurfEch) * (1 - mig[[.y]])
        #     }, basp = sim_BAsp[t-1,,drop = FALSE], banonsp = sim_BAnonSp,
        #     mesh = meshs, SurfEch = SurfEch, mig = migrate )
        #
        # if(regional){
        #     rec_reg <- map(Forest$species, sp_rec.species, sim_clim, TRUE)
        #
        #     reg_recrues <- imap(
        #         rec_reg,
        #         function(x, .y, basp, bareg, banonsp, mesh, SurfEch, mig){
        #             exec(x, basp[[.y]], bareg[[.y]], banonsp[.y], mesh[[.y]], SurfEch) * mig[[.y]]
        #         }, basp = sim_BAsp[t-1,,drop = FALSE], bareg = reg_ba, banonsp = reg_banonsp,
        #         mesh = meshs, SurfEch = SurfEch, migrate )
        # } else {
        #     reg_recrues <- map(recrues, ~ .x * 0)
        # }
        #
        # # X <- map2(X, recrues, `+`) # gain time
        # X <- sapply(names(X), function(n, x, y, z) x[[n]] + y[[n]] + z[[n]],
        #             X, recrues, reg_recrues, simplify = FALSE)
        #
        ## Save BA ####
        # compute new BA for selecting the right IPM and save values
        sim_BAsp[t, ] <- map_dbl(X, ~ sum(pi*(.x[.x>0]/2*1e-3)^2 / SurfEch))
        sim_BAstand[t, ] <- map2_dbl(X, stand_above_dth,
                                     ~ sum(pi*(.x[.x>.y]/2*1e-3)^2 / SurfEch))
        sim_BA[t] <- sum(sim_BAsp[t,])
        sim_BAnonSp <- map2_dbl( - sim_BAsp[t, ,drop = FALSE], sim_BA[t],  `+`)

        # Update X and extract values per ha
        if (t <= tlim) {
            tmp <- imap(X, function(x, .y, ba, bast, harv){
                c(ba[[.y]], bast[[.y]], length(x), harv)
            },
            ba = sim_BAsp[t,, drop = FALSE],
            bast = sim_BAstand[t,,drop = FALSE],
            harv = Harv )

            tmp <- do.call("c", tmp)
            sim_X[, t] <- tmp
        }


        # ## Stop loop if BA larger than LIPM largest BA ####
        # if (any(map2_lgl(sim_BA[t], BAsp, ~ ! between(.x, min(.y), max(.y))))) {
        #     warning("Maximum Basal Area reached for this simulation.")
        #     # TODO  say which species reached BA limit !
        #     the <- t
        #     break()
        # }
        #
        # ## Get sim IPM ####
        # # Is there a disturbance ?
        # if(run_disturb && t < equil_time){ # IDEA rewrite this ?
        #     if(t_disturb[t+1]){
        #         disturb_surv <- disturbance[disturbance$t == t+1, "IsSurv"]
        #     } else {
        #         disturb_surv <- TRUE
        #     }
        # }
        #
        # sim_ipm <- map(Forest$species, ~ get_step_IPM(
        #     x = .x$IPM, BA = sim_BA[t], climate = sim_clim, sim_corr = correction,
        #     IsSurv = disturb_surv
        # ))

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



