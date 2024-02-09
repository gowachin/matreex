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

# small example ####
# rm(list = ls())
base <- new.env()
base$X <- 1:5
base$nsp <- 0

foo <- function(){
    res <- new.env()
    res$X <- 1:5
    res$nsp <- 0
    return(res)
}

a <- foo()
b <- foo()

attach(base)
X
ls()

detach(base)
with(a,{
    X <- X +1
})
a$X <- a$X - 1
a$X
b$X

fuu <- function(x){
    attach(x)
}

rm(a, b, base, fuu, foo)

# function for mosaic ####

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

init_forest_env <- function(Mosaic,
                            index,
                            tlim,
                            run_disturb,
                            disturbance,
                            SurfEch = SurfEch
){

    # DEV
    # index <- "A"
    # tlim <- 10
    # SurfEch <- 0.03
    # EoDev

    res <- new.env()
    res$species <- Mosaic$forests[[index]]$species

    res$index <- index
    res$SurfEch <- SurfEch
    res$tlim <- tlim

    res$sp <- unname(Mosaic$forests[[index]]$info$species)
    nsp <- length(res$sp)

    res$meshs <- map(Mosaic$ipms, ~ .x$mesh)[res$sp]
    res$types <- map_chr(res$species, ~ .x$info["type"])
    res$stand_above_dth <- map2(res$meshs, res$species,  ~ .x > .y$harv_lim["dth"])
    res$BAsp <- map(Mosaic$ipms, ~ .x$BA)

    sim_X <- init_sim(nsp, tlim, res$meshs)
    # browser()
    res$sim_X <- as.data.frame(do.call("rbind", sim_X))
    res$sim_BAstand <- res$sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, res$sp)
    ))
    res$sim_BA <- rep(NA_real_, tlim)
    res$sim_BAnonSp <- rep(NA_real_, nsp)

    # Initiate the pops
    res$X <- map2(map(res$species, `[[`, "init_pop"),
              res$meshs,
              exec, SurfEch = SurfEch)
    res$Harv <- map(lengths(res$meshs), ~ rep(0, .x))
    res$ct <- map(res$meshs, Buildct, SurfEch = SurfEch)

    # disturbance
    # browser()
    if(run_disturb){
        res$disturbance <- dplyr::filter(disturbance, plot == res$index)
        t_disturb <- logical(tlim)
        t_disturb[res$disturbance$t] <- TRUE
        res$t_disturb <- t_disturb
    }


    # Harv cst
    res$harv_rules <- Mosaic$forests[[index]]$harv_rules
    res$favoured_species <- Mosaic$forests[[index]]$favoured_sp
    res$alpha <- Mosaic$forests[[index]]$harv_rule["alpha"]
    res$Pmax <- Mosaic$forests[[index]]$harv_rule["Pmax"]
    res$dBAmin <- Mosaic$forests[[index]]$harv_rule["dBAmin"]
    res$disturb <- FALSE
    res$disturb_surv <- TRUE

    return(res)
}

# x <- res

save_step_env <- function(x, t){

    # Dev
    # t <- 1
    # EoDev

    x$sim_BAsp[t, ] <- map2_dbl(x$X, x$ct, ~ .x %*% .y )
    x$standX <- map2(x$X, x$stand_above_dth, `*`)
    x$sim_BAstand[t, ] <- map2_dbl(x$standX, x$ct, `%*%`)
    x$sim_BA[t] <- sum(x$sim_BAsp[t,])
    # browser()
    #
    # microbenchmark::microbenchmark(
    #     simpl = x$sim_BA[t] - unlist(x$sim_BAsp[t, ,drop = TRUE]),
    #     map = map2_dbl( - x$sim_BAsp[t, ,drop = FALSE], x$sim_BA[t],  `+`)
    # )
    #
    # str(map2_dbl( - x$sim_BAsp[t, ,drop = FALSE], x$sim_BA[t],  `+`))
    # str( - unlist(x$sim_BAsp[t, ,drop = TRUE]) + x$sim_BA[t])

    x$sim_BAnonSp <- x$sim_BA[t] - unlist(x$sim_BAsp[t, ,drop = TRUE])
    # x$sim_BAnonSp <- map2_dbl( - x$sim_BAsp[t, ,drop = FALSE], x$sim_BA[t],  `+`)



    tmp <- imap(x$X, function(x, .y, ba, bast, harv, surf){
        c(x / surf , ba[[.y]], bast[[.y]], sum(x) / surf,
          harv[[.y]] / surf, sum(harv[[.y]]) / surf)
    },
    ba = x$sim_BAsp[t,, drop = FALSE],
    bast = x$sim_BAstand[t,,drop = FALSE],
    harv = x$Harv, surf = x$SurfEch )
    tmp <- do.call("c", tmp)

    x$sim_X[, t] <- tmp

    if (any(map2_lgl(x$sim_BA[t], x$BAsp, ~ ! between(.x, min(.y), max(.y -1))))) {
        waring("Border Basal Area reached for this simulation.")
        break()
    }

    return(x)
}

get_step_env <- function(x, Mosaic, t, climate, correction){

    # Dev
    # t <- 1
    # x <- landscape[[1]]
    # EoDev

    x$sim_ipm <- map(Mosaic$ipms[x$sp], ~ get_step_IPM(
        x = .x, BA = x$sim_BA[t], climate = climate, sim_corr = correction,
        IsSurv = x$disturb_surv
    ))

    return(x)
}


# sp_rec ####
#' sp recruit
#'
#' Get species recruitment function
#'
#' @param x Species to get the recruitment function from
#' @param climatic Climate vector is needed for mu_gr object to build the
#' corresponding recruitment function.
#' @param regional TRUE/FALSE if we want to use a regional BA for fecundity.
#'
#' @name mosaic_rec
#'
#' @export
mosaic_rec <- function(x, climatic, regional, FecnCom){
    UseMethod("mosaic_rec")
}

#' @method mosaic_rec ipm
#' @export
mosaic_rec.ipm <- function(x, climatic, regional = FALSE,
                           FecnCom = c("Fecundity", "Competition")){

    if(FecnCom == "Fecundity"){
        res <- exp_recFun(params = x$fit$fec$params_m, list_covs = climatic, regional = regional)
    } else {
        res <- mosa_recFun(params = x$fit$rec$params_m, list_covs = climatic)
    }
    return(res)
}

#' @method mosaic_rec species
#' @export
mosaic_rec.species <- function(x, climatic, regional = FALSE,
                               FecnCom = c("Fecundity", "Competition")){

    res <- mosaic_rec(x = x$IPM, climatic, regional = regional, FecnCom = FecnCom)
    return(res)
}

#' Export recruitment function from estimated parameters.
#'
#' Rebuild the function to use BASp and BAnonSp for a species.
#'
#' @param params Estimated parameters for the fit of the model.
#' @param list_covs Climatic covariates values.
#'
#' @importFrom purrr map
#' @importFrom rlang expr call2 env_unbind
#'
#' @details
#' Each function has an environment binded with params and list_covs.
#' I can't remove it and it may be usefull later after all.
#'
#' @return
#' Function with 4 parameters : BATOTSP, BATOTNonSP, mesh and SurfEch
#'
#' @examples
#' params <- c(intercept = -0.864, BATOTSP = -0.018, sgddb = 286.813,
#' wai = -0.057, wai2 = 0.288 )
#' list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)
#'
#' foo <- exp_recFun(params, list_covs)
#' foo
#' foo(1, 2, 1:5, 0.03)
#'
#' @noRd
mosa_recFun <- function(params, list_covs){

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    invar <- invar[! grepl("logBATOTSP", invar)]

    exp_invar <- map(invar, multi, df2)
    add_invar <- map(exp_invar,
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    final_res <- list(
        expr(final <- exp(res) * distrib ),
        expr(return(final))
    )
    calls <- c(exp_invar, add_invar, final_res)

    empty <- function(distrib, BATOTSP, BATOTNonSP){}

    body(empty)[[2]] <- expr(res <- 0)
    for(i in seq_along(calls)){
        body(empty)[[i + 2]] <- calls[[i]]
    }

    # this is messy and is to remove binded env.
    env_unbind(env = environment(empty), c("i", "calls", "final_res",
                                           "add_invar", "exp_invar",
                                           "invar", "df2", "SurfEch"),
               inherit = FALSE)
    # empty
    return(empty)
}


recrut_env <- function(landscape, t, sim_clim){

    # Dev
    # t <- 1]
    # plot <- landscape[[2]]
    # params <- plot$species[[1]]$IPM$fit$rec$params_m
    # list_covs <- sim_clim
    # EoDev
    all_recrues <- map(landscape, function(plot, t){
        fec <- map(plot$species, mosaic_rec.species, sim_clim, TRUE, "Fecundity")

        recrues <- imap(
            fec,
            function(x, .y, basp, mesh, SurfEch){
                if(basp[[.y]] == 0){ # if species is absent, no recruitment
                    return(mesh[[.y]] * 0)
                }
                exec(x, 0, basp[[.y]], 0, mesh[[.y]], SurfEch)
            }, basp = plot$sim_BAsp[t-1,,drop = FALSE],
            mesh = plot$meshs, SurfEch = plot$SurfEch)

        return(recrues)
    }, t = t) %>% flatten()

    # regroup all recrues
    nms_sp <- unique(names(all_recrues))
    rec <- vector("list", length(nms_sp))
    names(rec) <- nms_sp

    for(s in seq_along(all_recrues)){
        spi <- names(all_recrues)[s]
        if(is.null(rec[[spi]])){
            rec[[spi]] <- all_recrues[[s]]
        } else {
            rec[[spi]] <- rec[[spi]] + all_recrues[[s]]
        }
    }
    # divide total
    rec <- map(rec, ~ .x / length(landscape) )

    landscape <- map(landscape, function(plot, rec, t){
        comp <- map(plot$species, mosaic_rec.species, sim_clim, TRUE, "Competition")

        recrues <- imap(
            comp,
            function(x, .y, recr, basp, banonsp, recr){
                exec(x, recr[[.y]], basp[[.y]], banonsp[.y])
            }, recr = rec, basp = plot$sim_BAsp[t-1,,drop = FALSE],
            banonsp = plot$sim_BAnonSp)

        plot$X <- sapply(names(plot$X), function(n, x, y) x[[n]] + y[[n]],
                    plot$X, recrues, simplify = FALSE)


        return(plot)
    }, rec = rec, t = t)

    return(landscape)

}


growth_mortal_env <- function(x, t = 1, harvest, run_disturb, verbose){

    # Dev
    # t <- 2
    # x <- landscape[[1]]
    # harvest <- "default"
    # run_disturb <- FALSE
    # EoDev

    ## t size distrib ####
    x$X <- map2(x$X, x$sim_ipm, ~ drop( .y %*% .x ) )# Growth

    # browser()
    ## Disturbance ####
    if(run_disturb && x$t_disturb[t]){ # TODO edit the distribution table
        x$disturb <- TRUE

        if (verbose) {
            message(sprintf(
                "Plot %s time %i | Disturbance : %s I = %.2f",
                x$index, t, x$disturbance[x$disturbance$t == t, "type"],
                x$disturbance[x$disturbance$t == t, "intensity"]
            )
            )
        }

        qmd <- QMD(size = unlist(x$meshs), n = unlist(x$X))
        # TODO remove unborn size from X before computations

        # compute percentage of coniferous (relative share in number of stems)
        total_stem <- purrr::reduce(x$X, sum, .init = 0)
        sp_stem <- map_dbl(x$X, ~ sum(.x) / total_stem)
        perc_coni <- sum(sp_stem[names(x$types[x$types == "Coniferous"])])

        x$Disturb <- imap(
            # FIXME The correct distubr_fun() require the mesh from species...damn
            map(x$species, `[[`, "disturb_fun"),
            function(f, .y, X, sp, disturb, ...){
                exec(f, X[[.y]], sp[[.y]], disturb, ...)
            }, X = x$X, sp = x$species,
            disturb = x$disturbance[x$disturbance$t == t, ],
            qmd = qmd, perc_coni = perc_coni
        )

        x$X <- map2(x$X, x$Disturb, `-`)

    }


    ## Harvest ####
    if(!x$disturb && t %% x$harv_rules["freq"] == 0 &&
       harvest %in% c("Uneven", "Favoured_Uneven")){
        ### Uneven ####
        BAstandsp <- map2_dbl(x$X, x$species, getBAstand, SurfEch)
        BAstand <- sum(BAstandsp)
        BAcut <- getBAcutTarget(BAstand, targetBA, Pmax, dBAmin )

        sfav <- sum(x$favoured_sp)
        if( harvest == "Favoured_Uneven" && (sfav == 0 || sfav == length(x$favoured_sp))){
            print('!!!!!!!!!!!!!!!!!!!!!  WARNING  !!!!!!!!!!!!!!!!!!!!!')
            # warning("No species are favoured in the forest object, harvest mode 'Favoured_Uneven' is replaced with 'Uneven'")
            harvest <- "Uneven"
        }

        if(harvest == "Uneven"){
            pi <- BAstandsp / BAstand
            Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
            targetBAcut <- Hi * BAstandsp
        } else { # Favoured_Uneven
            p_fav <- sum(BAstandsp[x$favoured_sp])/BAstand
            if(p_fav > 0.5){
                Hi <- BAcut / BAstand
            }  else {
                pi <- ifelse(x$favoured_sp, p_fav, 1-p_fav)
                Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
            }
            targetBAcut <- Hi * BAstandsp
        }

        x$Harv <- imap(
            map(x$species, `[[`, "harvest_fun"),
            function(f, .y, X, sp, bacut, ct, ...){
                exec(f, X[[.y]], sp[[.y]],
                     targetBAcut = bacut[[.y]],
                     ct = ct[[.y]], ...)
            }, X = x$X, sp = x$species, bacut = targetBAcut,
            ct = x$ct, t = t
        )

        x$X <- map2(x$X, x$Harv, `-`)
    } else if(!x$disturb && harvest == "Even"){
        ### Even ####
        if(t %% final_harv == 0){
            x$Harv <- x$X
            x$X <- map2(map(x$species, `[[`, "init_pop"),
                        x$meshs,
                        exec, SurfEch = x$SurfEch)
        } else if(t %% x$harv_rule["freq"] == 0){
            x$Harv <- imap(
                map(x$species, `[[`, "harvest_fun"),
                function(f, .y, X, sp, tRDI, tKg, ct, ...){
                    exec(f, X[[.y]], sp[[.y]],
                         targetRDI = tRDI[[.y]],
                         targetKg = tKg[[.y]],
                         ct = ct[[.y]],
                         ...)
                }, X = x$X, sp = x$species, tRDI = targetRDI,
                tKg = targetKg, ct = x$ct, t = t, SurfEch = x$SurfEch
            )
            x$X <- map2(x$X, x$Harv, `-`)

        } else {
            x$Harv <- map(x$meshs, ~ rep(0, length(.x)))
        }
    } else if (!x$disturb && t %% x$harv_rules["freq"] == 0 && harvest == "default") {
        ### Nothing ####
        x$Harv <- imap(
            map(x$species, `[[`, "harvest_fun"),
            function(f, .y, X, sp, ct, ...){
                exec(f, X[[.y]], sp[[.y]], ct = ct[[.y]], ...)
            }, X = x$X, sp = x$species, ct = x$ct, t = t, SurfEch = x$SurfEch
        )

        x$X <- map2(x$X, x$Harv, `-`)
    } else if(x$disturb){
        x$Harv <- x$Disturb
        x$disturb <- FALSE
    } else {
        x$Harv <- map(x$meshs, ~ rep(0, length(.x)))
    }

    # Is there a disturbance ?
    if(run_disturb ){ # IDEA rewrite this ?
        if(x$t_disturb[t+1]){
            x$disturb_surv <- x$disturbance[x$disturbance$t == t+1, "IsSurv"]
        } else {
            x$disturb_surv <- TRUE
        }
    }
    # print(t) ; print(x$index);  print(x$disturb)

    ## Return ####
    return(x)

}


# Simulations ####
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


    # tlim = 100
    # harvest = "default"
    # climate = NULL
    # disturbance = NULL
    # correction = "none"
    # SurfEch = 0.03
    # verbose = TRUE

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

        assertList(disturbance, len = length(Mosaic$forests))

        disturbance <- dplyr::bind_rows(disturbance, .id = "plot")
        # TODO idiot proof disturbance
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
                                       run_disturb = run_disturb,
                                       disturbance = disturbance,
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
            .x, t = t, harvest = harvest, run_disturb, verbose
        ))
        ### Recruitment ####
        landscape <- recrut_env(landscape, t, sim_clim)
        ## Save BA ####
        landscape <- map(landscape, ~ save_step_env(.x, t = t))
        ## Get sim IPM ####
        landscape <-  map(landscape,
                          ~ get_step_env(.x, Mosaic, t= t, climate = sim_clim, correction))
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
    landsim <- map(landscape, function(plot){

        tmp <- new_deter_sim(as.matrix(plot$sim_X), mesh = plot$meshs)
        return(tree_format(tmp))
    })
    names(landsim) <- nms_plot

    if (verbose) {
        message("Simulation ended after time ", ifelse(is.na(the), t-1, the))
        tmp <- Sys.time() - start
        message("Time difference of ", format(unclass(tmp), digits = 3),
                " ", attr(tmp, "units"))
    }

    # Return ####
    final <- dplyr::bind_rows(landsim, .id = "plot")

    return(final)
}
