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

init_forest_env <- function(Mosaic,
                            index,
                            tlim,
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
    res$sim_X <- do.call("rbind", sim_X)
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
    x$sim_BAnonSp <- map2_dbl( - x$sim_BAsp[t, ,drop = FALSE], x$sim_BA[t],  `+`)

    tmp <- imap(x$X, function(x, .y, ba, bast, harv){
        c(x / SurfEch , ba[[.y]], bast[[.y]], sum(x) / SurfEch,
          harv[[.y]] / SurfEch, sum(harv[[.y]]) / SurfEch )
    },
    ba = x$sim_BAsp[t,, drop = FALSE],
    bast = x$sim_BAstand[t,,drop = FALSE],
    harv = x$Harv )
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


recrut_env <- function(landscape, t){

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
            function(x, .y, rec, basp, banonsp, recr){
                exec(x, recr[[.y]], basp[[.y]], banonsp[.y])
            }, recr = rec, basp = plot$sim_BAsp[t-1,,drop = FALSE],
            banonsp = plot$sim_BAnonSp)

        plot$X <- sapply(names(plot$X), function(n, x, y) x[[n]] + y[[n]],
                    plot$X, recrues, simplify = FALSE)


        return(plot)
    }, rec = rec, t = t)

    return(landscape)

}


growth_mortal_env <- function(x, t = 1, harvest, run_disturb){

    # Dev
    # t <- 2
    # x <- landscape[[1]]
    # harvest <- "default"
    # run_disturb <- FALSE
    # EoDev

    ## t size distrib ####
    x$X <- map2(x$X, x$sim_ipm, ~ drop( .y %*% .x ) )# Growth

    ## Disturbance ####
    if(FALSE && run_disturb && t_disturb[t]){ # TODO edit the distribution table
        x$disturb <- TRUE

        if (verbose) {
            message(sprintf(
                "Plot %s time %i | Disturbance : %s I = %.2f",
                x$index, t, disturbance[disturbance$t == t, "type"],
                disturbance[disturbance$t == t, "intensity"]
            )
            )
        }

        qmd <- QMD(size = unlist(x$meshs), n = unlist(x$X))
        # TODO remove unborn size from X before computations

        # compute percentage of coniferous (relative share in number of stems)
        total_stem <- purrr::reduce(x$X, sum, .init = 0)
        sp_stem <- map_dbl(x$X, ~ sum(.x) / total_stem)
        perc_coni <- sum(sp_stem[names(x$types[x$types == "Coniferous"])])

        Disturb <- imap(
            # FIXME The correct distubr_fun() require the mesh from species...damn
            map(x$species, `[[`, "disturb_fun"),
            function(f, .y, X, sp, disturb, ...){
                exec(f, X[[.y]], sp[[.y]], disturb, ...)
            }, X = x$X, sp = x$species,
            disturb = disturbance[disturbance$t == t, ],
            qmd = qmd, perc_coni = perc_coni
        )

        x$X <- map2(x$X, Disturb, `-`)

    }


    browser()
    print(x$disturb)
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
    } else if(disturb){
        x$Harv <- x$Disturb
        x$disturb <- FALSE
    } else {
        x$Harv <- map(x$meshs, ~ rep(0, length(.x)))
    }

    # Is there a disturbance ?
    if(run_disturb ){ # IDEA rewrite this ?
        if(t_disturb[t+1]){
            x$disturb_surv <- disturbance[disturbance$t == t+1, "IsSurv"]
        } else {
            x$disturb_surv <- TRUE
        }
    }

    print(t) ; print(x$index);  print(x$disturb)

    ## Return ####
    return(x)

}

