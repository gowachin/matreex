# small example ####
rm(list = ls())
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
X
with(a,{
    X <- X +1
})
a$X <- a$X - 1
a$X
b$X

fuu <- function(x){
    attach(x)
}

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

    res$species <- Mosaic$forests[[index]]$species
    res <- new.env()

    res$index <- index
    res$SurfEch <- SurfEch
    res$tlim <- tlim

    res$sp <- unname(Mosaic$forests[[index]]$info$species)
    nsp <- length(res$sp)

    res$meshs <- map(Mosaic$ipms, ~ .x$mesh)[res$sp]
    res$types <- map_chr(species, ~ .x$info["type"])
    res$stand_above_dth <- map2(res$meshs, species,  ~ .x > .y$harv_lim["dth"])
    res$BAsp <- map(Mosaic$ipms, ~ .x$BA)

    sim_X <- init_sim(nsp, tlim, res$meshs)
    res$sim_X <- do.call("rbind", sim_X)
    res$sim_BAstand <- res$sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, res$sp)
    ))
    res$sim_BA <- rep(NA_real_, tlim)
    res$sim_BAnonSp <- rep(NA_real_, nsp)

    # Initiate the pops
    res$X <- map2(map(species, `[[`, "init_pop"),
              res$meshs,
              exec, SurfEch = SurfEch)
    res$Harv <- map(lengths(res$meshs), ~ rep(0, .x))
    res$ct <- map(res$meshs, Buildct, SurfEch = SurfEch)

    # Harv cst
    res$alpha <- Forest$harv_rule["alpha"]
    res$Pmax <- Forest$harv_rule["Pmax"]
    res$dBAmin <- Forest$harv_rule["dBAmin"]
    res$disturb <- FALSE
    res$disturb_surv <- TRUE

    return(res)
}

x <- res

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

growth_mortal_env <- function(x, t = 1, harvest, run_disturb){

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
            }, x$X = X, sp = x$species,
            disturb = disturbance[disturbance$t == t, ],
            qmd = qmd, perc_coni = perc_coni
        )

        x$X <- map2(x$X, Disturb, `-`)

    }

    ## Harvest ####
    if(!x$disturb && t %% Forest$harv_rule["freq"] == 0 &&
       harvest %in% c("Uneven", "Favoured_Uneven")){
        ### Uneven ####
        BAstandsp <- map2_dbl(x$X, Forest$species, getBAstand, SurfEch)
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
            if(p_fav > 0.5){
                Hi <- BAcut / BAstand
            }  else {
                pi <- ifelse(Forest$favoured_sp, p_fav, 1-p_fav)
                Hi <- BAcut / BAstand * ((pi ^ (alpha - 1)) / sum(pi ^ alpha))
            }
            targetBAcut <- Hi * BAstandsp
        }

        x$Harv <- imap(
            map(Forest$species, `[[`, "harvest_fun"),
            function(f, .y, X, sp, bacut, ct, ...){
                exec(f, X[[.y]], sp[[.y]],
                     targetBAcut = bacut[[.y]],
                     ct = ct[[.y]], ...)
            }, X = x$X, sp = x$species, bacut = targetBAcut,
            ct = x$ct, t = t
        )

        x$X <- map2(x$X, x$Harv, `-`)
    } else if(!disturb && harvest == "Even"){
        ### Even ####
        if(t %% final_harv == 0){
            x$Harv <- x$X
            x$X <- map2(map(x$species, `[[`, "init_pop"),
                      x$meshs,
                      exec, SurfEch = x$SurfEch)
        } else if(t %% Forest$harv_rule["freq"] == 0){
            x$Harv <- imap(
                map(Forest$species, `[[`, "harvest_fun"),
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
    } else if (!disturb && t %% Forest$harv_rule["freq"] == 0 && harvest == "default") {
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

    ## Return ####
    return(x)

}

recrut_env <- function(x){


    x$rec_sp <- map(x$species, ~ .x$prout)

    rec <- map(x$species, sp_rec.species, sim_clim)

    res <- NULL

    return(res)
}
