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
                            equil_time,
                            SurfEch = SurfEch
){

    # DEV
    index <- "B"
    tlim <- 10
    equil_time <- 12
    SurfEch <- 0.03
    # EoDev

    species <- Mosaic$forest[[index]]$species
    res <- new.env()

    res$index <- index
    res$SurfEch <- SurfEch
    res$tlim <- tlim

    res$sp <- unname(Mosaic$forest[[index]]$info$species)
    nsp <- length(res$sp)

    res$meshs <- map(Mosaic$ipms, ~ .x$mesh)
    res$stand_above_dth <- map2(res$meshs, species,  ~ .x > .y$harv_lim["dth"])
    res$BAsp <- map(Mosaic$ipms, ~ .x$BA)

    sim_X <- init_sim(nsp, tlim, res$meshs)
    res$sim_X <- do.call("rbind", sim_X)
    res$sim_BAstand <- res$sim_BAsp <- as.data.frame(matrix(
        ncol = nsp, nrow = tlim + 2, dimnames = list(NULL, res$sp)
    ))
    res$sim_BA <- rep(NA_real_, equil_time)
    res$sim_BAnonSp <- rep(NA_real_, nsp)

    # Initiate the pops
    res$X <- map2(map(species, `[[`, "init_pop"),
              res$meshs,
              exec, SurfEch = SurfEch)
    res$Harv <- map(lengths(res$meshs), ~ rep(0, .x))
    res$ct <- map(res$meshs, Buildct, SurfEch = SurfEch)

    return(res)
}

x <- res

save_step_env <- function(x, t){

    # Dev
    t <- 1
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
        stop(paste(
            "Border Basal Area reached for this simulation.",
            "This maximum is reached before iteration, check init_pop functions"
        ))
    }

    return(x)
}
