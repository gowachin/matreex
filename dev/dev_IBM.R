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
