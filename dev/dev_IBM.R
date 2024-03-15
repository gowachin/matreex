load_all()
document()
data("fit_Picea_abies")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
climate

params <- fit_Picea_abies$sv$params_m
list_covs <- data.frame(wai = 0.5)
# list_covs <- climate

exp_allFun <- function(params, list_covs){

    (df2 <- format_fit(params, list_covs))

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(df2$K[! (df2$var1 %in% invar | df2$var2 %in% invar)])

    # multi("BATOTcomp", df2)
    # multi("size:sgdd", df2)
    # multi("size:wai", df2)
    # multi("wai2", df2)
    #
    # df <- df2
    # x <- "wai2"


    exp_invar <- map(invar, multi, df2)
    exp_invar
    invar <- unlist(map(exp_invar, ~ attributes(.x)$var))
    invar <- unique(invar[!is.na(invar)])
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    final_res <- list(
        expr(return(res))
    )
    calls <- c(exp_invar, add_invar, final_res)

    empty <- function(){}
    arguments <- setNames(rep(alist(x=), length(invar) + 1), c(invar, "..."))
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
    return(empty)
}

