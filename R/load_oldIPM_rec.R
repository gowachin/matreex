#' Create call from old name
#'
#' Get variable name and value in fitted derivated object and return a call
#'
#' This function is used to load old IPM and modify the recruitment function
#' to use more arguments (BAtot and BAsp) in case of multiple species
#'
#' @param x Name of the argument. Can be preceded with log
#' (log variable is computed) or end with 2 (square)chr.
#' @param df2 Dataframe with x in var1 column and a column named K with numeric
#' values. Data.frame
#'
#' @examples
#' df <- data.frame(var1 = "BATOTSP", var2 = NA_character_, params = -0.0179,
#' value.x = 1, value.y = 1, K = -0.0179)
#' multi(x = "BATOTSP", df)
#' #> BATOTSP_in <- -0.0179 * BATOTSP
#'
#' # Square variable
#' df <- data.frame(var1 = "BATOTSP2", var2 = NA_character_, params = -0.0179,
#' value.x = 1, value.y = 1, K = -0.0179)
#' multi(x = "BATOTSP2", df)
#' #> BATOTSP2_in <- -0.0179 * BATOTSP ^ 2
#'
#' @importFrom rlang ensym call2
#' @import checkmate
#' @noRd
multi <- function(x, df){

    # assertString(x)
    # assertDataFrame(df, min.rows = 1, ncols = 6)

    new_name <- paste0(x, "_in")

    if(grepl(":", x)){
        y <- sub("^.*:", "", x)
        x <- sub(":.*$", "", x)
        selec <- df$var1 == x & df$var2 == y & !is.na(df$var2)
        if(all(!selec)){
            stop("Can't write expression if param is not present in formula !")
        }
    } else {
        selec <- NULL
        y <- NA_character_
        if(! x %in% df$var1){
            stop("Can't write expression if param is not present in formula !")
        }
    }

    if(grepl("log", x)){
        xvar <- sub("^log", "", x)
        px <- xvar
        xvar <- call2("log", ensym(xvar))
    } else if(grepl("2$", x)){
        xvar <- sub("2$", "", x)
        px <- xvar
        xvar <- call2("^", ensym(xvar), 2)
    } else {
        px <- x
        xvar <- ensym(x)
    }

    if(grepl("log", y)){
        yvar <- sub("^log", "", y)
        py <- yvar
        foo <- "log"
    } else if(grepl("2$", y)){
        yvar <- sub("2$", "", y)
        py <- yvar
        foo <- "^"
    } else {
        foo <- ""
        py <- y
        yvar <- y
    }
    if(!is.na(y) && df$true.y[selec]){
        yvar <- df$value.y[selec]
    } else {
        py <- yvar
        yvar <- ensym(yvar)
    }
    yvar <- switch(foo,
                   log = call2("log", yvar),
                   '^' = call2("^", yvar, 2),
                   yvar )

    if(!is.na(y)){
        evar <- call2("*", xvar, yvar)
        tmp <- call2("*", df$params[selec], evar)
    } else {
        tmp <- call2("*", df$params[df$var1 == x & is.na(df$var2)], xvar)
    }
    res <- call2("<-", ensym(new_name), tmp)
    # res
    attributes(res)$var <- c(px, py)
    return(res)
}


#' Format fitted derivated object to table
#'
#' Format values to be used in other functions
#'
#' @param params Estimated parameters for the fit of the model.
#' @param list_covs Climatic covariates values.
#'
#' @importFrom stats setNames
#'
#' @noRd
format_fit <- function(params, list_covs){

    nms <- names(params)
    unparams <- unname(params)
    invar <- nms[!nms %in% names(list_covs)]

    lc <- c(drop(as.matrix(list_covs)),
            setNames(rep(1, length(invar)), invar))

    x <- sub(":.*$", "", nms)
    y <- sub("^[[:alnum:]]*:?", "", nms)
    y[nchar(y) == 0] <- NA_character_

    vx <- unname(lc[x])
    vy <- unname(lc[y])
    fy <- is.na(vy) | !y %in% names(list_covs)
    vy[fy] <- 1
    K <- unparams * vx * vy

    # HACK List is Fast and Furious. Graou !
    res <- list(var1 = x, var2 = y, params = unparams,
                value.x = vx, value.y = vy, true.y = !fy, K = K
    )
    class(res) <- "data.frame"
    attributes(res)$row.names <- 1:length(x)

    return(res)
}

#' Export recruitment function from estimated parameters.
#'
#' Rebuild the function to use BASp and BAnonSp for a species.
#'
#' @param params Estimated parameters for the fit of the model.
#' @param list_covs Climatic covariates values.
#' @param regional The recruitment is now dependant on the regional basal area
#' by requiring another entry BAFecSP.
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
exp_recFun <- function(params, list_covs, regional = FALSE){

    # params <- c(intercept = -0.864, BATOTSP = -0.018, BATOTNonSP = 42,
    #             logBATOTSP = 12, sgddb = 286.813,
    # wai = -0.057, wai2 = 0.288 )
    # list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)


    if(regional){
        names(params) <- sub("^logBATOTSP$", "logBAFecSP", names(params))
    }

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(df2$K[! (df2$var1 %in% invar | df2$var2 %in% invar)])


    exp_invar <- map(invar, multi, df2)
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    SurfEch <- NULL # hack rm the note in devtools::check() about unbinded
    final_res <- list(
        expr(mesh <- length(mesh)),
        expr(distrib <- c(rep(1/2, 2), numeric(mesh - 2)) ),
        expr(final <- exp(res) * SurfEch / 0.03 * distrib ),
        expr(return(final))
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
    arguments <- setNames(rep(alist(x=), length(invar) + 2 + !regional),
                          # I really tried to do it nicely but...R is for regrets.
                          c("BATOTSP", if(regional){"BAFecSP"}else{ NULL},
                            "BATOTNonSP", "mesh", "SurfEch"))
                          # c(invar, "mesh", "SurfEch"))
    arguments$SurfEch <- 0.03
    formals(empty) <- arguments

    setNames(rep(alist(x=), length(invar) + 2), c(invar, "mesh", "SurfEch"))
    setNames(vector("pairlist", length(invar) + 2), c(invar, "mesh", "SurfEch"))

    body(empty)[[2]] <- call2("<-", expr(intercept), inter)
    body(empty)[[3]] <- expr(res <- 0)
    for(i in seq_along(calls)){
        body(empty)[[i + 3]] <- calls[[i]]
    }

    # this is messy and is to remove binded env.
    env_unbind(env = environment(empty), c("i", "calls", "final_res",
                                           "add_invar", "exp_invar",
                                           "inter", "invar", "df2", "SurfEch"),
               inherit = FALSE)
    # empty
    return(empty)
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
#' Function with 1 parameter : size
#'
#' @examples
#' params <- c(intercept = -0.864, size = -0.018, sgddb = 286.813,
#' wai = -0.057, wai2 = 0.288 )
#' list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)
#'
#' foo <- exp_sizeFun(params, list_covs)
#' foo
#' foo(1, 2, 1:5, 0.03)
#'
#' @noRd
exp_sizeFun <- function(params, list_covs){

    df2 <- format_fit(params, list_covs)

    # use this when we have sgdd:wai for example in params and not in list_covs
    # NOTE does not cover sgdd2:wai yet...but not needed from all data(fit_species)
    nms <- names(list_covs)
    tmp <- paste(nms, rep(nms, each = length(nms)), sep = ":")
    sel <- tmp %in% paste0(nms, ":", nms)
    tmp[sel] <- sub(":.*", "2", tmp[sel])
    combination <- tmp

    exp_list_covs <- unique(c(nms, combination))

    invar <- names(params)[!names(params) %in% exp_list_covs]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(with(df2, K[! var1 %in% invar | var2 %in% invar]))

    exp_invar <- map(invar, multi, df2)
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    final_res <- list( expr(return(res)) )
    calls <- c(exp_invar, add_invar, final_res)
    empty <- function(size){}

    body(empty)[[2]] <- call2("<-", expr(intercept), inter)
    body(empty)[[3]] <- expr(res <- 0)
    body(empty)[seq_along(calls)+3] <- calls
    # this is messy and is to remove binded env.
    env_unbind(env = environment(empty), c("calls", "final_res",
                                           "add_invar", "exp_invar",
                                           "inter", "invar", "df2"),
               inherit = FALSE)
    # empty
    return(empty)
}


#' Export vital function from estimated parameters.
#'
#' Rebuild the function to use size, BAsp and all for a species.
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
#' Objectif is to use it to replace previous exp_...Fun() with a generic one.
#'
#' @return
#' Function with many parameter : size, BATOTcomp, BATOTSP etc. Use ...
#'
#' @examples
#' params <- c(intercept = -0.864, size = -0.018, sgddb = 286.813,
#' wai = -0.057, wai2 = 0.288 )
#' list_covs <- data.frame(wai = -0.187, sgddb = 0, waib = 1.23, wai2 = 0.34)
#'
#' foo <- exp_sizeFun(params, list_covs)
#' foo
#' foo(1, 2, 1:5, 0.03)
#'
#' @noRd
exp_allFun <- function(params, list_covs){

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(df2$K[! (df2$var1 %in% invar | df2$var2 %in% invar)])


    exp_invar <- map(invar, multi, df2)
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

