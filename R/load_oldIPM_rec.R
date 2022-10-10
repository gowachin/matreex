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

    assertCharacter(x, any.missing = FALSE, len = 1)
    assertDataFrame(df, min.rows = 1, ncols = 6)

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
        xvar <- call2("log", ensym(xvar))
    } else if(grepl("2$", x)){
        xvar <- sub("2$", "", x)
        xvar <- call2("^", ensym(xvar), 2)
    } else {
        xvar <- ensym(x)
    }

    if(grepl("log", y)){
        yvar <- sub("^log", "", y)
        foo <- "log"
    } else if(grepl("2$", y)){
        yvar <- sub("2$", "", y)
        foo <- "^"
    } else {
        foo <- ""
        yvar <- y
    }
    if(!is.na(y) && df$value.y[selec] != 1){
        yvar <- df$value.y[selec]
    } else {
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
    res

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
    invar <- nms[!nms %in% names(list_covs)]

    lc <- c(drop(as.matrix(list_covs)),
            setNames(rep(1, length(invar)), invar))

    x <- sub(":.*$", "", nms)
    y <- sub("^[[:alnum:]]*:?", "", nms)
    y[nchar(y) == 0] <- NA_character_
    res <- data.frame(var1 = x, var2 = y, params = params,
                      value.x = lc[x], value.y = lc[y],
                      row.names = NULL)
    value.x <- NULL # HACK rm the note in devtools::check() about unbinded
    res <- within(res, {
        value.y[is.na(value.y)] <- 1
        K <- params * value.x * value.y
    })

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
#' @importFrom dplyr filter
#' @importFrom rlang expr call2 env_unbind .data
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
exp_recFun <- function(params, list_covs){

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(filter(df2, ! .data$var1 %in% invar | .data$var2 %in% invar)$K)


    exp_invar <- map(invar, multi, df2)
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    SurfEch <- NULL # HACK rm the note in devtools::check() about unbinded
    final_res <- list(
        expr(mesh <- length(mesh)),
        expr(distrib <- c(rep(1/2, 2), numeric(mesh - 2)) ),
        expr(final <- exp(res) * SurfEch / 0.03 * distrib ),
        expr(return(final))
    )
    calls <- c(exp_invar, add_invar, final_res)

    empty <- function(BATOTSP, BATOTNonSP, mesh, SurfEch = 0.03){}

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
#' @importFrom purrr map cross map_chr simplify_all
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
    combination <- cross(list(nms, nms)) %>%
        simplify_all() %>%
        map_chr(~ if(.x[1] == .x[2]){
            paste0(.x[1],"2")
        } else {
            paste0(.x, collapse = ":")
        })

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
    empty
    return(empty)
}

