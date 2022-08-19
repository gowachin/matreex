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

    if(grepl("log", x)){
        evar <- sub("^log", "", x)
        evar <- call2("log", ensym(evar))
    } else if(grepl("2$", x)){
        evar <- sub("2$", "", x)
        evar <- call2("^", ensym(evar), 2)
    } else {
        evar <- ensym(x)
    }

    tmp <- call2("*", df[df$var1 == x, "K"], evar)
    new_name <- paste0(x, "_in")
    res <- call2("<-", ensym(new_name), tmp)

    return(res)
}


#' Format fitted derivated object to table
#'
#' Format values to be used in other functions
#'
#' @param params Estimated parameters for the fit of the model.
#' @param list_covs Climatic covariates values.
#'
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr pivot_longer everything separate
#' @importFrom tibble rownames_to_column
#'
#' @noRd
format_fit <- function(params, list_covs){

    invar <- names(params)[!names(params) %in% names(list_covs)]

    lc <- pivot_longer(list_covs, cols = everything())
    lc <- rbind(lc, data.frame(name= c(invar, NA),  value=1))
    p <- as.data.frame(params) %>% rownames_to_column(var = "var") %>%
        separate(var, c("var1", "var2"), sep = "\\:", fill = "right")
    res <- left_join(p, lc, by=c('var1'='name')) %>%
        left_join(lc, by=c('var2'='name')) %>%
        mutate(K = params * value.x * value.y)

    return(res)
}



#' Export recruitment function from estimated parameters.
#'
#' @param params Estimated parameters for the fit of the model.
#' @param list_covs Climatic covariates values.
#'
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @importFrom rlang expr call2
#'
#' @return
#' Function with 4 parameters : BATOTSP, BATOTNonSP, mesh and SurfEch
#'
#' @noRd
exp_recFun <- function(params, list_covs){

    df2 <- format_fit(params, list_covs)

    invar <- names(params)[!names(params) %in% names(list_covs)]
    invar <- invar[! grepl("ntercept", invar)]
    inter <- sum(filter(df2, ! var1 %in% invar | var2 %in% invar)$K)


    exp_invar <- map(invar, multi, df2)
    add_invar <- map(c(list(expr(intercept <- 1)),exp_invar),
                     ~ call2("<-", expr(res), call2("+", expr(res), .x[[2]] )))

    final_res <- list(
        expr(mesh <- length(mesh)),
        expr(distrib <- c(rep(1/2, 2), numeric(mesh - 2)) ),
        expr(final <- exp(res) * SurfEch / 0.03 * distrib ),
        expr(return(final))
    )
    calls <- c(exp_invar, add_invar, final_res)

    empty <- function(BATOTSP, BATOTNonSP, mesh, SurfEch = 0.003){}

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


