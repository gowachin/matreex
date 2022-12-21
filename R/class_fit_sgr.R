#' Constructor for fit_sgr class
#'
#' Only used in the matreex package
#'
#' @param sv_params Named vector of survival parameters fitted for this species
#' and climatic condition. Minimal parameters are intercept and size.
#' @param sv_family Family list contain the details of the models used when
#' fitting by functions such as glm. Only the inverse of the link function
#' (linkinv) is used in the model for now. list.
#' @param gr_params Named vector of growth parameters fitted for this species and
#' climatic condition. Minimal parameters are intercept and size.
#' @param gr_sigma Standard deviation of the residuals for the growth fitted
#' model.
#' @param rec_params Named vector of growth parameters fitted for this species and
#' climatic condition. Minimal parameters are intercept, BATOTSP and BATOTNonSP.
#' @param species Name of the species to run simulation on. Single char.
#' @param max_dbh Maximum diameter of the fitted dataset. Single dbh.
#'
#' @keywords internal
#' @export
new_fit_sgr <- function(sv_params, sv_family,
                        gr_params, gr_sigma,
                        rec_params,
                        species, max_dbh){
    fit <- list(
        sv = list(params_m = sv_params, family = sv_family),
        gr = list(params_m = gr_params, sigma = gr_sigma),
        rec = list(params_m = rec_params),
        info = c(species = species, max_dbh = max_dbh)
    )

    class(fit) <- "fit_sgr"

    return(fit)
}

#' validator for IPM class.
#'
#' @param x IPM class object
#'
#' @import checkmate
#'
#' @noRd
validate_fit_sgr <- function(x){

    assertClass(x, "fit_sgr")
    values <- unclass(x)
    names <- attr(x, "names")

    # check names of the object ####
    assertCharacter(names)
    if(any(names != c("sv", "gr", "rec", "info"))){
        stop(paste0("fit_sgr class must be composed of elements sv, gr, rec and info"))
    }

    # check all values ####
    assertList(values$sv, types = c("numeric", "list"), any.missing = FALSE,
               len = 2)
    tmp <- names(values$sv$params_m)
    assertCharacter(tmp, any.missing = FALSE)
    if(! all(c("intercept", "size") %in% tmp) ){
        stop("Survival model should at least depend on an intercept and the size variable.")
    }
    assertList(values$gr, types = c("numeric", "numeric"), any.missing = FALSE,
               len = 2)
    tmp <- names(values$gr$params_m)
    assertCharacter(tmp, any.missing = FALSE)
    if(! all(c("intercept", "size") %in% tmp) ){
        stop("Growth model should at least depend on an intercept and the size variable.")
    }

    tmp <- names(values$rec$params_m)
    assertCharacter(tmp, any.missing = FALSE)
    if(! all(c("intercept", "BATOTSP", "BATOTNonSP") %in% tmp) ){
        stop("Recruitment model should at least depend on an intercept, BATOTSP and BATOTNonSP variable.")
    }
    # check infos ####
    assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("species", "max_dbh"))){
        stop("fit_sgr class must have info of elements species and max_dbh")
    }

    invisible(x)
}

#' Create a new fitted models from
#'
#' Species are defined by an IPM which is a transition matrix from size between
#' t and t+1, recruitment and harvest functions. Each species has these items
#' defined for a given climate.
#' An additionnal vector of harvest parameers is required with minimal size to
#' harvest (dth), size above wich harvest is constant (dha).
#'
#' @inheritParams new_fit_sgr
#'
#' @export
fit_sgr <- function(sv_params, sv_family,
                    gr_params, gr_sigma,
                    rec_params,
                    species, max_dbh){

    res <- validate_fit_sgr(new_fit_sgr(
        sv_params = sv_params, sv_family = sv_family,
        gr_params= gr_params, gr_sigma = gr_sigma,
        rec_params = rec_params,
        species = species, max_dbh = max_dbh
    ))

    return(res)
}


#' load old fitted models
#'
#' @param species Name of the species to run simulation on. Single char.
#' @param path Place to save the resulting file. Single Char.
#' @param replicat Numeric for the simulation to select. By default, the 42th.
#' @param mean Should the return element be a mean of all models or a single value.
#' FALSE by default, TRUE will ignore replicat argument
#'
#' @import checkmate
#' @import here
#' @importFrom purrr map_dbl
#'
#' @export
old_fit2fit <- function(species, path = here(), replicat = 42, mean = FALSE){

    assertCharacter(species, len = 1)
    assertCharacter(path, len = 1)

    f_fit <- here(path, "output", species, "fit_sgr_all.Rds")
    fit <- readRDS(assertFileExists(f_fit)) # 0.16 sec
    if(grepl(" ", species)){
        species <- sub(" ", "_", species)
    }

    if(mean){
        max_dbh <- max(map_dbl(fit, ~ .x$maxDBH))
        res_fit <- mean_oldfit(fit, species, max_dbh)
    } else {
        assertNumber(replicat, lower = 1, upper = length(fit))
        res_fit <- fit[[replicat]]
        res_fit <- fit_sgr(res_fit$sv$params_m, res_fit$sv$family,
                           res_fit$gr$params_m, res_fit$gr$sigma,
                           res_fit$rec$params_m,
                           species = species, max_dbh = res_fit$maxDBH)
    }

    return(res_fit)
}


#' @noRd
#'
#' @param fit list of fitted models for survival and growth. Each element is a
#' list with sv$params_m, sp$family, gr$params_m and gr$sigma. All other
#' elements are discarded.
#' @param species single chr with the species name.
#' @param max_dbh maximum diameter of fitted data.
#'
#' @return a fit_sgr element.$
#'
#' @import checkmate
#' @importFrom purrr map reduce map_dbl
#' @importFrom dplyr full_join group_by summarize pull
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom rlang .data
#' @noRd
mean_oldfit <- function(fit, species, max_dbh){

    assertList(fit, any.missing = FALSE)
    assertCharacter(species, len = 1)
    assertNumber(max_dbh, lower = 0)

    sv_params <- map(fit, ~ data.frame(var = names(.x$sv$params_m),
                                       value = .x$sv$params_m)) %>%
        reduce(full_join, by = "var") %>%
        # replace(is.na(.data), 0) %>%
        pivot_longer(-.data$var) %>%
        replace_na(list(value = 0)) %>%
        group_by(.data$var) %>% summarize(mean = mean(.data$value)) %>%
        pull(.data$mean, .data$var)

    sv_family <- list(family = fit[[1]]$sv$family$family,
                      link = fit[[1]]$sv$family$link,
                      linkinv = fit[[1]]$sv$family$linkinv)
    class(sv_family) <- "family"

    gr_params <- map(fit, ~ data.frame(var = names(.x$gr$params_m),
                                       value = .x$gr$params_m)) %>%
        reduce(full_join, by = "var") %>%
        # replace(is.na(.data), 0) %>%
        pivot_longer(-.data$var) %>%
        replace_na(list(value = 0)) %>%
        group_by(.data$var) %>% summarize(mean = mean(.data$value)) %>%
        pull(.data$mean, .data$var)

    gr_sigma <- mean(map_dbl(fit, ~ .x$gr$sigma))

    rec_params <- map(fit, ~ data.frame(var = names(.x$rec$params_m),
                                       value = .x$rec$params_m)) %>%
        reduce(full_join, by = "var") %>%
        # replace(is.na(.data), 0) %>%
        pivot_longer(-.data$var) %>%
        replace_na(list(value = 0)) %>%
        group_by(.data$var) %>% summarize(mean = mean(.data$value)) %>%
        pull(.data$mean, .data$var)

    res_fit <- fit_sgr(sv_params, sv_family, gr_params, gr_sigma, rec_params,
                       species = species, max_dbh = max_dbh)

    return(res_fit)
}
