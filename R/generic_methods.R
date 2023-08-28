
# sp_name ####
#' List the species in matreex class object
#'
#' @param x matreex class object. Used on ipm and species at this time.
#' @name sp_name
#'
#' @export
sp_name <- function(x){
    UseMethod("sp_name")
}

#' @method sp_name ipm
#' @export
sp_name.ipm <- function(x){
    return(unname(x$info["species"]))
}

#' @method sp_name mu_gr
#' @export
sp_name.mu_gr <- function(x){
    return(unname(x$info["species"]))
}

#' @method sp_name species
#' @export
sp_name.species <- function(x){
    return(unname(x$info["species"]))
}

#' @method sp_name deter_sim
#' @export
sp_name.deter_sim <- function(x){
    nms <- rownames(x)
    nms <- nms[grepl(".*.BAsp", nms)]
    res <- sub("(.*)\\.(BAsp)", "\\1", nms)

    return(res)
}

# Climatic ####
#' Get climatic values in matreex class object
#'
#' @param x matreex class object. Used on ipm and species at this time.
#' @name climatic
#'
#' @export
climatic <- function(x){
    UseMethod("climatic")
}

#' @method climatic ipm
#' @export
climatic.ipm <- function(x){

    res <- unname(x$info["clim_lab"])
    return(res)
}

#' @method climatic mu_gr
#' @export
climatic.mu_gr <- function(x){

    res <- unname(x$info["clim_lab"])
    return(res)
}

#' @method climatic species
#' @export
climatic.species <- function(x){

    res <- unname(x$info["clim_lab"])
    return(res)
}



# Delay ####
#' Add delay states to a system
#'
#' Add delay states to an existing pop_state class object or an ipm class
#' object.
#'
#' @details
#' This function is a method that call \code{delay.ipm} internal function.
#'
#' @param x an object that require a delay addition
#' @param delay the number of time delay to add. single int. Its possible to
#' give negative values to remove delay from previously delayed objects. However
#' final delay must be positive (0 is acceptable).
#'
#' @name delay
#'
#' @export
delay <- function(x, delay = 0){
    UseMethod("delay")
}

#' @method delay numeric
#' @export
delay.numeric <- function(x, delay = 0){
    assert_integerish(delay, len = 1)

    if(delay == 0){
        return(x)
    }

    if(delay > 0){
        x <- c(rep(0, delay), x)
    } else {
        x <- x[-c(1:-delay)]
    }

    return(x)
}

#' @method delay ipm
#' @export
delay.ipm <- function(x, delay = 0){

    assert_integerish(delay, len = 1)
    if(delay == 0){
        return(x)
    }


    prev_d <- as.numeric(x$info["delay"])
    new_d <- delay + prev_d
    if(new_d < 0){
        stop(sprintf(
            "Negative delay is not possible for ipm objects. Minimal value here is %.f",
            prev_d))
    }
    x$info["delay"] <- as.character(new_d)

    if(as.logical(x$info["compress"])){
        x$IPM <- lapply(x$IPM, `*`, 1e-7)
        x$info["compress"] <- "FALSE"
    }
    x$IPM <- lapply(x$IPM, delay.dtCMatrix, delay)
    x$mesh <- delay(x$mesh, delay)

    # prev_d <- as.numeric(x$info["delay"])
    # x$info["delay"] <- as.character(delay + prev_d)

    return(validate_ipm(x))
}

#' @method delay species
#' @export
delay.species <- function(x, delay = 0){

    assert_integerish(delay, len = 1)
    if(delay == 0){
        return(x)
    }
    x$IPM <- delay(x$IPM, delay)
    return(x)
}

#' @method delay forest
#' @export
delay.forest <- function(x, delay = 0){

    assert_integerish(delay, len = 1)
    if(delay == 0){
        return(x)
    }

    x <- lapply(x$species, delay.species, delay)

    return(
        validate_forest(new_forest( species = x ))
    )
}

#' Delay dtCMatrix
#'
#' Adding a topleft corner to a matrix filled with 0.
#' @inheritParams delay
#' @method delay dtCMatrix
#' @examples
#' x <- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L), p = c(0L, 3L,  5L, 6L),
#'          Dim = c(3L, 3L), x = c(0.1, 0.2, 0.3, 0.1, 0.2, 0.1), uplo = "L", diag = "N")
#' x <- delay(x, 2)
#' delay(x, -1)
#' @importFrom Matrix Matrix
#' @export
delay.dtCMatrix <- function(x, delay = 0){

    assert_integerish(delay, len = 1)
    if(delay == 0){
        return(x)
    }

    size <- dim(x)[1]
    if(delay > 0){

        P <- cbind(matrix(0, nrow = size+delay, ncol = delay),
                   rbind(matrix(0, ncol=size, nrow=delay), x))
        # NOTE this move from dtCMatrix to dgCMatrix, is it an issue ?

        for (i in 1:delay) {
            P[i + 1, i] <- 1
        }

    } else {
        P <- as.matrix(x)[-c(1:-delay), -c(1:-delay)]
        P <- Matrix(P, sparse = TRUE)
    }

    return(P)
}


#' Delay dgCMatrix
#'
#' Adding a topleft corner to a matrix filled with 0.
#' @inheritParams delay
#' @method delay dgCMatrix
#' @examples
#' x <- new("dtCMatrix", i = c(0L, 1L, 2L, 1L, 2L, 2L), p = c(0L, 3L,  5L, 6L),
#'          Dim = c(3L, 3L), x = c(0.1, 0.2, 0.3, 0.1, 0.2, 0.1), uplo = "L", diag = "N")
#' x <- delay(x, 2)
#' delay(x, -1)
#' @export
delay.dgCMatrix <- delay.dtCMatrix


# Correction ####
#' Add correction states to a system
#'
#' Apply correction to a matrix by setting to 0 all right and bottom values
#'
#' @param x an object that require a correction addition
#' @param correction Type of correction, match between cut and none (default).
#' chr.
#' @name correction
#'
#' @export
correction <- function(x, correction = "none"){
    UseMethod("correction")
}

#' @method correction ipm
#' @export
correction.ipm <- function(x, correction = "none"){

    if(as.logical(x$info["compress"])){
        x$IPM <- lapply(x$IPM, `*`, 1e-7)
        x$info["compress"] <- "FALSE"
    }
    if(correction == "cut"){
        for (i in seq_along(x$IPM)) {
            x$IPM[[i]][, length(x$mesh)] <- 0
            x$IPM[[i]][length(x$mesh), ] <- 0
        }
    }

    return(validate_ipm(x))
}

#' @method correction species
#' @export
correction.species <- function(x, correction = "none"){

    x$IPM <- correction(x$IPM, correction)
    return(x)
}

#' @method correction forest
#' @export
correction.forest <- function(x, correction = "none"){

    x$species <- lapply(x$species, correction.species, correction)
    return(x)
}


#' @method correction mu_gr
#' @export
correction.mu_gr <- function(x, correction = "none"){
    # nothing to modify here !
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
#' @name sp_rec
#'
#' @export
sp_rec <- function(x, climatic, regional = FALSE){
    UseMethod("sp_rec")
}

#' @method sp_rec mu_gr
#' @export
sp_rec.mu_gr <- function(x, climatic, regional = FALSE){

    res <- exp_recFun(params = x$fit$rec$params_m, list_covs = climatic,
                      regional = regional) # TODO test this ?
    return(res)
}

#' @method sp_rec ipm
#' @export
sp_rec.ipm <- function(x, climatic, regional = FALSE){

    res <- exp_recFun(params = x$fit$rec$params_m, list_covs = climatic,
                      regional = regional)
    return(res)
}

#' @method sp_rec species
#' @export
sp_rec.species <- function(x, climatic, regional = FALSE){

    res <- sp_rec(x = x$IPM, climatic, regional = regional)
    return(res)
}


# get_maxdbh ####
#' Get the max_dbh of a fitted species
#'
#' @param x matreex class object. Used on fit_sgr for now.
#' @name get_maxdbh
#'
#' @export
get_maxdbh <- function(x){
    UseMethod("get_maxdbh")
}

#' @method get_maxdbh fit_sgr
#' @export
get_maxdbh.fit_sgr <- function(x){
    return(as.numeric(x$info["max_dbh"]))
}

# Summary ####
#' IPM summary
#'
#' Additionnal argument :
#' \describe{
#'  \item{fline}{skip the first line. FALSE by default.}
#' }
#'
#' @param object Ipm object.
#' @param ... Dots are dots.
#'
#' @method summary ipm
#' @export
summary.ipm <- function(object, ...){

    info <- object$info
    int <- object$int
    dim <- dim(object$IPM[[1]])
    BA <- object$BA
    mesh <- object$mesh

    dots <- list(...)
    fline <- dots$fline
    fline <- ifelse(is.null(fline), FALSE, fline)

    res <- sprintf(
        paste0(
            ifelse(fline, "%s%s",
                   "IPM object for species %s at climate '%s' \n\n"
            ),
            "Integation was done on a mesh from %.2f to %.2f with %.0f cells for BA between %.0f and %.0f. \n",
            "Gauss-Legendre was used on %.0f cells with %.0f x %.0f levels (size at t * t+1) \n",
            "Midbin was used on %.0f cells with %.0f levels. \n",
            "The correction was %s, the IPM does %scontain survival and the recruitment delay is %.0f. \n"
        ),
        ifelse(fline, "", info["species"]),
        ifelse(fline, "", info["clim_lab"]),
        mesh[1], tail(mesh, 1), length(mesh), BA[1], tail(BA, 1),
        int["gl_tresh"], int["gl1"], int["gl2"],
        int["mb_tresh"], int["mid_level"],
        info["correction"], ifelse(as.logical(info["surv"]), "", "not "),
        as.numeric(info["delay"])
    )

    message(res)

    invisible(object)
}

#' mu_gr summary
#'
#' This function is not written yet.
#'
#' @param object Ipm object.
#' @param ... Dots are dots.
#'
#' @method summary mu_gr
#' @export
summary.mu_gr <- function(object, ...){

    # TODO !
    res <- ""
    message(res)
    invisible(object)
}

#' species summary
#'
#' Additionnal argument :
#' \describe{
#'  \item{sum_ipm}{Print IPM summary. TRUE by default.}
#' }
#'
#' @param object Ipm object.
#' @param ... Dots are dots.
#'
#' @method summary species
#' @export
summary.species <- function(object, ...){

    info <- object$info
    harv_lim <- object$harv_lim
    rdi_coef <- object$rdi_coef
    disturb_coef <- object$disturb_coef

    dots <- list(...)
    sum_ipm <- dots$sum_ipm

    res <- sprintf(
        paste0(
            "IPM object for species %s at climate '%s' \n\n",
            "For Uneven harvest, dth = %.1f, dha = %.1f and hmax = %.2f.\n",
            "rdi_coef are %sdefined for this species. Even harvest %spossible\n",
            "disturb_coef are %sdefined for this species. Disturbance models %spossible\n"
        ),
        info["species"],  info["clim_lab"],
        harv_lim["dth"], harv_lim["dha"], harv_lim["hmax"],
        ifelse(is.null(rdi_coef), "not ", ""),
        ifelse(is.null(rdi_coef), "not ", ""),
        ifelse(is.null(disturb_coef), "not ", ""),
        ifelse(is.null(disturb_coef), "not ", "")
    )

    message(res)

    if(is.null(sum_ipm) || sum_ipm){
        summary(object$IPM, fline = TRUE)
    }
    invisible(object)
}
