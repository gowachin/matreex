
# sp_name ####
#' List the species in treeforce class object
#'
#' @param x treeforce class object. Used on ipm and species at this time.
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
#' Get climatic values in treeforce class object
#'
#' @param x treeforce class object. Used on ipm and species at this time.
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
#' @param delay the number of time delay to add. single positive int.
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
    assertCount(delay)

    x <- c(rep(0, delay), x)
    return(x)
}

#' @method delay ipm
#' @export
delay.ipm <- function(x, delay = 0){

    assertCount(delay)
    if(delay == 0){
        return(x)
    }
    if(as.logical(x$info["compress"])){
        x$IPM <- lapply(x$IPM, `*`, 1e-7)
        x$info["compress"] <- "FALSE"
    }
    x$IPM <- lapply(x$IPM, delay.dtCMatrix, delay)
    x$mesh <- delay(x$mesh, delay)

    prev_d <- as.numeric(x$info["delay"])
    x$info["delay"] <- as.character(delay + prev_d)

    return(validate_ipm(x))
}

#' @method delay species
#' @export
delay.species <- function(x, delay = 0){

    assertCount(delay)
    if(delay == 0){
        return(x)
    }
    x$IPM <- delay(x$IPM, delay)
    return(x)
}

#' @method delay forest
#' @export
delay.forest <- function(x, delay = 0){

    assertCount(delay)
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
#'          Dim = c(3L, 3L), x = c(1, 2, 3, 1, 2, 1), uplo = "L", diag = "N")
#' delay(x, 2)
#' @export
delay.dtCMatrix <- function(x, delay = 0){

    assertCount(delay)
    if(delay == 0){
        return(x)
    }

    size <- dim(x)[1]
    P <- cbind(matrix(0, nrow = size+delay, ncol = delay),
               rbind(matrix(0, ncol=size, nrow=delay), x))
    # NOTE this move from dtCMatrix to dgCMatrix, is it an issue ?

    for (i in 1:delay) {
        P[i + 1, i] <- 1
    }

    return(P)
}


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

