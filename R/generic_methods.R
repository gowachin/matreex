
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

#' @rdname sp_name
#' @export
sp_name.ipm <- function(x){
    return(unname(x$info["species"]))
}

#' @rdname sp_name
#' @export
sp_name.species <- function(x){
    return(unname(x$info["species"]))
}

#' @rdname sp_name
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

#' @rdname climatic
#' @export
climatic.ipm <- function(x){

    res <- as.numeric(
        unname(x$info["climatic"])
    )
    return(res)
}

#' @rdname climatic
#' @export
climatic.species <- function(x){

    res <- as.numeric(
        unname(x$info["climatic"])
    )
    return(res)
}


# Delay ####
#' Add delay states to a system
#'
#' Add delay states to an existing pop_state class object or an ipm class
#' object.
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

#' @rdname delay
#'
#'@export
delay.numeric <- function(x, delay = 0){
    assertCount(delay)

    x <- c(rep(0, delay), x)
    return(x)
}

#' @rdname delay
#'
#' @examples
#' mesh <- seq(1, 10, by = 1)
#' s <- state_init(mesh)
#' y <- delay(s, 5)
#'
#' @export
delay.pop_state <- function(x, delay = 0){
    NextMethod(x)
}


#' @rdname delay
#'
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

    return(validate_ipm(x))
}

#' @rdname delay
#'
#' @export
delay.species <- function(x, delay = 0){

    assertCount(delay)
    if(delay == 0){
        return(x)
    }
    x$IPM <- delay(x$IPM, delay)
    return(x)
}

#' @rdname delay
#'
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

#' @rdname delay
#'
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

#' @rdname correction
#'
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

#' @rdname correction
#'
#' @export
correction.species <- function(x, correction = "none"){

    x$IPM <- correction(x$IPM, correction)
    return(x)
}

#' @rdname correction
#'
#' @export
correction.forest <- function(x, correction = "none"){

    x$species <- lapply(x$species, correction.species, correction)
    return(x)
}

