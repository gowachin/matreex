#' List the species in treeforce class object
#'
#' @name species
#'
#' @export
species <- function(x){
    UseMethod("species")
}

#' @rdname species
#' @export
species.ipm <- function(x){
    return(unname(x$info["species"]))
}


#' Get climatic values in treeforce class object
#'
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
