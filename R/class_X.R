#' is.state
#'
#' State class is the set of all possible individual-level states and number of
#' individual per states in a population.
#'
#' @examples
#' mesh <- seq(1, 10, by = 1)
#' is.pop_state(state_init(mesh))
#'
#' @noRd
is.pop_state <- function(x){
    assertNumeric(x, lower = 0)

    if(length(attributes(x)$mesh) != length(x)){
        stop("mesh can't have a different length than state object x")
    }

    return(inherits(x, "pop_state"))
}

#' pop_state print method
#'
#' @examples
#' mesh <- seq(1, 10, by = 1)
#' print(state_init(mesh))
#'
#' @param x pop_state class (it is a numeric vector with attributes.)
#' @param ... Ignored
#'
#' @keywords internal
#' @method print pop_state
#' @export
print.pop_state <- function(x, ...){
    print(x[seq_along(x)])
    return(invisible(x))
}

#' pop_state initiation
#'
#' Initiate a population state by drawing from the mesh of the related IPM
#'
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#'
#' @examples
#' mesh <- seq(1, 10, by = 1)
#' state_init(mesh)
#'
#' @importFrom stats rbinom runif
#'
#' @export
state_init <- function(mesh){
    assertNumeric(mesh, lower = 0)

    x <- exp(runif(1,-.005,.005)*mesh)
    x[which(rbinom(length(mesh),1,runif(1,.6,.9)) == 1)] <- 0

    # attr(x, "mesh") <- mesh
    # class(x) <- c("pop_state", class(x))

    return(x)
}




