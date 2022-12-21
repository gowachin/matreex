#' get Quadratic Mean Diameter
#'
#' Compute the quadratic mean diameter for a given size distribution.
#'
#' @param size Size class vector of the following distribution. dbl vector
#' @param n population state distribution at time t. dbl vector.
#'
#' @details
#' Both input share the same length.
#'
#' @export
QMD <- function(size, n){

    assertNumeric(size)
    assertNumeric(n, len = length(size))

    # IDEA Does not count delay in qmd ?
    n <- n[size > 0]
    size <- size[size > 0]

    res <- sqrt(sum(size * size * n) / sum(n))

    return(res)
}

