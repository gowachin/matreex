#' get Quadratic Mean Diameter
#'
#' Compute the quadratic mean diameter for a given size distribution.
#'
#' @export
QMD <- function(size, n){

    assertNumeric(size)
    assertNumeric(n, len = length(size))

    res <- sqrt(sum(size * size * n) / sum(n))

    return(res)
}
