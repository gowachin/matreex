#' Add bandeau matrix as sub diagonal of a matrix
#'
#' @param matrix Matrix to fill in diagonal
#' @param fill Matrix to modify and paste inside matrix
#' @param dist Distance to the diagonal of the matrix. Default is 0 and
#' @param new This logical allow to initiate matrix with dim equal to ncol(fill)
#' and fill it with 0. This remove a duplicate on modify and diminish memory
#' usage during simulation.
#'
#' @return
#' Modifief matrix
#'
#' @examples
#' matrix <- matrix(9, ncol = 6, nrow = 6)
#' fill <- matrix(1, ncol = 6, nrow = 2)
#' fill[,2:4 ] <- 1:6
#' sub_diag(matrix, fill, dist = 0)
#' sub_diag(matrix, fill, dist = 2)
#' sub_diag(matrix = NULL, fill, dist = 2, new = TRUE)
#' sub_diag(matrix, fill, dist = 6)
#'
#' @import checkmate
#'
#' @noRd
sub_diag <- function(matrix, fill, dist = 0, new = FALSE) {

    assertMatrix(fill)
    nc <- ncol(fill)
    if(new){
        # matrix <- matrix(0, nrow = ncol(fill), ncol = ncol(fill))
        matrix <- numeric(nc^2)
        dim(matrix) <- c(nc, nc)
    } else {
        assertMatrix(matrix, ncols = nc, min.rows = 1)
    }
    assertCount(dist)

    nrm <- nrow(matrix)
    if(dist+1 > nrm){ # nothing to do
        return(matrix)
    }

    nr <- nrow(fill) - 1
    nrm_dist <- nrm-dist

    for (k in 1:(nc - dist)) {
        sel <- (k:min(k+nr, nrm_dist)) + dist
        matrix[sel, k] <- fill[seq_along(sel), k]
    }

    return(matrix)
}
