#' Add bandeau matrix as sub diagonal of a matrix
#'
#' @param matrix Matrix to fill in diagonal
#' @param fill Matrix to modify and paste inside matrix
#' @param dist Distance to the diagonal of the matrix. Default is 0 and
#'
#' @return
#' Modifief matrix
#'
#' @examples
#' matrix <- matrix(0, ncol = 6, nrow = 6)
#' fill <- matrix(1, ncol = 6, nrow = 2)
#' fill[,2:4 ] <- 1:6
#' sub_diag(matrix, fill, dist = 0)
#' sub_diag(matrix, fill, dist = 2)
#' sub_diag(matrix, fill, dist = 6)
#'
#' @import checkmate
#'
#' @keywords Internal
sub_diag <- function(matrix, fill, dist = 0) {

    assertMatrix(matrix)
    assertMatrix(fill, ncols = ncol(matrix), min.rows = 1)
    assertCount(dist)

    nrm <- nrow(matrix)
    if(dist+1 > nrm){ # nothing to do
        return(matrix)
    }

    nc <- ncol(matrix)
    nr <- nrow(fill)

    for (k in 1:(nc - dist)) {
        sel <- (k:min(k+nr-1, nrm-dist)) + dist
        matrix[sel, k] <- fill[1:length(sel), k]
    }

    return(matrix)
}


# clim <- c(sgdd = 2688, wai = -0.2)
# expand_clim <- function(clim){
#
#     nms <- names(clim)
#     sq <- clim ^ 2 ; names(sq) <- paste0(nms, 2)
#     b <- 1 / clim ; names(b) <- paste0(nms, "b")
#
#     res <- c(clim, sq, b)
# }


