
#' BAstand
#'
#' Compute BA standing for a species
#'
#' @param X Size distribution at time t
#' @param species The species class object of interest to get mesh and harv_lim
#' values from.
#' @param SurfEch Value of plot size surface in ha
#'
#' @return
#' Basal area above dth for a species
#'
#' @noRd
getBAstand <- function(X, species, SurfEch = 0.03){

    mesh <- species$IPM$mesh
    dth <- species$harv_lim["dth"]

    X[mesh <= dth] <- 0

    BAst <- drop(X %*% Buildct(mesh, SurfEch))

    return(BAst)
}

#' Compute the basal area to be cut to reach targetBA
#' @param BAst Standing basal area
#' @param targetBA the basal area to reach after maximum cut
#' @param Pmax the maximum proportion of BAcut / BA
#' @param dBAmin the minimum BA to perform cut
#'
#' @noRd
getBAcutTarget <- function(BAst,
                           targetBA = 20,
                           Pmax = 0.25,
                           dBAmin = 3) {
    assertNumber(BAst, lower = 0)
    assertNumber(targetBA, lower = 0)
    assertNumber(Pmax, lower = 0, upper = 1)
    assertNumber(dBAmin, lower = 0)

    if (BAst == 0 || (BAst - targetBA) < dBAmin) {
        BAcut <- 0
    } else {
        Htarget <- (BAst - targetBA) / BAst
        H <- min(Pmax, Htarget)
        BAcut <- BAst * H
    }

    return(BAcut)
}

#' Compute the harvest curve for uneven stand (Pcut)
#'
#' @param x is the size distribution
#' @param h is the simple harvest function (power=1)
#' @param ct is the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
#' @param hmax is the maximum harvest rate for a size class
#' @param targetBAcut Basal Area to cut.
#'
#' @details
#' Two cases:
#' \describe{
#'  \item{}{BA above dha is enough to reach target and we only cut above
#'          dbh such as we get target}
#'  \item{}{BA above dha is not enough, we optimize power k of h
#'          function to get target}
#' }
#'
#' @importFrom stats optimize
#'
#' @noRd
getPcutUneven <- function(x,
                          h,
                          ct,
                          hmax = 1,
                          targetBAcut = 0) {
    assertNumeric(x, lower = 0, len = length(ct))
    assertNumeric(ct, lower = 0, len = length(x))
    assertNumber(hmax, lower = 0, upper = 1)
    assertNumber(targetBAcut, lower = 0)

    BAha <- ct %*% (hmax * (h == 1) * x)

    getCostfunc <- function(k, x, h, ct, hmax = 1, targetBAcut) {
        BAcut <- ct %*% (h^k * x * hmax)
        res <- (BAcut - targetBAcut)^2
        return(res)
    }

    if (BAha >= targetBAcut) {
        Pcut <- (cumsum(rev(ct * x * hmax)) - targetBAcut) <= 0
        Pcut <- hmax * rev(Pcut)
    } else {
        kopt <- optimize(
            f = getCostfunc, lower = 0.1, upper = 10, x = x, h = h,
            ct = ct, hmax = hmax, targetBAcut = targetBAcut
        )$minimum
        Pcut <- h^kopt * hmax
    }
    # if (verbose){
    #     message('BA initial : ', round(ct %*% x, digits=1),
    #             ' / H target : ', round(BAcuttarget/(ct%*%x), digits=2))
    #     message('BA cut : ', round(ct %*% (x*Pcut), digits=1),
    #             ' / Target : ', round(BAcuttarget, digits=1))
    # }
    return(Pcut)
}


#' Perform the harvest Uneven
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and harv_lim
#' values from.
#' @param targetBAcut Basal Area to cut.
#' @param ct define this # TODO
#'
#' @export
Uneven_harv <- function(x,
                        species,
                        targetBAcut,
                        ct) {

    mesh <- species$IPM$mesh
    dth <- species$harv_lim["dth"]
    h <- (mesh - dth) / (species$harv_lim["dha"] - dth)
    h <- pmax(0, h)
    h <- pmin(1, h)

    # we only consider tree above dth in the algorithm
    x1 <- x
    x1[species$IPM$mesh <= dth] <- 0

    Pcut <- getPcutUneven(
        x = x1,
        h = h,
        ct = ct,
        hmax = species$harv_lim["hmax"],
        targetBAcut = targetBAcut
    )
    return(x * Pcut)
}
