
#' Compute the basal area to be cut to reach targetBA
#' @param x the size distribution
#' @param ct the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
#' @param targetBA the basal area to reach after maximum cut
#' @param Pmax the maximum proportion of BAcut / BA
#' @param dBAmin the minimum BA to perform cut
#'
#' @noRd
getBAcutTarget <- function(x,
                           ct,
                           targetBA = 20,
                           Pmax = 0.25,
                           dBAmin = 3) {
    assertNumeric(x, lower = 0, len = length(ct))
    assertNumeric(ct, lower = 0, len = length(x))
    assertNumber(targetBA, lower = 0)
    assertNumber(Pmax, lower = 0, upper = 1)
    assertNumber(dBAmin, lower = 0)

    BA <- as.numeric(ct %*% x)
    if (BA == 0 || (BA - targetBA) < dBAmin) {
        BAcut <- 0
    } else {
        Htarget <- (BA - targetBA) / BA
        H <- min(Pmax, Htarget)
        BAcut <- BA * H
    }

    return(BAcut)
}

#' Compute the harvest curve for uneven stand (Pcut)
#'
#' @param x is the size distribution
#' @param h is the simple harvest function (power=1)
#' @param ct is the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
#' @param hmax is the maximum harvest rate for a size class
#' @param targetBA is the basal area to reach after maximum cut
#' @param Pmax is the maximum proportion of BAcut / BA
#' @param dBAmin is the minimum BA to perform cut
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
                          targetBA = 20,
                          Pmax = 0.25,
                          dBAmin = 3) {
    assertNumeric(x, lower = 0, len = length(ct))
    assertNumeric(ct, lower = 0, len = length(x))
    assertNumber(hmax, lower = 0, upper = 1)
    assertNumber(targetBA, lower = 0)
    assertNumber(Pmax, lower = 0, upper = 1)
    assertNumber(dBAmin, lower = 0)

    BAcuttarget <- getBAcutTarget(
        x = x, ct = ct, targetBA = targetBA, Pmax = Pmax,
        dBAmin = dBAmin
    )

    BAha <- ct %*% (hmax * (h == 1) * x)

    getCostfunc <- function(k, x, h, ct, hmax = 1, BAcuttarget) {
        BAcut <- ct %*% (h^k * x * hmax)
        res <- (BAcut - BAcuttarget)^2
        return(res)
    }

    if (BAha >= BAcuttarget) {
        Pcut <- (cumsum(rev(ct * x * hmax)) - BAcuttarget) <= 0
        Pcut <- hmax * rev(Pcut)
    } else {
        kopt <- optimize(
            f = getCostfunc, lower = 0.1, upper = 10, x = x, h = h,
            ct = ct, hmax = hmax, BAcuttarget = BAcuttarget
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
#' @param harv_rule Vector for harvest rules at the scale of the forest.
#' \describe{
#'   \item{Pmax}{maximum proportion of BAcut / BA}
#'   \item{dBAmin}{the minimum BA to perform cut}
#'   \item{freq}{Frequence at which the harvest will be executed.}
#' }
#' @param targetBA the basal area to reach after maximum cut
#' @param ct define this # TODO
#' @param t time t of simulation
#'
#' @export
Uneven_harv <- function(x,
                        species,
                        harv_rule,
                        targetBA,
                        ct,
                        t) {
    if (t %% harv_rule["freq"] != 0) { # Not the right time to cut !
        return(rep(0, length(x)))
    }

    mesh <- species$IPM$mesh
    dth <- species$harv_lim["dth"]
    h <- (mesh - dth) / (species$harv_lim["dha"] - dth)
    h <- pmax(0, h)
    h <- pmin(1, h)

    # we only consider tree above dth in the algorithm
    x1 <- x
    x1[species$IPM$mesh <= dth] <- 0

    Pcut <- getPcutUneven(
        x = x1, h = h, ct = ct,
        hmax = species$harv_lim["hmax"],
        targetBA = targetBA,
        Pmax = harv_rule["Pmax"],
        dBAmin = harv_rule["dBAmin"]
    )
    return(x * Pcut)
}
