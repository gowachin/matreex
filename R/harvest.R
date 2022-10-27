# UNEVEN ####
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
#' @param ... Variables used in this case of Uneven harvest
#' \describe{
#'  \item{targetBAcut}{Basal Area to cut.}
#'  \item{ct}{is the vector to compute BA with x (ct = Buildct(mesh, SurfEch))}
#' }
#'
#' @family functions that defines harvest rules for a species.
#'
#' @export
Uneven_harv <- function(x,
                        species,
                        ...) {

    dots <- list(...)
    targetBAcut <- dots$targetBAcut
    ct <- dots$ct

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
# EVEN ####
#' Compute RDI for a species
#'
#' @param x Size distribution of a species at time t
#' @param RDI_int RDI curve intercept for this species
#' @param RDI_slo RDI curve slope for this species
#' @param meshcm2  Vector of all values between L and U of the IPM
#' (L being the smallest size and U the largest.). Length of this vector is
#' the number of size class in the IPM. This vector is in cm and squared. num.
#' @param tx Total effectif of individual.
#'
#' @details
#' Mesh needs to be in cm squared.
#' Please do this before this function to reduce computation time.
#'
#' @noRd
RDI <- function(x, RDI_int, RDI_slo, meshcm2, tx){
    res <- sum(x) / exp(RDI_int + RDI_slo / 2 * log(meshcm2 %*% x / tx))

    return(drop(res))
}

#' Compute the harvest curve for even stand (Pcut)
#'
#' @param x is the size distribution
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param RDIcoef is a one line dataframe with RDI coefficient for one species
#' @param targetKg is the target Kg
#' @param targetRDI is the target RDI
#' @param SurfEch Value of plot size surface in ha
#'
#' #' @details
#' If RDI < targetRDI, this function will return a Pcut of 0 for all sizes.
#'
#' @importFrom stats optim
#'
#' @noRd
getPcutEven <- function(x,
                        mesh,
                        RDIcoef,
                        targetRDI=0.65,
                        targetKg=0.9,
                        SurfEch = 0.03){

    assertNumeric(x, lower = 0)
    assertNumeric(mesh, lower = 0, len = length(mesh) )
    assertNumber(targetRDI, lower = 0)
    assertNumber(targetKg, lower = 0)

    x <- x / SurfEch
    N <- length(mesh)
    tx <- sum(x)
    meshcm2 <- (mesh / 10) ^ 2
    RDI_int <- RDIcoef[["intercept"]]
    RDI_slo <- RDIcoef[["slope"]]

    # Early return if RDI < targetRDI
    rdi <- RDI(x, RDI_int, RDI_slo, meshcm2, tx)
    if (rdi < targetRDI) {
      # cat(sprintf("rdi: %.2f\n", rdi))
      return(rep(0, N))
    }

    getCostFunc <- function(Par, meshcm2, N, tx, RDI_int, RDI_slo, x,
                            targetRDI, targetKg) {
      Dg2 <- meshcm2 %*% x / tx
      hmax <- Par[1]
      k <- Par[2]
      Pc <- hmax * (1:N)^(-k)
      if(sum(x*Pc) <= 0){
          return(.9) # escape if no harvest
      }
      # cat(sprintf("\nHarv: %.2f\n", sum(x*Pc)))
      Dgcut2 <- meshcm2 %*% (x * Pc) / sum(x * Pc)
      Kg <- Dgcut2 / Dg2
      RDI <- RDI(x = x * (1 - Pc), RDI_int, RDI_slo, meshcm2, tx)
      # cat(sprintf(
          # "hmax %.5f, k %.5f, error : %.5f \nKg : %.2f / Target :%.2f \nRDI : %.2f / Target :%.2f \n",
          # hmax, k,  (Kg-targetKg)^2 + (RDI-targetRDI)^2, Kg, targetKg, RDI, targetRDI)
      # )
      return((Kg - targetKg)^2 + (RDI - targetRDI)^2)
    }

    ParOpt <- optim(
      fn = getCostFunc,
      meshcm2 = meshcm2, N = N, tx = tx,
      RDI_int = RDI_int, RDI_slo = RDI_slo,
      x = x, targetRDI = targetRDI, targetKg = targetKg,

      par = c(0.1, 0.1), method = "L-BFGS-B",
      upper = c(.99, 1), lower = c(0, 1e-10)
    )$par

    Pcut <- ParOpt[1] * (1:N)^(-ParOpt[2])

    return(Pcut)
}


#' Perform the harvest Even
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param ... Variables used in this case of Uneven harvest
#' \describe{
#'  \item{targetKg}{is the target Kg}
#'  \item{targetRDI}{is the target RDI}
#'  \item{ct}{is the vector to compute BA with x (ct = Buildct(mesh, SurfEch))}
#'  \item{SurfEch}{Value of plot size surface in ha}
#' }
#'
#' @family functions that defines harvest rules for a species.
#'
#' @export
Even_harv <- function(x,
                      species,
                      ...) {

    dots <- list(...)
    targetRDI <- dots$targetRDI
    targetKg <- dots$targetKg
    ct <- dots$ct
    SurfEch <- ifelse(is.null(dots$SurfEch), 0.03, dots$SurfEch)

    Pcut <- getPcutEven(
        x = x,
        mesh = species$IPM$mesh,
        RDIcoef = species$rdi_coef,
        targetRDI = targetRDI, targetKg = targetKg,
        SurfEch = SurfEch
    )

    return(x * Pcut)
}
