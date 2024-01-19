# UNEVEN ####
#' BAstand
#'
#' Compute BA standing for a species per ha.
#'
#' @param X Size distribution at time t
#' @param species The species class object of interest to get mesh and harv_lim
#' values from.
#' @param SurfEch Value of plot size surface in ha
#'
#' @return
#' Basal area per ha above dth for a species
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
RDI <- function(x, RDI_int, RDI_slo, meshcm2){

    # print(RDI_int)
    # print(RDI_slo)
    # print(tx)
    res <- sum(x) / exp(RDI_int + RDI_slo / 2 * log(meshcm2 %*% x / sum(x)))

    return(drop(res))
}


#' @noRd
RDI_sp <- function(x, species, SurfEch){

    meshcm2 <- (species$IPM$mesh / 10) ^2
    lag <- meshcm2 > 0

    x <- x[lag]
    meshcm2 <- meshcm2[lag]
    x <- x / SurfEch
    tx <- sum(x)
    RDIcoef <- species$rdi_coef
    RDI_int <- RDIcoef[["intercept"]]
    RDI_slo <- RDIcoef[["slope"]]

    res <- tx / exp(RDI_int + RDI_slo / 2 * log(meshcm2 %*% x / tx))

    return(drop(res))
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
#' @details
#' Even harvest is done at the stand level with the getPcutEven function
#' to acomodate with the multispecific case.
#'
#'
#' @export
Even_harv <- function(x,
                      species,
                      ...) {

    return(x)
}


#' Compute the harvest curve for even stand (Pcut)
#'
#' @param x is the size distribution
#' @param sp species list for the rdi intercept and slope values.
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param RDIcoef is a one line dataframe with RDI coefficient for one species
#' @param targetKg is the target Kg
#' @param targetRDI is the target RDI
#' @param SurfEch Value of plot size surface in ha
#'
#' @details
#' During simulation, if RDI < targetRDI, this function is not called and
#' Pcut will be 0 for all sizes.
#'
#' @importFrom stats optim
#'
#' @noRd
getPcutEven <- function(x,
                        sp, meshs,
                        targetRDI=0.65,
                        targetKg=0.9,
                        SurfEch = 0.03){

    # assertNumeric(x, lower = 0)
    assertNumber(targetRDI, lower = 0)
    assertNumber(targetKg, lower = 0)

    RDIcoef = map(sp, ~.x$rdi_coef)

    meshcm2 <- map(meshs, ~ (.x / 10) ^ 2)
    lag <- map(meshcm2, ~ .x > 0)
    del <- map_dbl(lag, ~ sum(!.x))

    x <- map2(x, lag, ~ .x[.y])
    meshcm2 <- map2(meshcm2, lag, ~ .x[.y])

    x <- map(x, ~ .x / SurfEch)
    N <- lengths(meshcm2)

    RDI_int <- map_dbl(sp, ~ .x$rdi_coef[["intercept"]])
    RDI_slo <- map_dbl(sp, ~ .x$rdi_coef[["slope"]])

    # err <- 0
    # eRDI <- 0
    # eKg <- 0

    getCostFunc <- function(Par, meshcm2, N, RDI_int, RDI_slo, x,
                            targetRDI, targetKg) {
        Dg2 <- sum(map2_dbl(meshcm2, x, `%*%`)) / sum(map_dbl(x, sum))
        hmax <- Par[1]
        k <- Par[2]
        Pc <- map(N, ~ hmax * (1:.x)^(-k))
        if( sum(map2_dbl(x, Pc, ~ sum(.x * .y))) <= 0){
            return(.9) # escape if no harvest
        }
        tmp <- map2(x, Pc, `*`)
        Dgcut2 <- sum(map2_dbl(meshcm2, tmp, `%*%`)) / sum(map_dbl(tmp, sum))
        Kg <- Dgcut2 / Dg2
        rdi_sp <- imap(x, function(x, .y, Pc, RDI_int, RDI_slo, meshcm2){
            RDI(x = x * (1 - Pc[[.y]]), RDI_int[[.y]], RDI_slo[[.y]],
                meshcm2[[.y]])
        },
        Pc = Pc, RDI_int = RDI_int, RDI_slo = RDI_slo, meshcm2 = meshcm2)
        RDI <- sum(unlist(rdi_sp))

        # cat(sprintf(
        # "error : %.5f \nKg : %.3f / Target :%.2f \nRDI : %.3f / Target :%.2f \n",
        # (Kg-targetKg)^2 + (RDI-targetRDI)^2, Kg, targetKg, RDI, targetRDI)
        # )
        # cat(sprintf(
        # "par 1 : %.5f | par 2 : %.5f \n",
        # Par[1], Par[2])
        # )
        # err <<- (Kg - targetKg)^2 + (RDI - targetRDI)^2
        # eRDI <<- RDI
        # eKg <<- Kg

        return((Kg - targetKg)^2 + (RDI - targetRDI)^2)
    }

    ParOpt <- optim(
        fn = getCostFunc,
        meshcm2 = meshcm2, N = N,
        RDI_int = RDI_int, RDI_slo = RDI_slo,
        x = x, targetRDI = targetRDI, targetKg = targetKg,

        par = c(0.1, 0.1), method = "L-BFGS-B",
        upper = c(.99, 1), lower = c(0, 1e-10)
    )$par

    # cat(sprintf(
    # "error : %.5f \nKg : %.4f / Target :%.2f | RDI : %.4f / Target :%.2f \n",
    # err, eKg, targetKg, eRDI, targetRDI)
    # )
    # cat(sprintf(
    #     "final Par 1 : %.5f Par 2 : %.5f \n",
    #     ParOpt[1], ParOpt[2])
    # )

    Pcut <-map(N, ~ ParOpt[1] * (1:.x)^(-ParOpt[2]))
    Pcut <- map2(Pcut, del, ~ delay.numeric(.x, .y))

    return(Pcut)
}

#' Compute RDI and Kg from simulations output
#'
#' @param sim simulation output table.
#' @param rdi_c rdic coefficients for all species. Named vector. If NULL
#' (default, the values are taken from matreex::rdi_coef)
#'
#' @importFrom dplyr filter select mutate group_by group_modify left_join summarise bind_rows
#' @importFrom tidyr pivot_longer pivot_wider replace_na
#'
#' @export
sim_rdikg <- function(sim, rdi_c = NULL){

    # DEV
    # forest = Jasper_for_Even_dev
    # sim = Jasper_sim_f20_dev
    # rdi_c = NULL
    #********

    sp <- unique(sim$species)

    if(is.null(rdi_c)){
        rdi_c <- matreex::rdi_coef
        if(! all(sp %in% rdi_c$species)){
            stop("A simulated species is not present in the mackage rdi coefficient. Please provide the rdi_coef argument.")
        }
    }
    rdi_c <- filter(rdi_c, species %in% sp)

    sim <- sim %>%
        filter(size > 0, ! equil) %>%
        select(- equil, -mesh) %>%
        mutate(size = (size / 10)^2)
    # Compute rdi
    rdi_sp <- sim %>%
        filter(var == "n") %>%
        left_join(rdi_c, by = "species") %>%
        group_by(species, time) %>%
        group_modify(~ data.frame(rdi = RDI(
            x = .x$value,
            RDI_int = unique(.x$intercept), RDI_slo = unique(.x$slope),
            meshcm2 = .x$size)))# %>%

    rdi_val <- rdi_sp %>%
        group_by(time) %>%
        summarise(rdi = sum(rdi), species = "All") %>%
        select(species, time, rdi) %>%
        bind_rows(rdi_sp) %>%
        pivot_longer(rdi, names_to = "var")

    # Compute Kg
    tmp <- sim %>%
        filter(var %in% c("n", "h")) %>%
        tidyr::pivot_wider(names_from = var, values_from = value)
    kg_val <- tmp %>%
        group_by(time, species) %>%
        mutate(X = n + h) %>%
        summarise(
            surfx = drop(size %*% X),
            tx = sum(X),
            surfcut = drop(size %*% h),
            tcut = sum(h),
            .groups = "drop") %>%
        group_by(time) %>%
        summarise(
            Dg2 = sum(surfx) / sum(tx),
            Dgcut2 = sum(surfcut) / sum(tcut),
            Kg = Dgcut2 / Dg2)  %>%
        replace_na(list(Kg = 0)) %>%
        mutate(species = "All") %>%
        select(species, time, Dg2, Dgcut2, Kg) %>%
        pivot_longer(cols = c("Dg2", "Dgcut2", "Kg"), names_to = "var")

    # output all
    res <- bind_rows(rdi_val, kg_val)

    return(res)
}
