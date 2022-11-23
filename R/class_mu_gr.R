#' Constructor of mu_gr class
#'
#' @param mu_gr Matrix of computed mesh values along a mu gradient.
#' @param BA values of BA used to get the range of mu for this species. num.
#' @param mesh vector of mesh variables. m is the number of bins, L is the
#' minimum size and U the maximum size. h will be defined in the function as
#' \eqn{h <- (U - L) / m}.
#' @param mu_tab Mu values for which the values are computed. Each row of mu_gr
#' matrix is for a given mu.
#' @param fit fit_sgr model for survival growth and recruitment for the given
#' species.
#' @param species Name of the species to run simulation on. Single char.
#' @param correction IPM correction wanted for this species. Single char.
#' @param int Integration information about levels and dimension of integration.
#'
#' @keywords internal
#' @export
new_mu_gr <- function(mu_gr, BA, mesh, mu_tab, mu_step,
                       fit, species, correction, surv = TRUE, int){

    mu_gr <- list(mu_gr = mu_gr, BA = BA, mesh = mesh,
                   mu_tab = mu_tab,
                   sv = fit$sv, gr = fit$gr, rec = fit$rec,
                   info = c(species = species, correction = correction,
                            clim_lab = "mu_gr", step = mu_step, surv = surv),
                   int = int)

    class(mu_gr) <- "mu_gr"

    return(mu_gr)
}


#' validator for species class.
#'
#' @param x species class object
#'
#' @import checkmate
#'
#' @noRd
validate_mu_gr <- function(x){

    assertClass(x, "mu_gr")
    values <- unclass(x)
    names <- attr(x, "names")

    # check names of the object ####
    assertCharacter(names)
    if(any(names != c("mu_gr", "BA", "mesh", "mu_tab", "sv", "gr", "rec",
                      "info", "int"))){
        stop(paste0("mu_gr class must be composed of elements mu_gr, BA, mesh, ",
                    "mu_tab, sv, gr, rec, info and int"))
    }

    # check all values ####
    assertMatrix(values$mu_gr, any.missing = FALSE)
    assertNumeric(values$BA)
    assertNumeric(values$mesh)
    assertNumeric(values$int)

    # check infos ####
    # assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("species", "correction", "clim_lab",
                                   "step", "surv"))){
        stop(paste0("species class must have info of elements species, correction",
                    ", clim_lab, step and surv"))
    }

    invisible(x)
}


#' Build IPM for a given species and climate
#'
#' Integrate IPM for growth and survival function at a specific climate for a
#' species on a basal area variation.
#'
#' @param species The species names to be registered in the object
#' @param climate Climate table for the species.
#' Optionnal if the species is defined in the package.
#' The format is climatic variable
#' in column and different climate by row. An example is in the package with
#' \code{cliamte_species}.
#' @param fit Fitted model for growth and survival of the species and climate
#' given. Functions will depend on size and basal area.
#' @param mesh vector of mesh variables. m is the number of bins, L is the
#' minimum size and U the maximum size. h will be defined in the function as
#' \eqn{h <- (U - L) / m}.
#' @param BA Vector of basal area to integrate on. Integrating on 0 is important
#' so use it. Integrating above 200 is absurd.
#' @param correction Correction to apply to the IPM matrix for eviction. Choices
#' constant (default), ceiling, sizeExtremes and none.
#' @param stepMu Step between each mu in the species range. This value has effect
#' simulation. Default is 1e-3. Values below 1e-5 takes longer than classical
#' matrix integration.
#' @param level Number of point to use for integration in a cell during
#' Gauss-Legendre integration. This value will be divided by 3 since size t is
#' integrated at level = 3 and size t+1 at level = level/3. single int
#' (default 420).
#' @param diag_tresh Threshold for Gauss-Legendre integration, which a distance
#' to the diagonal. Number of cell integrated is the number of cell for which
#' size t+1 - size t is inferior to this threshold. single dbl (default 50).
#' @param midbin_tresh Number of cells external to the GL integration to
#' integrate with the mid bin method.
#' @param mid_level Number of point to use for integration in a cell during
#' mid bin integration.
#' @param year_delta Number of year between 2 obersavtion when using this model.
#' default 1, single int. NOTE : value for dev usage only !
#' @param IsSurv Adding survival to the IPM. Set to FALSE is useful to test for
#' eviction of the model. TRUE by default.
#' @param verbose Print message. FALSE by defaul
#'
#'
#' @details
#' The check between climate variables and fitted variable will assert if all
#' variables in the model are provided expect variables derived from "size"
#' (size, size2, logsize), "intercept" and "BATOTcomp". An error will be
#' triggered if the climate variable is missing.
#'
#' One can desactivate each kind of integration with some treshold values.
#' A negative value in diag_tresh (ex: -1) will cancel the Gauss-Legendre
#' integration and a midbin_tresh null value (ex: 0) will cancel the midbin
#' integration.
#'
#' This is a working function to test faster integration but it integrate a
#' ba value
#'
#' @import cli
#' @import checkmate
#' @importFrom stats dnorm
#'
#' @export
make_mu_gr <- function(species,
                       fit,
                       climate = NULL,
                       mesh = c(m = 700, L = 90, U = 1500),
                       BA = 0:200,
                       correction = c("constant", "none", "ceiling", "sizeExtremes"),
                       stepMu=1e-3, # NEW
                       level = c(3, 140), # 1 not used
                       diag_tresh = 50,
                       midbin_tresh = 25,
                       mid_level = 5,
                       year_delta = 1, # NOTE reu 3/10 pour le cas ou n est superieur a 1
                       IsSurv = TRUE,
                       verbose = FALSE) {

    # Idiot Proof ####
    assertCharacter(species, len = 1)
    if(is.null(climate)){
        if(! species %in% treeforce::fit_species && missing(climate)){
            stop(paste0("This species is not listed in species for which ",
                        "treeforce package has climate."))
        }
        sp <- NULL # hack to bind value.
        climate <- subset(treeforce::climate_species,
                          sp == species,
                          select = -sp
        )
    }
    assertDataFrame(climate, nrows = 3)
    # assert fit and all required climate in fit
    nms <- unique(c(names(fit$sv$params_m), names(fit$gr$params_m)))
    nms <- unique(unlist(strsplit(nms, ":", )))
    nms <- nms[!nms %in%
                   c("", "size", "logsize", "size2", "intercept", "BATOTcomp")]
    if (!all(nms %in% names(climate))) {
        stop(paste0(
            "Missing climate variables used in fit model.\n",
            "Missing : ", paste(nms[!nms %in% names(climate)], collapse = " ")
        ))
    }
    assertNumeric(mesh, len = 3, lower = 1, upper = 3000)
    names <- names(mesh)
    assertCharacter(names)
    if (any(!names %in% c("m", "L", "U"))) {
        stop("mesh must be consitued of m, L and U")
    }
    assertIntegerish(BA, lower = 0, upper = 200)
    correction <- match.arg(correction)
    assertIntegerish(level, lower = 1, any.missing = FALSE, len = 2)
    assertNumber(diag_tresh)
    assertNumber(midbin_tresh)
    assertCount(mid_level)
    assertCount(year_delta)
    assertLogical(IsSurv, len = 1)
    assertLogical(verbose, len = 1)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    start <- Sys.time()
    # Precomput constant ####
    #  we integrate on 2 dimension along the size at t and level/3 along size at t+1
    U <- mesh["U"]
    L <- mesh["L"]
    m <- unname(mesh["m"])
    h <- (U - L) / m

    N_int <- ceiling(diag_tresh / h)
    # minor midbin_tresh if mesh is shorter.
    if(midbin_tresh + N_int > m){
        midbin_tresh <- m - N_int
    }

    # get mu range
    mu_range <- getRangemu(climate = climate, fit = fit, BA = BA,
                           mesh = seq(L, U, by = h))
    if (verbose) {
        message("Mu range done")
    }
    mu_mesh <- seq(0, (N_int + midbin_tresh)*h, by=h)
    mu_tab <- seq(floor(mu_range["min"]), ceiling(mu_range["max"]), by=stepMu)

    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h / 2, h / 2, level[[1]]) # For x integration
    weights1 <- out1$weights / sum(out1$weights) # equivalent to divided by h
    out2 <- gaussQuadInt(-h / 2, h / 2, level[[2]]) # For x1 integration
    if(N_int > 0){
        WMat <- t(build_weight_matrix(out2$weights, N_int))
    }
    # Create matrix for GL integration on dimension level
    mesh_x1B <- as.vector(outer(
        out2$nodes, seq(0, (N_int - 1) * h, length.out = N_int), "+"
    ))
    mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
    mesh_x1B <- mesh_x1B - out1$nodes[2]
    mesh_x1C <- mesh_x1B - out1$nodes[3]

    meshMid <- seq(mu_mesh[(N_int + (N_int > 0))] - h/2,
                   mu_mesh[N_int + (N_int > 0) + midbin_tresh -1 ] + h/2,
                   by=h/mid_level)[-1]
    ca <- factor(rep(1:midbin_tresh, each = mid_level))
    ca <- .Internal(split(1:length(meshMid), ca))
    # empty matrix
    mu_gr <- matrix(0, ncol = N_int + (N_int > 0) + midbin_tresh - 1,
                     nrow = length(mu_tab))
    # Functions ####
    ### Growth
    svlink <- fit$sv$family$linkinv
    sig_gr <- fit$gr$sigma

    # function used in gauss legendre integration
    temp <- function(d_x1_x, mu, sig, year_delta = 1){
        out <- numeric(length(d_x1_x))
        sel <- d_x1_x > 0
        tmp <- d_x1_x[sel]
        # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
        out[sel] <- dnorm(log(tmp/ year_delta), mu[sel], sig) / tmp * year_delta
        return(out)
    }

    if (verbose) {
        message("Launching mu computation loop")
        if(N_int == 0){
            message("GL integration won't occur because of negative treshold")
        } else {
            message("GL integration occur on ", N_int, " cells")
        }
        if(midbin_tresh == 0){
            message("midbin integration won't occur because of treshold at 0")
        } else {
            message("midbin integration occur on ", midbin_tresh, " cells")
        }
        cli_progress_bar("Integration", total = length(mu_tab))
    }

    # Loop ####
    for (mu in seq_along(mu_tab)) {
        if (verbose) {
            cli_progress_update()
        }

        ## Gauss-Legendre ####
        #     tt <- vapply(meshQuad2, FUN = temp, mu = mu_tab[mu], sig = sig_gr, # old way
        #                  year_delta = year_delta, numeric(1))
        #     tt <- structure(tt, dim = dim(meshQuad2))
        #     ValQuad <- tt %*% out2$weights
        if(N_int > 0){
            P_LG <- outer(mesh_x1A, mu_tab[mu], "temp", sig_gr, year_delta) * weights1[1] +
                outer(mesh_x1B, mu_tab[mu], "temp", sig_gr, year_delta) * weights1[2] +
                outer(mesh_x1C, mu_tab[mu], "temp", sig_gr, year_delta) * weights1[3]

            ValQuad <- WMat %*% P_LG
        }
        ## Midbin
        tt <- sapply(meshMid, temp, mu = mu, sig = sig_gr) * h / mid_level
        ValMid <- unlist(lapply(ca, function(i) sum(tt[i])))

        mu_gr[mu, ] <- c(ValQuad, ValMid)
    }


    if (verbose) {
        message("Loop done.")
        tmp <- Sys.time() - start
        message(
            "Time difference of ", format(unclass(tmp), digits = 3),
            " ", attr(tmp, "units")
        )
    }

    # res <- validate_mu_gr( # TODO validate_mu_gr
    res <- new_mu_gr(
        mu_gr = mu_gr, BA = BA,
        mesh = seq(L + h / 2, U - h / 2, length.out = m),
        mu_tab = mu_tab, mu_step = stepMu,
        fit = fit, species = species, correction = correction,
        surv = IsSurv,
        int = c(gl1 = level[[1]], gl2 = level[[2]],
                gl_tresh = N_int, mb_tresh = midbin_tresh,
                mid_level = mid_level, year_delta = year_delta)
    )
    # )
    return(res)
}


#' Get mu range for a species
#'
#' @param climate climate table for the species.
#' @param fit Fitted model for growth and survival of the species and climate
#' given. Functions will depend on size and basal area.
#' @param mesh vector of mesh variables. m is the number of bins, L is the
#' minimum size and U the maximum size. h will be defined in the function as
#' \eqn{h <- (U - L) / m}.
#' @param BA Vector of basal area to integrate on. Integrating on 0 is important
#' so use it. Integrating above 200 is absurd.
#'
#' @keywords internal
#' @export
getRangemu <- function(climate,
                       fit,
                       BA = 0:200,
                       mesh = seq(90, 900, by = 10)) {

    assertIntegerish(BA, lower = 0, upper = 200)
    assertNumeric(mesh, lower = 0)

    climate_species <- climate
    N <- NULL # hack to bind global value

    fres <- data.frame(min = 1:3, max = 1:3)
    for (Nc in 1:3) {
        climate <- subset(climate_species, N == Nc, select = -N)
        climate <- drop(as.matrix(climate)) # we need it as a vector.
        list_covs <- c(climate, BATOTcomp = 0)
        res <- matrix(
            ncol = 2, nrow = length(BA),
            dimnames = list(NULL, c("min", "max"))
        )

        for (iBA in seq_along(BA)) {
            list_covs["BATOTcomp"] <- BA[iBA]

            grFun <- exp_sizeFun(fit$gr$params_m, list_covs)
            mu <- grFun(mesh)

            res[iBA, ] <- range(mu)
        }
        fres[Nc, ] <- c(min(res[, "min"]), max(res[, "max"]))
    }
    range <- c(min = min(fres$min), max = max(fres$max), sig = fit$gr$sigma)

    return(range)
}
