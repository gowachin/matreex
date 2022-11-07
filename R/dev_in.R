#' Build IPM for a given species and climate
#'
#' Integrate IPM for growth and survival function at a specific climate for a
#' species on a basal area variation.
#'
#' @inheritParams make_IPM
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
#' This is a working function to test faster integration.
#'
#' @import cli
#' @import checkmate
#' @importFrom stats dnorm
#'
#' @keywords internal
dev_make_IPM <- function(species,
                     climate,
                     clim_lab,
                     fit,
                     mesh = c(m = 700, L = 90, U = 1500),
                     BA = 0:100,
                     correction = c("constant", "none", "ceiling", "sizeExtremes"),
                     level = c(3, 140),
                     diag_tresh = 50,
                     midbin_tresh = 25,
                     mid_level = 5,
                     year_delta = 1, # NOTE reu 3/10 pour le cas ou n est superieur a 1
                     IsSurv = TRUE,
                     verbose = FALSE) {

    # Idiot Proof ####
    assertCharacter(species, len = 1)
    assertNumeric(climate, any.missing = FALSE)
    assertCharacter(clim_lab, len = 1)
    # assert fit and all required climate in fit
    nms <- unique(names(fit$sv$params_m), names(fit$gr$params_m))
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
    assertNumeric(BA, lower = 0, upper = 200)
    correction <- match.arg(correction)
    assertIntegerish(level, lower = 1, any.missing = FALSE, len = 2)
    assertNumeric(diag_tresh,any.missing = FALSE, len = 1)
    assertLogical(IsSurv, len = 1)
    assertCount(year_delta)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    start <- Sys.time()

    list_covs <- c(climate, BATOTcomp = 0)
    IPM <- vector("list", length(BA))

    # Precomput constant ####
    #  we integrate on 2 dimension along the size at t and level/3 along size at t+1
    U <- mesh["U"]
    L <- mesh["L"]
    m <- unname(mesh["m"])
    h <- (U - L) / m
    int_log_err <- 1:(m/2)

    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h / 2, h / 2, level[[1]]) # For x integration
    weights1 <- out1$weights / sum(out1$weights) # equivalent to divided by h
    out2 <- gaussQuadInt(-h / 2, h / 2, level[[2]]) # For x1 integration
    mesh_x <- seq(L + h / 2, U - h / 2, length.out = m)
    N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
    # N_int <- ceiling(diag_tresh / h) # IDEA equivalent and quick
    # minor midbin_tresh if mesh is shorter.
    if(midbin_tresh + N_int > m){
        midbin_tresh <- m - N_int
    }
    if(N_int > 0){
        WMat <- t(build_weight_matrix(out2$weights, N_int))
    }
    # vector for integration on dim 2
    # mesh_x <- as.vector(outer(mesh_x, out1$nodes, "+"))
    mesh_x <- outer(mesh_x, out1$nodes, "+") # HACK testing stuff

    # empty matrix
    e_P <- matrix(0, ncol = m, nrow = m)

    ## Functions ####
    ### Growth
    svlink <- fit$sv$family$linkinv
    sig_gr <- fit$gr$sigma
    ### Survival
    if(!IsSurv) {
        P_sv <- rep(1, m)
    }

    # Create matrix for GL integration on dimension level
    # mesh_x1B <- outer(out2$nodes, 0:(N_int - 1) * h, "+") # HACK faster and clean
    mesh_x1B <- as.vector(outer(
        out2$nodes, seq(0, (N_int - 1) * h, length.out = N_int), "+"
    ))

    mesh_x1 <- purrr::map(out1$nodes,  ~ mesh_x1B - .x)
    # mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
    # mesh_x1B <- mesh_x1B - out1$nodes[2] # IDEA why - here !!!
    # mesh_x1C <- mesh_x1B - out1$nodes[3]

    # function used in gauss legendre integration
    temp <- function(d_x1_x, mu, sig, year_delta) {
        out <- numeric(length(d_x1_x))
        sel <- d_x1_x > 0
        tmp <- d_x1_x[sel]
        # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
        out[sel] <- dnorm(log(tmp/ year_delta), mu[sel], sig) / tmp * year_delta
        return(out)
    }

    if (verbose) {
        message("Launching integration loop")
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
        cli_progress_bar("Integration", total = length(BA))
    }

    ## loggin ####
    int_log <- c(year_delta = year_delta, MaxError = 0,
                 GL_Nint = N_int, GL_level1 = level[[1]],  GL_level2 = level[[2]], GL_min = 0,
                 MB_Nint = midbin_tresh, MB_level = mid_level, MB_max = 0)
    MaxError <- numeric(length(BA))
    GL_min <- numeric(length(BA))
    MB_max <- numeric(length(BA))
    # Loop ####
    for (ba in seq_along(BA)) {
        if (verbose) {
            cli_progress_update()
        }

        # Update BA and matrix
        list_covs["BATOTcomp"] <- BA[ba]
        P <- e_P
        ## Functions ####
        grFun <- exp_sizeFun(fit$gr$params_m, list_covs)
        svFun <- exp_sizeFun(fit$sv$params_m, list_covs)
        # # HACK testing stuff
        # mu_gr <- as.vector(grFun(mesh_x)) # same line
        mu_gr <- grFun(mesh_x)
        if (IsSurv) {
            P_sv <- svlink(svFun(mesh_x))
            P_sv <- purrr::imap(weights1, function(.x, .y, P_sv){
                P_sv[, .y] * .x
                }, P_sv = P_sv) %>%
                reduce(`+`)

            # P_sv <- P_sv[1:m] * weights1[1] +
            #     P_sv[(m + 1):(2 * m)] * weights1[2] +
            #     P_sv[(2 * m + 1):(3 * m)] * weights1[3]
            # # HACK testing stuff
            # P_sv <- svlink(svFun(mesh_x))
            # P_sv <- P_sv[, 1] * weights1[1] + # works quickier with matrix
            #     P_sv[, 2] * weights1[2] +
            #     P_sv[, 3] * weights1[3]
            # P_sv <- colSums(t(P_sv) * weights1) # alt works with integration > 3
            P_sv <- 1 - P_sv
            # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
            P_sv <- P_sv ^ year_delta
        }

        # Integration ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ## Gauss-Legendre Int ####
        if(N_int > 0){
            # old_P_LG <- outer(mesh_x1A, mu_gr[1:m], "temp", sig_gr, year_delta) * weights1[1] +
            #     outer(mesh_x1B, mu_gr[(m + 1):(2 * m)], "temp", sig_gr, year_delta) * weights1[2] +
            #     outer(mesh_x1C, mu_gr[(2 * m + 1):(3 * m)], "temp", sig_gr, year_delta) * weights1[3]

            P_LG <- purrr::imap(weights1,
                               function(.x, .y, mesh, mu, sig_gr, year_delta)
                               {
                                  outer(mesh[[.y]], mu[, .y], "temp", sig_gr, year_delta) * .x
                               }, mesh = mesh_x1, mu = mu_gr,
                               sig_gr = sig_gr, year_delta = year_delta) %>%
                reduce(`+`)

            P_LG <- WMat %*% P_LG
            P <- sub_diag(P, P_LG, dist = 0)

            GL_min[ba] <- min(P_LG) # loggin
        }
        ## Midbin Int ####
        ## ADD mid point integration for the rest of the triangular matrix
        if(midbin_tresh > 0){
            P_midint <- fun_mid_int(
                seq(L, U, length.out = m), h, grFun, sig_gr,
                # if no GL, N_int mesh cell is not integrated
                N_ini = N_int + (N_int > 0), N_int = midbin_tresh,
                # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
                Level = mid_level, year_delta = year_delta
            )
            P <- sub_diag(P, P_midint, dist = N_int)

            MB_max[ba] <- max(P_midint) # loggin
        }
        MaxError[ba] <- max(1 - colSums(P)[int_log_err]) # loggin
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ## Correction ####
        if (correction == "none") { # Based on IPMpack
            P <- t(t(P) * P_sv)
        } else if (correction == "constant") { # Based on IPMpack
            nvals <- colSums(P)
            P <- t((t(P) / nvals))
            P <- P * (matrix(rep(t(P_sv), m), ncol = m, nrow = m))
        } else if (correction == "sizeExtremes") { # Based on IPMpack
            selectsize_t <- (N_int + 0:(m - 1)) > m
            DiffNvals <- pmax(1 - colSums(P), 0)
            P[m, selectsize_t] <- P[m, selectsize_t] + DiffNvals[selectsize_t]
            P <- t(t(P) * P_sv)
        } else if (correction == "ceiling") {
            # Based on Williams et al. 2012 Ecology integral towards
            # infinity is not explicitely calculated
            P_sv_U <- svlink(svFun(U + out1$nodes))
            P_sv_U <- 1 - sum(P_sv_U * weights1)
            nvals <- colSums(P)
            P <- rbind(P, pmax((1 - nvals), 0))
            P <- cbind(P, c(rep(0, m), 1))
            P <- t(t(P) * c(P_sv, P_sv_U))
        }
        # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        ## Matrix and exp ####
        P <- Matrix(P, sparse = TRUE)
        IPM[[ba]] <- P
    }


    if (verbose) {
        message("Loop done.")
        tmp <- Sys.time() - start
        message(
            "Time difference of ", format(unclass(tmp), digits = 3),
            " ", attr(tmp, "units")
        )
    }
    # Format ####

    int_log["MaxError"] <- max(MaxError)
    int_log["MB_max"] <- max(MB_max)
    int_log["GL_min"] <- min(GL_min)

    if(int_log["MB_max"] > 1e-2){
        warning(paste(
            "At least one mid_bin integration has values above 1e-2.",
            "This is linked with insufficient Gauss-Legendre",
            "integration treshold. value :", diag_tresh, "mm."))
    }

    names(IPM) <- BA
    res <- validate_ipm(
        new_ipm(
            IPM = IPM, BA = BA, mesh = seq(L, U, length.out = m),
            climatic = climate, clim_lab = clim_lab, rec_params = fit$rec$params_m,
            species = species, compress = FALSE, int_log = int_log
        )
    )
    return(res)
}


#' Build IPM for a given species and climate
#'
#' Integrate IPM for growth and survival function at a specific climate for a
#' species on a basal area variation.
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
#' @keywords internal
make_mutrix <- function(species,
                        fit,
                        mesh = c(m = 700, L = 90, U = 1500),
                        BA = 0:200,
                        correction = c("constant", "none", "ceiling", "sizeExtremes"),
                        stepMu=1e-3, # NEW
                        level = c(3, 140), # 1 not used
                        diag_tresh = 50,
                        midbin_tresh = 25,
                        mid_level = 5,
                        year_delta = 1, # NOTE reu 3/10 pour le cas ou n est superieur a 1
                        verbose = FALSE) {

    # Idiot Proof ####
    assertCharacter(species, len = 1)

    if(! species %in% treeforce::fit_species){
        stop(paste0("This species is not listed in species for which ",
                    "treeforce package has climate."))
    }
    climate <- subset(treeforce::climate_species,
                              sp == species,
                              select = -sp
    )
    assertDataFrame(climate, nrows = 3)
    # assert fit and all required climate in fit
    nms <- unique(names(fit$sv$params_m), names(fit$gr$params_m))
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
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    start <- Sys.time()

    list_covs <- c(climate, BATOTcomp = 0)
    IPM <- vector("list", length(BA))

    # Precomput constant ####
    #  we integrate on 2 dimension along the size at t and level/3 along size at t+1
    U <- mesh["U"]
    L <- mesh["L"]
    m <- unname(mesh["m"])
    h <- (U - L) / m
    int_log_err <- 1:(m/2)

    # N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
    N_int <- ceiling(diag_tresh / h) # IDEA equivalent and quick
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

    meshQuad <- mu_mesh[1:N_int]
    outQuad <- gaussQuadInt(-h/2, h/2, level[[2]])
    meshQuad2 <- outer(meshQuad, outQuad$nodes, '+')

    meshMid <- seq(mu_mesh[(N_int+1+(N_int > 0))] - h/2,
                   mu_mesh[N_int + midbin_tresh +(N_int > 0)] + h/2,
                   by=h/mid_level)[-1]
    ca <- factor(rep(1:N_int, each = mid_level))
    ca <- .Internal(split(1:length(meshMid), ca))

    # empty matrix
    Mutrix <- matrix(0, ncol = N_int + midbin_tresh, nrow = length(mu_tab))

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
        tt <- vapply(meshQuad2, FUN = temp, mu = mu_tab[mu], sig = sig_gr,
                     year_delta = year_delta, numeric(1))
        tt <- structure(tt, dim = dim(meshQuad2))
        ValQuad <- tt %*% outQuad$weights
        ## Midbin
        tt <- sapply(meshMid, temp, mu = mu, sig = sig_gr) * h / mid_level
        ValMid <- unlist(lapply(ca, function(i) sum(tt[i])))

        Mutrix[mu, ] <- c(ValQuad, ValMid)
    }


    if (verbose) {
        message("Loop done.")
        tmp <- Sys.time() - start
        message(
            "Time difference of ", format(unclass(tmp), digits = 3),
            " ", attr(tmp, "units")
        )
    }
    # res <- validate_mutrix( # TODO validate_mutrix
    res <- new_mutrix(
            mutrix = Mutrix, BA = BA, mesh = mu_mesh, mu_tab = mu_tab,
            fit = fit, species = species, correction = correction
        )
    # )

    return(res)
}

#' Constructor of mutrix class
#'
#' @param mutrix Matrix of computed mesh values along a mu gradient.
#' @param BA values of BA used to get the range of mu for this species. num.
#' @param mesh mesh of Delta size t-1 and size t. This mesh is lower
#' than the full mesh because we do not inegrate the full matrix. num
#' @param mu_tab Mu values for which the values are computed. Each row of mutrix
#' matrix is for a given mu.
#' @param fit fit_sgr model for survival growth and recruitment for the given
#' species.
#' @param species Name of the species to run simulation on. Single char.
#' @param correction IPM correction wanted for this species. Single char.
#'
#' @export
new_mutrix <- function(mutrix, BA, mesh, mu_tab, fit, species, correction){

    mutrix <- list(mutrix = mutrix, BA = BA, mesh = mesh, mu_tab = mu_tab,
                   fit = fit,
                   info = c(species = species, correction = correction))

    class(mutrix) <- "mutrix"

    return(mutrix)
}


#' Get mu range for a species
#'
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

  fres <- data.frame(min = 1:3, max = 1:3)
  for (Nc in 1:3) {
    climate <- subset(climate_species, N == Nc, select = -N)
    climate <- drop(as.matrix(climate)) # we need it as a vector.
    list_covs <- c(climate, BATOTcomp = 0)
    # browser()
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
