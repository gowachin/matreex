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
    int <- c(gl1 = level[[1]], gl2 = level[[2]], gl_tresh = N_int, gl_min = 0,
             mb_tresh = midbin_tresh, mid_level = mid_level, mb_max = 0,
             year_delta = year_delta, max_error = 0)
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
                purrr::reduce(`+`)

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

    int["max_error"] <- max(MaxError)
    int["mb_max"] <- max(MB_max)
    int["gl_min"] <- min(GL_min)

    if(int["mb_max"] > 1e-2){
        warning(paste(
            "At least one mid_bin integration has values above 1e-2.",
            "This is linked with insufficient Gauss-Legendre",
            "integration treshold. value :", diag_tresh, "mm."))
    }

    if(correction == "ceiling"){
        out_mesh <- c(out_mesh, U)
    }

    names(IPM) <- BA
    res <- validate_ipm(
        new_ipm(
            IPM = IPM, BA = BA, mesh = seq(L + h / 2, U - h / 2, length.out = m),
            climatic = climate, clim_lab = clim_lab, fit = fit,
            species = species, correction = correction,
            compress = FALSE, int = int, survival = IsSurv
        )
    )
    return(res)
}

