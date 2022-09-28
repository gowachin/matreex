#' Build IPM for a given species and climate
#'
#' `r lifecycle::badge("experimental")`
#'
#' Integrate IPM for growth and survival function at a specific climate for a
#' species on a basal area variation.
#'
#' @param species The species names to be registered in the object
#' @param climate Named vector of the environmental variables used in the fitted
#' model. Was data_plot_pred before.
#' @param clim_lab Label for climatic used. This values will be matched when
#' simulating multiple species together.
#' @param fit Fitted model for growth and survival of the species and climate
#' given. Functions will depend on size and basal area.
#' @param mesh vector of mesh variables. m is the number of bins, L is the
#' minimum size and U the maximum size. h will be defined in the function as
#' \eqn{h <- (U - L) / m}.
#' @param BA Vector of basal area to integrate on. Integrating on 0 is important
#' so use it. Integrating above 200 is absurd.
#' @param correction Correction to apply to the IPM matrix for eviction. Choices
#' constant (default), ceiling, sizeExtremes and none.
#' @param level TODO
#' @param diag_tresh Treshold
#' @param midbin_tresh Number of cells external to the GL integration to
#' integrate with the mid bin method.
#' @param IsSurv Adding survival to the IPM. Set to FALSE is usefull to test for
#' eviction of the model. TRUE by default.
#' @param verbose Print message. FALSE by default
#'
#' @details
#' The check between climate variables and fitted variable will assert if all
#' variables in the model are provided expect variables derived from "size"
#' (size, size2, logsize), "intercept" and "BATOTcomp". An error will be
#' triggered if the climate variable is missing.
#'
#' @import cli
#' @import checkmate
#' @importFrom stats dnorm
#'
#' @details
#' One can desactivate each kind of integration with some treshold values.
#' A negative value in diag_tresh (ex: -1) will cancel the Gauss-Legendre
#' integration and a midbin_tresh null value (ex: 0) will cancel the midbin
#' integration.
#'
#' @export
make_IPM <- function(species,
                     climate,
                     clim_lab,
                     fit,
                     mesh = c(m = 700, L = 90, U = 1500),
                     BA = 0:100,
                     correction = c("constant", "none", "ceiling", "sizeExtremes"),
                     level = 420,
                     diag_tresh = 50,
                     midbin_tresh = 100,
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
    assertCount(level)
    assertNumeric(diag_tresh,any.missing = FALSE, len = 1)
    assertLogical(IsSurv, len = 1)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    start <- Sys.time()

    list_covs <- c(climate, BATOTcomp = 0)
    IPM <- vector("list", length(BA))

    # Precomput constant ####
    #  we integrate on 2 dimension along the size at t and level/3 along size at t+1
    inlevel <- floor(level / 3)
    U <- mesh["U"]
    L <- mesh["L"]
    m <- mesh["m"]
    h <- (U - L) / m

    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h / 2, h / 2, 3) # For x integration
    weights1 <- out1$weights / sum(out1$weights) # equivalent to divided by h
    out2 <- gaussQuadInt(-h / 2, h / 2, inlevel) # For x1 integration
    mesh_x <- seq(L + h / 2, U - h / 2, length.out = m)
    N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
    # minor midbin_tresh is mesh is shorter.
    if(midbin_tresh + N_int > m){
        midbin_tresh <- m - N_int
    }
    if(N_int > 0){
        WMat <- t(build_weight_matrix(out2$weights, N_int))
    }
    # vector for integration on dim 2
    mesh_x <- as.vector(outer(mesh_x, out1$nodes, "+"))
    # empty matrix
    e_P <- matrix(0, ncol = m, nrow = m)

    ## Functions ####
    ### Growth
    svlink <- fit$sv$family$linkinv
    sig_gr <- fit$gr$sigma
    ### Survival
    if(!IsSurv) {
        P_sv <- rep(1, length(mesh_x) / 3)
    }

    # Create matrix for GL integration on dimension level
    mesh_x1B <- as.vector(outer(
        out2$nodes,
        seq(0, (N_int - 1) * h, length.out = N_int),
        "+"
    ))
    mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
    mesh_x1B <- mesh_x1B - out1$nodes[2]
    mesh_x1C <- mesh_x1B - out1$nodes[3]

    # function used in gauss legendre integration
    temp <- function(d_x1_x, mu, sig) {
        out <- numeric(length(d_x1_x))
        sel <- d_x1_x > 0
        tmp <- d_x1_x[sel]
        out[sel] <- dnorm(log(tmp), mu[sel], sig) / tmp
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
        mu_gr <- grFun(mesh_x)
        if (IsSurv) {
            P_sv <- svlink(svFun(mesh_x))
            m3 <- length(mesh_x) / 3
            P_sv <- P_sv[1:m3] * weights1[1] +
                P_sv[(m3 + 1):(2 * m3)] * weights1[2] +
                P_sv[(2 * m3 + 1):(3 * m3)] * weights1[3]
            P_sv <- 1 - P_sv
        }

        # Integration ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        ## Gauss-Legendre Int ####
        if(N_int > 0){
            P_LG <- outer(mesh_x1A, mu_gr[1:m], "temp", sig_gr) * weights1[1] +
                outer(mesh_x1B, mu_gr[(m + 1):(2 * m)], "temp", sig_gr) * weights1[2] +
                outer(mesh_x1C, mu_gr[(2 * m + 1):(3 * m)], "temp", sig_gr) * weights1[3]

            P_LG <- WMat %*% P_LG
            P <- sub_diag(P, P_LG, dist = 0)
        }
        ## Midbin Int ####
        ## ADD mid point integration for the rest of the triangular matrix
        if(midbin_tresh > 0){
            P_midint <- fun_mid_int(
                seq(L, U, length.out = m), h, grFun, sig_gr,
                N_ini = N_int + 1, N_int = midbin_tresh, Level = 100
            )
            P <- sub_diag(P, P_midint, dist = N_int)
        }
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
    names(IPM) <- BA
    res <- validate_ipm(
        new_ipm(
            IPM = IPM, BA = BA, mesh = seq(L, U, length.out = m),
            climatic = climate, clim_lab = clim_lab,
            species = species, compress = FALSE
        )
    )
    return(res)
}


#' Function to integrate growth  with midbin function.
#'
#' @param mesh Mesh to integrate on, a vector of size m between U and L. dbl.
#' @param h mesh cell size. single dbl.
#' @param gr Growth function from the fitted model. Argument is size vector.
#' @param sig_gr Sigma of the fitted growth function.
#' @param N_ini Distance to diagonal to start the integration. single dbl.
#' @param N_int Number of cells to integrate on. single dbl.
#' @param Level TODO
#'
#' @importFrom stats dnorm
#'
#' @noRd
fun_mid_int <- function(mesh, h, gr, sig_gr, N_ini, N_int, Level=100){

    mu_mesh <- gr(mesh)
    dx1 <- seq(N_ini*h, (N_ini+N_int)*h, by=h/Level)[-1]
    inf <- dx1 > 0
    dx1i <- dx1[inf]
    ldx1 <- log(dx1i)

    rep0 <- numeric(length = length(dx1))

    ca <- factor(rep(0:N_int, c(Level/2, rep(Level, N_int-1), Level/2)))
    ca <- .Internal(split(1:length(dx1), ca))

    P_incr <- matrix(NA_real_, ncol= length(mesh), nrow=N_int)

    for (k in seq_along(mu_mesh)){
        out <- rep0 # 0.004ms
        out[inf] <- dnorm(ldx1, mu_mesh[k], sig_gr) / dx1i   # 220ms
        # C code not speedable or loose precision
        # https://stackoverflow.com/questions/27425658/speed-up-dnorm-function
        g <- out * h / Level
        res <- unlist(lapply(ca, function(i) sum(g[i])))
        P_incr[, k] <- res[1:N_int] # 380ms
    }

    return(P_incr)
}



#' Quadratique integration
#'
#' @param L Minimum mesh size. single dbl.
#' @param U Maximum mesh size. single dbl.
#' @param order TODO, default is 7. int.
#'
#' @importFrom statmod gauss.quad
#' @import checkmate
#' @noRd
gaussQuadInt <- function(L,
                         U,
                         order = 7) {

    assertNumber(L)
    assertNumber(U, lower = L)
    assertCount(order)


    out <- gauss.quad(order) # GL is the default
    w <- out$weights
    x <- out$nodes
    weights <- 0.5 * (U - L) * w
    nodes <- 0.5 * (U + L) + 0.5 * (U - L) * x
    return(list(weights = weights, nodes = nodes))
}


#' Build matrix of weight for integration.
#'
#' @param weight TODO
#' @param N_int
#'
#' @import checkmate
#'
#' @noRd
build_weight_matrix <- function(weight, N_int) {

    assertNumeric(weight)
    assertCount(N_int)

    N_sub <- length(weight)
    ct <- ceiling(
        matrix(rep(1:(N_sub * N_int), N_int), ncol = N_int, nrow = N_int * N_sub)
        / N_sub
    )

    ct2 <- t(matrix(rep(1:N_int, N_int * N_sub),
                    ncol = N_int * N_sub, nrow = N_int
    ))
    ind <- 0 * ct
    ind[ct == ct2] <- 1
    for (k in 1:N_int) {
        ind[which(ind[, k] == 1), k] <- weight
    }

    return(ind)
}