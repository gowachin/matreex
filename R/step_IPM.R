
#' Extract an IPM matrix
#'
#' Extract the IPM matrix needed during simulation for a given BA and optionally
#' climate.
#'
#' Methods are set for IPM objects and mu_gr. Only the second one require
#' to give the climate and simulation correction to apply to the matrix.
#'
#' @param x IPM or mu_gr class object.
#' @param ... Variables used depending on the class of x.
#' \describe{
#'  \item{BA}{Total basal area to get the IPM for during simulation.}
#'  \item{climate}{Climate for which the IPM is needed. Only used for mu_gr.}
#'  \item{sim_corr}{Simulation correction applied to the IPM. "cut" or "none"}
#'  \item{IsSurv}{Does this step IPM require survival. If missing or NULL,
#'  the value will be taken from x$info object}
#' }
#'
#' @export
get_step_IPM <- function(x, ...){
    UseMethod("get_step_IPM")
}

#' @method get_step_IPM ipm
#' @export
get_step_IPM.ipm <- function(x, ...){

    dots <- list(...)
    ipm <- x
    BA <- dots$BA
    climate <- dots$climate
    IsSurv <- dots$IsSurv
    if(is.null(IsSurv)){
        IsSurv <- TRUE
    }

    # Idiot Proof ####
    assertClass(ipm, "ipm")
    assertNumber(BA, lower = 0, upper = 200)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    BAsp <- ipm$BA

    low_id <- which(BAsp == max(BAsp[BAsp <= BA]))
    high_id <- which(BAsp == min(BAsp[BAsp > BA]))

    delta <- diff(BAsp[c(low_id, high_id)])
    w <- (BA - BAsp[low_id]) / delta

    res <- ipm$IPM[[low_id]] * (1  - w) + ipm$IPM[[high_id]] * w

    if(! IsSurv && as.logical(ipm$info["surv"])){
        m <- length(ipm$mesh)
        U <- ipm$mesh[[m]]
        L <- ipm$mesh[[1]]
        h <- (U - L) / m
        x_level <- ipm$int["gl1"]
        year_delta <- ipm$int["year_delta"]
        # build weight for GL integration on the two dim
        out1 <- gaussQuadInt(-h / 2, h / 2, x_level) # For x integration
        weights1 <- out1$weights / sum(out1$weights) # equivalent to divided by h
        mesh_x <- ipm$mesh
        mesh_sv <- outer(mesh_x, out1$nodes, "+")

        list_covs <-  c(climate, BATOTcomp = BA)
        svFun <- exp_sizeFun(ipm$fit$sv$params_m, list_covs)
        svlink <- ipm$fit$sv$family$linkinv
        P_sv <- svlink(svFun(mesh_sv))
        P_sv <- colSums(t(P_sv) * weights1) # alt works with integration > 3
        P_sv <- 1 - P_sv
        # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
        P_sv <- P_sv ^ year_delta

        res <- res / (matrix(rep(t(P_sv), m), ncol = m, nrow = m))
    }

    return(res)
}

#' @method get_step_IPM mu_gr
#' @export
get_step_IPM.mu_gr <- function(x, ...){

    dots <- list(...)
    mu_gr <- x
    BA <- dots$BA
    climate <- dots$climate
    sim_corr <- dots$sim_corr
    IsSurv <- dots$IsSurv

    # Idiot Proof ####
    assertClass(mu_gr, "mu_gr")
    assertNumber(BA, lower = 0, upper = 200)
    assertNumeric(climate, any.missing = FALSE)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    list_covs <- c(climate, BATOTcomp = BA)
    IPM <- vector("list", length(BA))

    # # Precomput constant ####
    m <- length(mu_gr$mesh)
    U <- mu_gr$mesh[[m]]
    L <- mu_gr$mesh[[1]]
    h <- (U - L) / m
    N_int <- nrow(mu_gr$mu_gr)
    if(is.null(IsSurv)){
        IsSurv <- as.logical(mu_gr$info["surv"])
    }
    correction <- mu_gr$info["correction"]
    x_level <- mu_gr$int["gl1"]
    year_delta <- mu_gr$int["year_delta"]
    #
    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h / 2, h / 2, x_level) # For x integration
    weights1 <- out1$weights / sum(out1$weights) # equivalent to divided by h
    mesh_x <- mu_gr$mesh
    mesh_sv <- outer(mesh_x, out1$nodes, "+")
    # empty matrix
    # P <- matrix(0, ncol = m, nrow = m)

    ## Functions ####
    ### Growth
    sig_gr <- mu_gr$fit$gr$sigma
    ### Survival
    svlink <- mu_gr$fit$sv$family$linkinv
    if(!IsSurv) {
        P_sv <- rep(1, m)
    }
    ## Functions ####
    grFun <- exp_sizeFun(mu_gr$fit$gr$params_m, list_covs)
    svFun <- exp_sizeFun(mu_gr$fit$sv$params_m, list_covs)
    mu_growth <- grFun(mesh_x)
    if (IsSurv) {
        P_sv <- svlink(svFun(mesh_sv))
        P_sv <- colSums(t(P_sv) * weights1) # alt works with integration > 3
        P_sv <- 1 - P_sv
        # NOTE reu 3/10 pour le cas ou year_delta est superieur a 1
        P_sv <- P_sv ^ year_delta
    }
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ## extract MU ####
    # IDEA round by mu step ???
    index <- findInterval(mu_growth, mu_gr$mu_tab)
    P_LG <- t(mu_gr$mu_gr[index, ])
    P <- sub_diag(matrix = NULL, P_LG, dist = 0, new = TRUE)
    # NEW allow to create new matrix in and don' duplicate on modify
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
    # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Format ####
    if(sim_corr == "cut"){
        P[, m] <- 0
        P[m, ] <- 0
    }
    ## Matrix and exp ####
    res <- Matrix(P, sparse = TRUE)

    return(res)
}
