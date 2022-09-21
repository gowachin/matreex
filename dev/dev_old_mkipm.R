#' Build IPM for a given species and climate
#'
#' Integrate IPM for growth and survival function at a specific climate for a
#' species on a basal area variation.
#'
#' @param species The species names to be registered in the object
#' @param climate Named vector of the environmental variables used in the fitted
#' model. Was data_plot_pred before.
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
#' @param diag_tresh TODO
#' @param midbin_tresh Number of cells external to the GL integration to
#' integrate with the mid bin method.
#' @param IsSurv Adding survival to the IPM. Set to FALSE is usefull to test for
#' eviction of the model. TRUE by default.
#'
#' @import cli
#' @import checkmate
#'
#'- @export
make_IPM <- function(
        species,
        climate,
        fit,
        mesh = c(m = 700, L = 90, U = 1500),
        BA = 0:100,
        correction = "constant",
        level = 420,
        diag_tresh = 50,
        midbin_tresh = 100,
        IsSurv = TRUE,
        verbose = FALSE){

    # Idiot Proof ####
    assertCharacter(species, len = 1)
    # TODO assert climate
    # TODO assert fit and all required climate in fit
    invar <- c(names(fit$sv$params_m),  names(fit$gr$params_m))

    assertNumeric(mesh, len = 3, lower = 1, upper = 3000)
    names <- names(mesh)
    assertCharacter(names)
    if(any(! names %in% c("m", "L", "U"))){
        stop("mesh must be consitued of m, L and U")
    }
    assertNumeric(BA, lower = 0, upper = 200)
    correction <- match.arg(correction,
                            c("none", "constant", "ceiling", "sizeExtremes"))
    assertCount(level)
    assertCount(diag_tresh)
    assertLogical(IsSurv, len = 1)

    start <- Sys.time()

    list_covs <- cbind(climate, BATOTcomp = 0)
    IPM <- vector("list", length(BA))

    # Precomput constant ####
    #  we integrate on 2 dimension along the size at t and level/3 along size at t+1
    inlevel <- floor(level/3)
    U <- mesh["U"]
    L <- mesh["L"]
    m <- mesh["m"]
    h <- (U -L)/m

    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h/2, h/2, 3) # For x integration
    weights1 <- out1$weights / sum(out1$weights) #equivalent to devided by h
    out2 <- gaussQuadInt(-h/2, h/2, inlevel) # For x1 integration
    mesh_x <- seq(L+h/2, U-h/2, length.out=m)
    N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)

    WMat <- t(build_weight_matrix(out2$weights,N_int))
    # vector for integration on dim 2
    mesh_x <- as.vector(outer(mesh_x, out1$nodes, '+'))
    # empty matrix
    e_P <- matrix(0, ncol = m, nrow = m)

    ## Functions
    svlink <- fit$sv$family$linkinv
    sig_gr <- fit$gr$sigma
    if (! IsSurv){
        P_sv <- rep(1,length(mesh_x)/3)
    }

    #Create matrix for GL integration on dimension level
    mesh_x1B <- as.vector(outer(out2$nodes,
                                seq(0,(N_int-1)*h,length.out=N_int),'+'))
    mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
    mesh_x1B <- mesh_x1B - out1$nodes[2]
    mesh_x1C <- mesh_x1B - out1$nodes[3]

    temp <- function(d_x1_x, mu, sig){
        out <- numeric(length(d_x1_x))
        sel <- d_x1_x>0
        tmp <- d_x1_x[sel]
        out[sel] <- dnorm(log(tmp), mu[sel], sig) / tmp
        return(out)
    }

    if(verbose){
        message("Launching integration loop")
        # cli_progress_bar("Integration", total = length(BA))
    }

    # Loop ####
    # for(ba in seq_along(BA)){
        ba <- 2

        if(verbose){
        #     cli_progress_update()
        }

        # Update BA and matrix
        list_covs["BATOTcomp"] <- BA[ba]
        P <- e_P
        ## Functions ####
        grFun <- exp_sizeFun(fit$gr$params_m, list_covs)
        svFun <- exp_sizeFun(fit$sv$params_m, list_covs)
        mu_gr <- grFun(mesh_x)

        ## Gauss-Legendre Int ####
        P_incr <- outer(mesh_x1A,mu_gr[1:m],'temp', sig_gr) * weights1[1] +
            outer(mesh_x1B,mu_gr[(m+1):(2*m)],'temp', sig_gr) * weights1[2] +
            outer(mesh_x1C,mu_gr[(2*m+1):(3*m)],'temp', sig_gr) * weights1[3]

        P_incr <- WMat %*% P_incr

        P <- sub_diag(P, P_incr, dist = 0)

        ## Midbin Int ####
        ## ADD mid point integration for the rest of the
        # triangular matrix (+100 points)
        P2 <- fun_mid_int_stripe(seq(L, U, length.out=m), h, grFun, sig_gr,
                                 N_ini = N_int+1, N_int = midbin_tresh, Level = 100)
        P <- sub_diag(P, P_incr, dist = N_int)
        # for (k in ((N_int+1):m)){
        #     ind_k <- k:min(k+midbin_tresh-1,m)
        #     P[ind_k,k] <- P2[1:min(m-k+1,midbin_tresh),k]
        # }
        # P[(N_int+1):(N_int+100), ] <- P[(N_int+1):(N_int+100), ] + P2

        ## Survival ####
        if (IsSurv){
            P_sv <- svlink(svFun(mesh_x))
            m <- length(mesh_x)/3
            P_sv <- P_sv[1:m] * weights1[1] +
                P_sv[(m+1):(2*m)] * weights1[2] +
                P_sv[(2*m+1):(3*m)] * weights1[3]
            P_sv <- 1 - P_sv
        }

        ## Correction ####
        if (correction ==  "none"){# Based on IPMpack
            P <- t(t(P) *P_sv)
        }

        if(correction == "constant"){# Based on IPMpack
            nvals <- colSums(P)
            P <- t((t(P)/nvals))
            P <- P * (matrix(rep(t(P_sv),m),ncol=m,nrow=m))
        }

        if (correction ==  "sizeExtremes"){# Based on IPMpack
            selectsize_t <- (N_int + 0:(m-1)) > m
            DiffNvals <- pmax(1- colSums(P), 0)
            P[m,selectsize_t] <- P[m, selectsize_t] + DiffNvals[selectsize_t]
            P <- t(t(P) *P_sv)
        }

        if (correction == "ceiling"){
            # Based on Williams et al. 2012 Ecology integral towards
            # infinity is not explicitely calculated
            P_sv_U <- svlink(svFun(U +out1$nodes))
            P_sv_U <- 1 - sum(P_sv_U * weights1)
            nvals <- colSums(P)
            P <- rbind(P,pmax((1-nvals),0))
            P <- cbind(P,c(rep(0,m),1))
            P <- t(t(P) * c(P_sv, P_sv_U))
        }

        ## Matrix and exp ####
        P <- Matrix(P,sparse=TRUE)
        IPM[[ba]] <- P
    # }


    if(verbose){
        message("Loop done.")
        tmp <- Sys.time() - start
        message("Time difference of ", format(unclass(tmp), digits = 3),
                " ", attr(tmp, "units"))
    }

    # Format ####
    names(IPM) <- BA
    res <- validate_ipm(
        new_ipm(IPM = IPM, BA = BA, mesh = seq(L, U, length.out=m), species = species,
                   climatic = "custom", compress = FALSE)
    )
    # })
    return(res)
}

# ### Build IPM for all Resamples
# make_FullIPM_iClim_resample <- function(spsel,
#                                         iClim,
#                                         save=FALSE,
#                                         status_fit_sgr=1,
#                                         status_data_plots_pred=1,
#                                         NbIPM=100,
#                                         NbModel=1,
#                                         m_size=700,
#                                         path='.'){
#     print(spsel)
#     print(iClim)
#     if (status_fit_sgr*status_data_plots_pred==1){
#         # data_plots_pred <- readRDS(file.path('output', path, paste0('data_plots_pred_',spsel,'.Rds')))
#         # fit_sgr <- readRDS(file.path('output', path, paste0('fit_sgr_all_',spsel,'.Rds')))
#         # HACK dev
#         data_plots_pred <- readRDS(file.path(path, 'output', spsel, 'plots_pred.Rds'))
#         fit_sgr <- readRDS(file.path(path, 'output', spsel, 'fit_sgr_all.Rds'))
#         M <- NULL
#         for (k in 1:NbModel){
#             print(paste0('Model : ',k))
#             M[[k]] <- make_FullIPM_iClim(iClim, data_plots_pred,
#                                          fit_sgr=fit_sgr[[k]], spsel, save=FALSE, NbIPM=NbIPM, m_size=m_size)
#         }
#         if (save==TRUE){
#             saveRDS(M,file.path('output', path, paste0('FullIPM_',spsel,'_Clim_',iClim,'.Rds')))
#             return(status_FullIPM_resample=1)
#         }else{
#             return(M)
#         }
#     }else{
#         return(status_FullIPM_resample=0)
#     }
# }
#
# ### Build IPM for 1 resample
# make_FullIPM_iClim <- function(iClim,
#                                data_plots_pred,
#                                fit_sgr,spsel,
#                                m_size = 700,
#                                level = 420,
#                                correction = "constant",
#                                diag_tresh = 50,
#                                IsSurv=TRUE,
#                                save=FALSE,
#                                NbIPM=100){
#     data_plots_pred <- filter(data_plots_pred, N==iClim)
#     data_plots_predBA <- data.frame(cbind(matrix(rep(as.matrix(data_plots_pred), each=NbIPM),
#                                                  nrow=NbIPM), 1:NbIPM))
#     names(data_plots_predBA) <- c(names(data_plots_pred), 'BATOTSP')
#     data_plots_predBA <- mutate(data_plots_predBA, BATOTcomp=BATOTSP) # Rajouter BATOTcomp ne sert a rien car pas utilise null part.
#     # # ligne alternative : (dt/8 in microsecond)
#     # data_plots_predBA <-cbind(data_plots_pred, data.frame(BATOTSP = 1:NbIPM,
#     #                                                       BATOTcomp = 1:NbIPM))
#
#     L <- 100*0.9
#     U <- 1.1* fit_sgr$maxDBH
#     h <- (U - L) / m_size
#     mesh_x <- seq(L+h/2,U-h/2,length.out=m_size)
#     N_int <- sum((mesh_x-min(mesh_x))<diag_tresh)
#     out2 <- gaussQuadInt(-h/2, h/2, floor(level/3))
#     WMat <- build_weight_matrix(out2$weights,N_int)
#     #:::::::::::::::::::::::::::::::::::::::::::
#     print("before make_IPM_GL_2_i")
#     LIPM <- lapply(1:NbIPM, make_IPM_GL_2_i,
#                    data_plots_pred=data_plots_predBA,
#                    fit_sgr=fit_sgr,
#                    spsel=spsel,
#                    m_size=m_size,
#                    level=level,
#                    correction=correction,
#                    diag_tresh=diag_tresh,
#                    WMat=WMat,
#                    IsSurv=IsSurv)
#
#     RecFun <- fun_recruitment_mean(fit_sgr$rec, data_plots_pred)
#     FullIPM <- list(LIPM=LIPM,meshpts=mesh_x,list_m=data_plots_pred,
#                     sv=fit_sgr$sv, gr=fit_sgr$gr, rec=fit_sgr$rec, RecFun=RecFun)
#
#     if (save==TRUE){
#         saveRDS(FullIPM,file=paste0('output/FullIPM_',spsel,'.Rds'))
#         return(status=1)
#     }else{
#         return(FullIPM)
#     }
# }
#
# # Build 1 BA (i) IPM
# make_IPM_GL_2_i <- function(i, data_plots_pred, fit_sgr,
#                             spsel= "Fagus sylvatica",
#                             m_size = 700,
#                             level = 420,
#                             correction = "ceiling",
#                             diag_tresh = 50,
#                             WMat,
#                             IsSurv=TRUE,
#                             midPoint=TRUE){
#     minSize <- 100*0.9
#     maxSize <- 1.1* fit_sgr$maxDBH
#     h <- (maxSize - minSize) / m_size
#     meshpts <- seq(minSize+h/2, maxSize-h/2, length.out = m_size)
#     list_m <- filter(data_plots_pred, BATOTSP==i)
#     print("before mk_P_GL_2")
#     if (!missing(WMat)){
#         P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
#                        g_res = fit_sgr$gr, s_res = fit_sgr$sv,
#                        list_covs = list_m, diag_tresh = diag_tresh,
#                        level = level, correction = correction,WMat=WMat,IsSurv=IsSurv, midPoint=midPoint)
#     }else{
#         P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
#                        g_res = fit_sgr$gr, s_res = fit_sgr$sv,
#                        list_covs = list_m, diag_tresh = diag_tresh,
#                        level = level, correction = correction,IsSurv=IsSurv, midPoint=midPoint)
#     }
#     if(correction == "ceiling") {meshpts <- c(meshpts, maxSize)}
#     return(P = P)
# }
#
#
#
# #:::::::::::::::::::::::::::::::::::::::::::::::::
# #::::::::::::::::::::::::::::::::::::::::::::::::
# # Gauss-Legendre 2 dimensions integration
# # P_int percentage of the distance to the diagonal in cells where the GL integration is applied
# mk_P_GL_2 <- function(m,
#                       L,
#                       U,
#                       g_res,
#                       s_res,
#                       list_covs,
#                       diag_tresh= 50,
#                       level=420,
#                       correction="none",
#                       WMat,
#                       IsSurv=TRUE,
#                       midPoint=TRUE){
#
#
#     if(! correction %in% c("constant", "ceiling", "sizeExtremes", "none")){
#         stop("correction must be in constant, ceiling, sizeExtremes, or none")
#     }
#     level <- floor(level/3)
#     # we integrate on 2 dimension along the size at t and level/3 along size at t+1
#     h <- (U - L) / m
#     # build weight for GL integration on the two dim
#     out1 <- gaussQuadInt(-h/2, h/2, 3) # For x integration
#     weights1 <- out1$weights / sum(out1$weights) #equivalent to devided by h
#     out2 <- gaussQuadInt(-h/2, h/2, level) # For x1 integration
#     mesh_x <- seq(L+h/2, U-h/2, length.out=m)
#     N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
#     # mesh_x_t <- mesh_x #  USELESS
#     # vector for integration on dim 2
#     mesh_x <- as.vector(outer(mesh_x, out1$nodes, '+'))
#     ### Growth
#     mu_gr <- fun_growth_mean(g_res, list_covs, mesh_x)
#     sig_gr <- g_res$sigma
#     ### Survival
#     if (IsSurv==FALSE){
#         P_sv <- rep(1,length(mesh_x)/3)
#     }else{
#         P_sv <- fun_surv_mean(s_res, list_covs, mesh_x, weights1)
#     }
#     print("before matrix")
#     ### Create matrix for GL integration on dimension level
#     mesh_x1B <- as.vector(t(outer(seq(0,(N_int-1)*h,length.out=N_int),out2$nodes,'+')))
#     mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
#     mesh_x1B <- mesh_x1B - out1$nodes[2]
#     mesh_x1C <- mesh_x1B - out1$nodes[3]
#
#     temp <- function(d_x1_x, mu, sig=sig_gr){ # TODO : default value in parent env !
#         out <- rep(0, length.out = length(d_x1_x))
#         out[d_x1_x>0] <- dnorm(log(d_x1_x[d_x1_x>0]), mu[d_x1_x>0], sig) * 1/d_x1_x[d_x1_x>0]
#         return(out)
#     }
#     P_incr <- outer(mesh_x1A,mu_gr[1:m],'temp') * weights1[1] +
#         outer(mesh_x1B,mu_gr[(m+1):(2*m)],'temp') * weights1[2] +
#         outer(mesh_x1C,mu_gr[(2*m+1):(3*m)],'temp') * weights1[3]
#     if(missing(WMat)){
#         WMat <- build_weight_matrix(out2$weights,N_int)
#     }
#     P_incr <- t(WMat) %*% P_incr
#     P <- matrix(0,ncol=m,nrow=m)
#     for (k in (1:m)){
#         ind_k <- k:min(k+N_int-1,m)
#         P[ind_k,k] <- P_incr[1:min(m-k+1,N_int),k]
#     }
#
#     print("before mid point")
#     ## ADD mid point integration for the rest of the triangular matrix (+100 points)
#     if (midPoint==TRUE){
#         P2 <- fun_mid_int(m ,L ,U , g_res, s_res ,N_ini=N_int+1, N_int=100 ,list_covs)
#         P[(N_int+1):(N_int+100), ] <- P[(N_int+1):(N_int+100), ] + P2
#     }
#     ##
#     gc()
#
#     print("before correction")
#     if(correction == "constant"){
#         # Based on IPMpack
#         nvals <- colSums(P)
#         P <- t((t(P)/nvals))
#         P <- P * (matrix(rep(t(P_sv),m),ncol=m,nrow=m))
#     }
#
#     if (correction ==  "sizeExtremes"){# Based on IPMpack
#         selectsize_t <- (N_int + 0:(m-1)) > m
#         DiffNvals <- pmax(1- colSums(P), 0)
#         P[m,selectsize_t] <- P[m, selectsize_t] + DiffNvals[selectsize_t]
#         P <- t(t(P) *P_sv)
#     }
#     if (correction ==  "none"){# Based on IPMpack
#         P <- t(t(P) *P_sv)
#     }
#     if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology integral towards infinity is not explicitely calculated
#         P_sv_U<-   fun_surv_mean(s_res, list_covs, U +out1$nodes , weights1)
#         nvals <- colSums(P)
#         P <- rbind(P,pmax((1-nvals),0))
#         P <- cbind(P,c(rep(0,m),1))
#         P <- t(t(P) * c(P_sv, P_sv_U))
#     }
#     P <- Matrix(P,sparse=TRUE)
#     return(P)
# }
#
#
# mk_P_GL_2_stripe <- function(m, L, U,
#                       g_res,
#                       s_res,
#                       list_covs,
#                       diag_tresh= 50,
#                       level=420,
#                       correction=c("none", "constant", "ceiling", "sizeExtremes"),
#                       WMat,
#                       IsSurv=TRUE, # to rm
#                       midPoint=TRUE){
#
#     profvis({
#         level=420 # HACK dev
#     correction <- match.arg(correction,
#                             c("none", "constant", "ceiling", "sizeExtremes") # HACK dev
#                             )
#
#     grFun <- exp_sizeFun(g_res$params_m, list_covs)
#     sig_gr <- g_res$sigma
#     svFun <- exp_sizeFun(s_res$params_m, list_covs)
#     svlink <- s_res$family$linkinv
#
#     level <- floor(level/3)
#     # we integrate on 2 dimension along the size at t and level/3 along size at t+1
#     h <- (U - L) / m
#     # build weight for GL integration on the two dim
#     out1 <- gaussQuadInt(-h/2, h/2, 3) # For x integration
#     weights1 <- out1$weights / sum(out1$weights) #equivalent to devided by h
#     out2 <- gaussQuadInt(-h/2, h/2, level) # For x1 integration
#     mesh_x <- seq(L+h/2, U-h/2, length.out=m)
#     N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
#
#     if(missing(WMat)){
#         WMat <- build_weight_matrix(out2$weights,N_int)
#     }
#
#     # vector for integration on dim 2
#     mesh_x <- as.vector(outer(mesh_x, out1$nodes, '+'))
#     ### Growth
#     mu_gr <- grFun(mesh_x)
#
#     print("before matrix")
#     ### Create matrix for GL integration on dimension level
#     mesh_x1B <- as.vector(t(outer(seq(0,(N_int-1)*h,length.out=N_int),out2$nodes,'+')))
#     mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
#     mesh_x1B <- mesh_x1B - out1$nodes[2]
#     mesh_x1C <- mesh_x1B - out1$nodes[3]
#
#     temp <- function(d_x1_x, mu, sig){
#         out <- numeric(length(d_x1_x))
#         sel <- d_x1_x>0
#         tmp <- d_x1_x[sel]
#         out[sel] <- dnorm(log(tmp), mu[sel], sig) / tmp
#         return(out)
#     }
#     P_incr <- outer(mesh_x1A,mu_gr[1:m],'temp', sig_gr) * weights1[1] +
#         outer(mesh_x1B,mu_gr[(m+1):(2*m)],'temp', sig_gr) * weights1[2] +
#         outer(mesh_x1C,mu_gr[(2*m+1):(3*m)],'temp', sig_gr) * weights1[3]
#
#     P_incr <- t(WMat) %*% P_incr
#     P <- matrix(0,ncol=m,nrow=m)
#     for (k in (1:m)){
#         ind_k <- k:min(k+N_int-1,m)
#         P[ind_k,k] <- P_incr[1:min(m-k+1,N_int),k]
#     }
#
#     print("before mid point")
#     ## ADD mid point integration for the rest of the triangular matrix (+100 points)
#     if (midPoint==TRUE){
#         P2 <- fun_mid_int_stripe(seq(L, U, length.out=m), h, grFun, sig_gr,
#                                  svFun, svlink,
#                                  N_ini = N_int+1, N_int = 100)
#         P[(N_int+1):(N_int+100), ] <- P[(N_int+1):(N_int+100), ] + P2
#     }
#
#     print("before correction")
#     ### Survival
#     if (IsSurv==FALSE){
#         P_sv <- rep(1,length(mesh_x)/3)
#     }else{
#         P_sv <- svlink(svFun(mesh_x))
#         m <- length(mesh_x)/3
#         P_sv <- P_sv[1:m] * weights1[1] +
#             P_sv[(m+1):(2*m)] * weights1[2] +
#             P_sv[(2*m+1):(3*m)] * weights1[3]
#         P_sv <- 1 - P_sv
#         # P_sv <- fun_surv_mean(s_res, list_covs, mesh_x, weights1)
#     }
#
#     if(correction == "constant"){
#         # Based on IPMpack
#         nvals <- colSums(P)
#         P <- t((t(P)/nvals))
#         P <- P * (matrix(rep(t(P_sv),m),ncol=m,nrow=m))
#     }
#
#     if (correction ==  "sizeExtremes"){# Based on IPMpack
#         selectsize_t <- (N_int + 0:(m-1)) > m
#         DiffNvals <- pmax(1- colSums(P), 0)
#         P[m,selectsize_t] <- P[m, selectsize_t] + DiffNvals[selectsize_t]
#         P <- t(t(P) *P_sv)
#     }
#     if (correction ==  "none"){# Based on IPMpack
#         P <- t(t(P) *P_sv)
#     }
#     if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology integral towards infinity is not explicitely calculated
#         # P_sv_U<-   fun_surv_mean(s_res, list_covs, U +out1$nodes , weights1)
#         P_sv_U <- svlink(svFun(U +out1$nodes))
#         P_sv_U <- 1 - sum(P_sv_U * weights1)
#         nvals <- colSums(P)
#         P <- rbind(P,pmax((1-nvals),0))
#         P <- cbind(P,c(rep(0,m),1))
#         P <- t(t(P) * c(P_sv, P_sv_U))
#     }
#     P <- Matrix(P,sparse=TRUE)
#     }) #, interval = 0.005)
#     return(P)
# }
#
#
# fun_mid_int <- function(m, L, U, gr, sv, N_ini, N_int, list_covs, Level=100){
#     require(dplyr)
#     require(data.table)
#     temp2 <- function(d_x1_x, mu, sig=sig_gr){
#         out <- rep(0, length.out = length(d_x1_x))
#         out[d_x1_x>0] <- dnorm(log(d_x1_x[d_x1_x>0]), mu, sig) * 1/d_x1_x[d_x1_x>0]
#         return(out)
#     }
#     mesh <- seq(L, U, length.out=m)
#     h <- (U - L)/m
#     sig_gr <- gr$sigma
#     mu_mesh <- fun_growth_mean(gr, list_covs, mesh)
#     sv_mesh <- fun_surv_mean(sv, list_covs, mesh)
#     dx1 <- seq(N_ini*h, (N_ini+N_int)*h, by=h/Level)[-1]
#     k <- 0
#     P_incr <- matrix(NA, ncol=m, nrow=N_int)
#     for (mui in mu_mesh){
#         k <- k + 1
#         svi <- sv_mesh[k]
#         g <- temp2(dx1, mui, sig_gr) * svi
#         v <- data.frame(g=g, ca=c(rep(0, Level/2), rep(1: (N_int-1), each=Level),
#                                   rep(N_int, each=Level/2)), dx=dx1)
#         v$ca[v$ca<0] <- 0
#         v <- group_by(v, ca) %>%
#             summarise(P=sum(g) * h / Level, dxm=sum(range(dx))/2) %>%
#             ungroup()
#         P_incr[, k] <- v$P[1:(N_int)]
#     }
#     return(P_incr)
# }


fun_mid_int_stripe <- function(mesh, h, gr, sig_gr, N_ini, N_int, Level=100){

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

build_weight_matrix <- function(weight,N_int){
    N_sub <- length(weight)
    ct <- ceiling(
        # matrix(rep(1:(N_sub*N_int),N_sub),ncol=N_int,nrow=N_int*N_sub) # HACK dev
        matrix(rep(1:(N_sub * N_int), N_int), ncol = N_int, nrow = N_int*N_sub)
        / N_sub
        )

    ct2 <- t(matrix(rep(1:N_int, N_int * N_sub),
                    ncol = N_int * N_sub, nrow = N_int))
    ind <- 0 * ct
    ind[ct == ct2] <- 1
    for (k in 1:N_int){
        ind[which(ind[,k]==1),k] <- weight
    }
    return(ind)
}



# Gauss-Legendre quadrature nodes and weights on interval (L,U)
gaussQuadInt <- function(L,
                         U,
                         order=7) {
    require(statmod)
    out <- gauss.quad(order) # GL is the default
    w <- out$weights; x <- out$nodes
    weights <- 0.5 * (U - L) * w
    nodes <- 0.5 * (U + L) + 0.5 * (U - L) * x
    return(list(weights = weights, nodes = nodes))
}

# fun_recruitment_mean <- function(r_rec,
#                                  list_covs,
#                                  IsFunc=TRUE){
#     params_rec <- r_rec$params_m
#     if (is.null(list_covs$BATOTNonSP)){list_covs <- mutate(list_covs,BATOTNonSP=0)}
#     params_i <- params_rec[!grepl("ntercept",names(params_rec)) &
#                                !grepl("BATOTSP",names(params_rec))] # Parameters without BA
#     params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
#     L1 <- sub('.*:','',names(params_i_inter))
#     L2 <- sub(':.*','',names(params_i_inter))
#     K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
#     params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
#     K_i <- params_rec[grepl("ntercept",names(params_rec))] +
#         sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
#         K_i_inter # K for plot i
#     ### Interaction between size and covariates has to be written as size:cov
#     params_size_inter <- params_rec[grepl('BATOTSP',names(params_rec)) &
#                                         !grepl('logBATOTSP',names(params_rec)) &
#                                         grepl(':',names(params_rec))]
#     L <- sub('.*:','',names(params_size_inter))
#     K_BA <- params_rec['BATOTSP'] + sum(unlist(list_covs[L]) * params_size_inter)
#     params_logsize_inter <- params_rec[grepl('logBATOTSP',names(params_rec)) &
#                                            grepl(':',names(params_rec))]
#     L <- sub('.*:','',names(params_logsize_inter))
#     K_logBA <- params_rec['logBATOTSP'] +
#         sum(unlist(list_covs[L]) * params_logsize_inter)
#     K_logBA[is.na(K_logBA)] <- 0
#     mu_rec <- function(BATOTSP){as.numeric(K_i + K_BA * BATOTSP + K_logBA * log(BATOTSP))}
#     if (IsFunc==FALSE){
#         mu_rec <- data.frame(K_i=K_i, K_BA=K_BA, K_logBA=K_logBA)
#     }else{
#         mu_rec <- function(BATOTSP){as.numeric(K_i + K_BA * BATOTSP + K_logBA * log(BATOTSP))}
#     }
#     return(mu_rec)
# }
#
# fun_growth_mean <- function(g_res,
#                             list_covs,
#                             mesh_x
#                             ){
#     params_gr <- g_res$params_m
#     nms <- names(params_gr)
#     params_i <- params_gr[!grepl("ntercept",nms) &
#                               !grepl("size",nms)] # Parameters without size
#     params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
#     L1 <- sub('.*:','',names(params_i_inter))
#     L2 <- sub(':.*','',names(params_i_inter))
#     K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
#     params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
#     K_i <- params_gr[grepl("ntercept",nms)] +
#         sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
#         K_i_inter # K for plot i
#     ### Interaction between size and covariates has to be written as size:cov
#     params_size_inter <- params_gr[grepl('size',nms) &
#                                        !grepl('logsize',nms) &
#                                        grepl(':',nms)]
#     L <- sub('.*:','',names(params_size_inter))
#     K_size <- params_gr['size'] + sum(unlist(list_covs[L]) * params_size_inter)
#     params_logsize_inter <- params_gr[grepl('logsize',nms) &
#                                           grepl(':',nms)]
#     L <- sub('.*:','',names(params_logsize_inter))
#     K_logsize <- params_gr['logsize'] +
#         sum(unlist(list_covs[L]) * params_logsize_inter)
#     K_logsize[is.na(K_logsize)] <- 0
#     mu_gr <- K_i + K_size * mesh_x + K_logsize * log(mesh_x)
#     return(mu_gr)
# }
#
# fun_surv_mean <- function(s_res,
#                           list_covs,
#                           mesh_x,
#                           weights1){
#     params_sv <- s_res$params_m
#     params_i <- params_sv[!grepl("ntercept",names(params_sv)) &
#                               !grepl("size",names(params_sv))] # Parameters without size
#     params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
#     L1 <- sub('.*:','',names(params_i_inter))
#     L2 <- sub(':.*','',names(params_i_inter))
#     K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
#     params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
#     K_i <- params_sv[grepl("ntercept",names(params_sv))] +
#         sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
#         K_i_inter # K for plot i
#     ### Interaction between size and covariates has to be written as size:cov
#     params_size_inter <- params_sv[grepl('size',names(params_sv)) &
#                                        !grepl('logsize',names(params_sv)) &
#                                        grepl(':',names(params_sv))]
#     L <- sub('.*:','',names(params_size_inter))
#     K_size <- params_sv['size'] + sum(unlist(list_covs[L]) * params_size_inter)
#     params_logsize_inter <- params_sv[grepl('logsize',names(params_sv)) &
#                                           grepl(':',names(params_sv))]
#     L <- sub('.*:','',names(params_logsize_inter))
#     K_logsize <- params_sv['logsize'] +
#         sum(unlist(list_covs[L]) * params_logsize_inter)
#     K_logsize[is.na(K_logsize)] <- 0 # usefull?
#     P_sv <- s_res$family$linkinv(K_i  + K_size * mesh_x + K_logsize * log(mesh_x))
#     if (!missing(weights1)){
#         m <- length(mesh_x)/3
#         P_sv <- P_sv[1:m] * weights1[1] +
#             P_sv[(m+1):(2*m)] * weights1[2] +
#             P_sv[(2*m+1):(3*m)] * weights1[3] # Avoid if no [1]
#     }
#     return(1-P_sv)
# }

