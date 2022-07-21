

R_0_SP <- function(spsel,
                   Harv=0.006,
                   correction='cut',
                   path='',
                   SurfEch=300*1e-4,
                   Ndelay=0,
                   Nclim=3,
                   IsServ=0,
                   BAintra=0.1,
                   SAVE=TRUE){
    if (IsServ==0){
        cl <- makeForkCluster(3)
        R_0 <- parLapply(cl, 1:Nclim, Compute_R_0, spsel=spsel, Harv=Harv,
                         correction=correction, path=path, SurfEch=SurfEch, Ndelay=Ndelay, BAintra=BAintra)
        stopCluster(cl)
    }else{
        R_0 <- lapply(1:Nclim, Compute_R_0, spsel=spsel, Harv=Harv,
                      correction=correction, path=path, SurfEch=SurfEch, Ndelay=Ndelay, BAintra=BAintra)
    }
    R_0 <- do.call(rbind, R_0) %>% mutate(Ndelay=Ndelay)
    if (SAVE==TRUE){
        saveRDS(file=file.path('output', path, paste0('R0_', spsel, '.Rds')), R_0)
    }else{
        return(R_0)
    }
}

#:::::::::::::::::::::::::::::::::::::::::::::::::
#:::::::::Compute Net reproductive rates :::::::::
#:::::::::::::::::::::::::::::::::::::::::::::::::
Compute_R_0 <- function(spsel,
                        iClim,
                        Harv=0.006,
                        correction='cut',
                        SurfEch=300*1e-4,
                        Ndelay=0,
                        BAintra=0.1,
                        path='',
                        NbIPM){
    FullIPM <- readRDS(file.path('output', path, paste0('FullIPM_',spsel,'_Clim_',iClim,'.Rds')))
    if (missing(NbIPM)){NbIPM <- length(FullIPM)}
    data_plots_pred <- readRDS(file.path('output', path, paste0('data_plots_pred_',spsel,'.Rds')))
    fit_sgr <- readRDS(file.path('output', path, paste0('fit_sgr_all_',spsel,'.Rds')))
    data_plots_pred <- filter(data_plots_pred, N==iClim)
    data_plots_predBA <- data.frame(cbind(matrix(rep(as.matrix(data_plots_pred), each=1), nrow=1), (BAintra)))
    names(data_plots_predBA) <- c(names(data_plots_pred), 'BATOTSP')
    data_plots_predBA <- mutate(data_plots_predBA, BATOTcomp=BATOTSP)
    print(data_plots_predBA)
    ct <- Buildct(FullIPM[[1]]$meshpts, SurfEch)
    b <- c(rep(1/2,2),rep(0,length(ct)-2))
    if (Ndelay > 0){
        b <- c(rep(0, Ndelay), b)
        ct <- t(c(rep(0, Ndelay), t(ct)))
    }
    R0 <- R30 <- R60 <- RR <- rep(NA, length(FullIPM))
    for (iM in 1:NbIPM){
        IPM0intra <- make_IPM_GL_2_i(BAintra, data_plots_predBA, fit_sgr[[iM]],
                                     spsel= spsel, m_size = 700, level = 420, correction = "constant",
                                     diag_tresh = 50, IsSurv=TRUE)
        Model <- FullIPM[[iM]]
        LIPM <- Model$LIPM
        if (max(LIPM[[30]]) > 1e6){LIPM[[30]] <- LIPM[[30]] * 1e-7;LIPM[[60]] <- LIPM[[60]] * 1e-7}
        Rec0 <- exp(Model$RecFun(BAintra)) * SurfEch/(300*1e-4)
        # TODO : en fait ici c'est les valeurs de params qui sont ranges dans un autre tableau !
        Rec30 <- exp(Model$RecFun(BAintra))*exp(Model$rec$Pval["BATOTNonSP",1] * 30 )*SurfEch/(300*1e-4)
        Rec60 <- exp(Model$RecFun(BAintra))*exp(Model$rec$Pval["BATOTNonSP",1] * 60 )*SurfEch/(300*1e-4)
        F0 <- b %*% ct * Rec0 / BAintra
        F30 <- b %*% ct * Rec30 / BAintra
        F60 <- b %*% ct * Rec60 / BAintra
        IPM0 <- as.matrix(IPM0intra * (1-Harv))
        IPM30 <- as.matrix(LIPM[[30]] * (1-Harv))
        IPM60 <- as.matrix(LIPM[[60]] * (1-Harv))
        IPM0[, length(ct)-Ndelay] <- 0
        IPM0[length(ct)-Ndelay,] <- 0
        IPM30[, length(ct)-Ndelay] <- 0
        IPM30[length(ct)-Ndelay,] <- 0
        IPM60[, length(ct)-Ndelay] <- 0
        IPM60[length(ct)-Ndelay,] <- 0
        if (Ndelay>0){
            # TODO : same as Add_Delay_Matrix without diag values !
            IPM0 <- cbind(matrix(0, nrow=length(ct), ncol=Ndelay),
                          rbind(matrix(0, ncol=length(ct)-Ndelay, nrow=Ndelay), IPM0))
            IPM30 <- cbind(matrix(0, nrow=length(ct), ncol=Ndelay),
                           rbind(matrix(0, ncol=length(ct)-Ndelay, nrow=Ndelay), IPM30))
            IPM60 <- cbind(matrix(0, nrow=length(ct), ncol=Ndelay),
                           rbind(matrix(0, ncol=length(ct)-Ndelay, nrow=Ndelay), IPM60))
        }
        R0[iM] <- eigen(F0 %*% solve(diag(length(ct)) - IPM0))$values[1]
        R30[iM] <- eigen(F30 %*% solve(diag(length(ct)) - IPM30))$values[1]
        R60[iM] <- eigen(F60 %*% solve(diag(length(ct)) - IPM60))$values[1]
    }
    OUT <- data.frame(RN0=c(R0,R30,R60), N=rep(1:length(FullIPM),3),
                      clim=iClim, sp=spsel, BA=rep(c(0,30,60), BAintra=BAintra, each=length(FullIPM)))
}


# Build 1 BA (i) IPM
make_IPM_GL_2_i <- function(i, data_plots_pred, fit_sgr,
                            spsel= "Fagus sylvatica",
                            m_size = 700,
                            level = 420,
                            correction = "ceiling",
                            diag_tresh = 50,
                            WMat,
                            IsSurv=TRUE,
                            midPoint=TRUE){
    minSize <- 100*0.9
    maxSize <- 1.1* fit_sgr$maxDBH
    h <- (maxSize - minSize) / m_size
    meshpts <- seq(minSize+h/2, maxSize-h/2, length.out = m_size)
    list_m <- filter(data_plots_pred, BATOTSP==i)
    if (!missing(WMat)){ # TODO : this makes no sense at all !
        P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
                       g_res = fit_sgr$gr, s_res = fit_sgr$sv,
                       list_covs = list_m, diag_tresh = diag_tresh,
                       level = level, correction = correction,WMat=WMat,IsSurv=IsSurv, midPoint=midPoint)
    }else{
        P <- mk_P_GL_2(m = m_size, L = minSize, U = maxSize,
                       g_res = fit_sgr$gr, s_res = fit_sgr$sv,
                       list_covs = list_m, diag_tresh = diag_tresh,
                       level = level, correction = correction,IsSurv=IsSurv, midPoint=midPoint)
    }
    if(correction == "ceiling") {meshpts <- c(meshpts, maxSize)}
    return(P = P)
}


# Gauss-Legendre 2 dimensions integration
# P_int percentage of the distance to the diagonal in cells where the GL integration is applied
mk_P_GL_2 <- function(m,
                      L,
                      U,
                      g_res,
                      s_res,
                      list_covs,
                      diag_tresh= 50,
                      level=420,
                      correction="none",
                      WMat,
                      IsSurv=TRUE,
                      midPoint=TRUE){
    if(! correction %in% c("constant", "ceiling", "sizeExtremes", "none")){
        stop("correction must be in constant, ceiling, sizeExtremes, or none")
    }
    level <- floor(level/3)
    # we integrate on 2 dimension along the size at t and level/3 along size at t+1
    h <- (U - L) / m
    # build weight for GL integration on the two dim
    out1 <- gaussQuadInt(-h/2, h/2, 3) # For x integration
    weights1 <- out1$weights / sum(out1$weights) #equivalent to devided by h
    out2 <- gaussQuadInt(-h/2, h/2, level) # For x1 integration
    mesh_x <- seq(L+h/2, U-h/2, length.out=m)
    N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
    mesh_x_t <- mesh_x
    # vector for integration on dim 2
    mesh_x <- as.vector(outer(mesh_x, out1$nodes, '+'))
    ### Growth
    mu_gr <- fun_growth_mean(g_res, list_covs, mesh_x)
    sig_gr <- g_res$sigma
    ### Survival
    if (IsSurv==FALSE){
        P_sv <- rep(1,length(mesh_x)/3)
    }else{
        P_sv <- fun_surv_mean(s_res, list_covs, mesh_x, weights1)
    }
    ### Create matrix for GL integration on dimension level
    mesh_x1B <- as.vector(t(outer(seq(0,(N_int-1)*h,length.out=N_int),out2$nodes,'+')))
    mesh_x1A <- mesh_x1B - out1$nodes[1] # to resacle the position on x
    mesh_x1B <- mesh_x1B - out1$nodes[2]
    mesh_x1C <- mesh_x1B - out1$nodes[3]

    temp <- function(d_x1_x, mu, sig=sig_gr){
        out <- rep(0, length.out = length(d_x1_x))
        out[d_x1_x>0] <- dnorm(log(d_x1_x[d_x1_x>0]), mu[d_x1_x>0], sig) * 1/d_x1_x[d_x1_x>0]
        return(out)
    }
    P_incr <- outer(mesh_x1A,mu_gr[1:m],'temp') * weights1[1] +
        outer(mesh_x1B,mu_gr[(m+1):(2*m)],'temp') * weights1[2] +
        outer(mesh_x1C,mu_gr[(2*m+1):(3*m)],'temp') * weights1[3]
    if(missing(WMat)){
        WMat <- build_weight_matrix(out2$weights,N_int)
    }
    P_incr <- t(WMat) %*% P_incr
    P <- matrix(0,ncol=m,nrow=m)
    for (k in (1:m)){
        ind_k <- k:min(k+N_int-1,m)
        P[ind_k,k] <- P_incr[1:min(m-k+1,N_int),k]
    }
    ## ADD mid point integration for the rest of the triangular matrix (+100 points)
    if (midPoint==TRUE){
        P2 <- fun_mid_int(m ,L ,U , g_res, s_res ,N_ini=N_int+1, N_int=100 ,list_covs)
        P[(N_int+1):(N_int+100), ] <- P[(N_int+1):(N_int+100), ] + P2
    }
    ##
    Test_Integration <- max(1 - colSums(P)[1:(m-N_int)])
    print(paste0('maximum error due to integration : ', Test_Integration))
    gc()
    if(correction == "constant"){
        # Based on IPMpack
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
    if (correction ==  "none"){# Based on IPMpack
        P <- t(t(P) *P_sv)
    }
    if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology integral towards infinity is not explicitely calculated
        P_sv_U<-   fun_surv_mean(s_res, list_covs, U +out1$nodes , weights1)
        nvals <- colSums(P)
        P <- rbind(P,pmax((1-nvals),0))
        P <- cbind(P,c(rep(0,m),1))
        P <- t(t(P) * c(P_sv, P_sv_U))
    }
    P <- Matrix(P,sparse=TRUE)
    return(P)
}

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

fun_growth_mean <- function(g_res,
                            list_covs,
                            mesh_x){
    params_gr <- g_res$params_m
    params_i <- params_gr[!grepl("ntercept",names(params_gr)) &
                              !grepl("size",names(params_gr))] # Parameters without size
    params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
    L1 <- sub('.*:','',names(params_i_inter))
    L2 <- sub(':.*','',names(params_i_inter))
    K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
    params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
    K_i <- params_gr[grepl("ntercept",names(params_gr))] +
        sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
        K_i_inter # K for plot i
    ### Interaction between size and covariates has to be written as size:cov
    params_size_inter <- params_gr[grepl('size',names(params_gr)) &
                                       !grepl('logsize',names(params_gr)) &
                                       grepl(':',names(params_gr))]
    L <- sub('.*:','',names(params_size_inter))
    K_size <- params_gr['size'] + sum(unlist(list_covs[L]) * params_size_inter)
    params_logsize_inter <- params_gr[grepl('logsize',names(params_gr)) &
                                          grepl(':',names(params_gr))]
    L <- sub('.*:','',names(params_logsize_inter))
    K_logsize <- params_gr['logsize'] +
        sum(unlist(list_covs[L]) * params_logsize_inter)
    K_logsize[is.na(K_logsize)] <- 0
    mu_gr <- K_i + K_size * mesh_x + K_logsize * log(mesh_x)
    return(mu_gr)
}

fun_surv_mean <- function(s_res,
                          list_covs,
                          mesh_x,
                          weights1){
    params_sv <- s_res$params_m
    params_i <- params_sv[!grepl("ntercept",names(params_sv)) &
                              !grepl("size",names(params_sv))] # Parameters without size
    params_i_inter <- params_i[grepl(":",names(params_i))] # Parameters with interaction
    L1 <- sub('.*:','',names(params_i_inter))
    L2 <- sub(':.*','',names(params_i_inter))
    K_i_inter <- sum(unlist(list_covs[L1])*unlist(list_covs[L2])*params_i_inter)
    params_i_no <- params_i[!(grepl(":",names(params_i)))] # Parameters without interaction
    K_i <- params_sv[grepl("ntercept",names(params_sv))] +
        sum(unlist(list_covs[names(params_i_no)]) * params_i_no) +
        K_i_inter # K for plot i
    ### Interaction between size and covariates has to be written as size:cov
    params_size_inter <- params_sv[grepl('size',names(params_sv)) &
                                       !grepl('logsize',names(params_sv)) &
                                       grepl(':',names(params_sv))]
    L <- sub('.*:','',names(params_size_inter))
    K_size <- params_sv['size'] + sum(unlist(list_covs[L]) * params_size_inter)
    params_logsize_inter <- params_sv[grepl('logsize',names(params_sv)) &
                                          grepl(':',names(params_sv))]
    L <- sub('.*:','',names(params_logsize_inter))
    K_logsize <- params_sv['logsize'] +
        sum(unlist(list_covs[L]) * params_logsize_inter)
    K_logsize[is.na(K_logsize)] <- 0 # usefull?
    P_sv <- s_res$family$linkinv(K_i  + K_size * mesh_x + K_logsize * log(mesh_x))
    if (!missing(weights1)){
        m <- length(mesh_x)/3
        P_sv <- P_sv[1:m] * weights1[1] +
            P_sv[(m+1):(2*m)] * weights1[2] +
            P_sv[(2*m+1):(3*m)] * weights1[3] # Avoid if no [1]
    }
    return(1-P_sv)
}

build_weight_matrix <- function(weight,N_int){
    N_sub <- length(weight)
    ct <- ceiling(matrix(rep(1:(N_sub*N_int),N_sub),ncol=N_int,nrow=N_int*N_sub)/N_sub)
    ct2 <- t(matrix(rep(1:N_int,N_int*N_sub),ncol=N_int*N_sub,nrow=N_int))
    ind <- 0 * ct
    ind[ct==ct2] <- 1
    for (k in 1:N_int){
        ind[which(ind[,k]==1),k] <- weight
    }
    return(ind)
}
