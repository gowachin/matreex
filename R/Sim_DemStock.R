
SaveSimDemNew <- function(spsel,
                          SurfEch=1,
                          Ndelay=0,
                          correction='cut',
                          Tlim=3e3,
                          Harv=0.006,
                          NbIPM=25,
                          BAsup=190,
                          NClim=3,
                          path='Main',
                          SAVE=FALSE,
                          IsServ=FALSE){
    DFT <- XET <- NULL
    if (IsServ==TRUE){
        OUT <- lapply(1:3, Launch_Simul, spsel=spsel, SurfEch=SurfEch, Ndelay=Ndelay,
                      correction=correction, Tlim=Tlim, Harv=Harv, NbIPM=NbIPM, BAsup=BAsup, path=path)
    }else{
        cl <- makeForkCluster(3)
        OUT <- parLapply(cl, 1:3, Launch_Simul, spsel=spsel, SurfEch=SurfEch, Ndelay=Ndelay,
                         correction=correction, Tlim=Tlim, Harv=Harv, NbIPM=NbIPM, BAsup=BAsup, path=path)
        stopCluster(cl)
    }
    for (iclim in 1:NClim){
        DFT <- rbind(DFT, OUT[[iclim]]$DFT)
        XET <- rbind(XET, OUT[[iclim]]$Xe)
    }
    OUT <- list(DFT=DFT, XET=XET)
    if (SAVE==TRUE){
        if (correction==''){
            Suff <- paste0('_Harv',gsub('[.]','',as.character(Harv)))
        }else if (correction=='cut'){
            Suff <- paste0('_Harv',gsub('[.]','',as.character(Harv)),'_Corrcut')
        }
        if (SurfEch==1){Suff <- paste0(Suff, '_Surf1')}
        if (Ndelay>0){Suff <- paste0(Suff, '_', Ndelay)}
        saveRDS(file=file.path('outputSimu', path, paste0('SimuDem_', spsel, Suff, '.Rds')), OUT)
    }else{
        return(OUT)
    }
}


Launch_Simul <- function(spsel,
                         iClim,
                         SurfEch=1,
                         Ndelay=0,
                         correction='cut',
                         Tlim=3e3,
                         Harv=0.006,
                         NbIPM=25,
                         BAsup=190,
                         path='Main'){
    DFT <- NULL
    XET <- NULL
    mod <- readRDS(file.path('output', path, paste0('fit_sgr_all_', spsel, '.Rds')))
    dataPl <- readRDS(file.path('output', path, paste0('data_plots_pred_', spsel, '.Rds')))
    list_covs <- filter(dataPl, N==iClim)
    XeAllN <- getSim(spsel=spsel, type='DistrEq', Harv=0.006, correction=correction, Type='NL', Ndelay=Ndelay, path=path) %>%
        filter(clim==iClim)
    BAeAllN <- getSim(spsel=spsel, type='BAeq', Harv=0.006, correction=correction, Type='NL', Ndelay=Ndelay, path=path) %>%
        filter(clim==iClim)
    Fipm <- readRDS(file.path('output', path, paste0('FullIPM_', spsel, '_Clim_', iClim, '.Rds')))
    for (iN in 1:NbIPM){
        fit_sgr <- mod[[iN]]
        Xe <- filter(XeAllN, N==iN)
        mesh <- Xe$mesh
        BAe <- filter(BAeAllN, N==iN)
        Ne <- sum(Xe$Xe) / XeAllN$SurfEch[1] * SurfEch
        if (is.na(BAe$BA)|BAe$BA<2|BAe$BA>199|is.na(Ne)|Ne<1){
            DFp <- data.frame(T=1:Tlim, BA=NA, Nt=NA, NHAt=NA, NewRec=NA, V=NA, N=iN, BAe=NA, NHAe=NA)
            Xep <- data.frame(Xe=NA, N=iN, clim=iClim, sp=spsel, Ndelay=Ndelay, SurfEch=SurfEch)
        }else{
            Model <- Fipm[[iN]]
            UVx <- CalcVx(Model, BAe$BA, Harv=Harv, correction='cut', SurfEch=SurfEch, Ndelay=Ndelay)
            Xini <- rmultinom(1, round(Ne), Xe$Xe/sum(Xe$Xe))
            Popini <- XtoPop(Xini, mesh)
            Out <- Stoch_simul(Popini, fit_sgr=fit_sgr, list_covs=list_covs,
                               SurfEch=SurfEch, Ndelay=Ndelay, Tlim=Tlim, Harv=Harv, Vx=UVx$Vx, mesh=mesh, BAsup=BAsup)
            DFp <- data.frame(T=1:Tlim, BA=Out$BA, Nt=Out$Nt, NHAt=Out$NHAt, NewRec=Out$NewRec, V=Out$V,
                              N=iN, BAe=BAe$BA, NHAe=Ne/SurfEch)
            Xep <- data.frame(Xe=Out$Xe) %>% mutate(N=iN, sp=spsel, Ndelay=Ndelay, SurfEch=SurfEch, clim=iClim)
        }
        XET <- rbind(XET, Xep)
        DFT <- rbind(DFT, DFp)
    }
    DFT <- mutate(DFT, clim=iClim, sp=spsel, Ndelay=Ndelay, SurfEch=SurfEch)
    OUT <- list(DFT=DFT, Xe=XET)
    return(OUT)
}


Stoch_simul <- function(Popini,
                        fit_sgr,
                        list_covs,
                        SurfEch=1,
                        Ndelay=0,
                        Tlim=3e3,
                        Harv=0.006,
                        IsDeath=TRUE,
                        Vx=NA,
                        mesh=NA,
                        BAsup=190){
    Pop <- Popini
    if (Ndelay > 0){NrecMem <- rep(0, Ndelay)}
    Theta_rec <- fit_sgr$rec$sigma
    Sig_gr <- fit_sgr$gr$sigma
    BAt <- Nt <- NewRec <- SizeMean <- Vt <- rep(0, Tlim)
    K_rec <- get_model_params(fit_sgr$rec, list_covs)
    K_gr <- get_model_params(fit_sgr$gr, list_covs)
    K_sv <- get_model_params(fit_sgr$sv, list_covs)
    ListCovs <- unique(c(as.character(K_rec$Var),
                         as.character(K_gr$Var), as.character(K_sv$Var)))
    for (t in 1:Tlim){
        ListBA <- pi * (Pop/2000)^2 / SurfEch
        BA <- sum(ListBA)
        size <- Pop
        if (length(Pop) > 5e4 | BA>190){
            print('Population too large')
            return(list(BA=BAt, Nt=Nt, NHAt=Nt/SurfEch, NewRec=NewRec, SizeMean=SizeMean, V=Vt))
        }
        if (length(Pop)<=0){
            print("Population went extinct")
            return(list(BA=BAt, Nt=Nt, NHAt=Nt/SurfEch, NewRec=NewRec, SizeMean=SizeMean, V=Vt))
        }

        Covs <- data.frame(size=Pop, logsize=log(Pop), BATOTcomp=BA,
                           BATOTSP= BA, logBATOTSP=log(BA), BATOTNonSP=0, Intercept=1)
        Covs <- as.data.table(Covs)
        if (t==1){if(prod(names(Covs) %in% ListCovs)==0){stop('Covs data.frame miss a covariate')}}
        BAt[t] <- BA
        Nt[t] <- length(Pop)
        if (!is.na(Vx[1])){Vt[t] <- PoptoV(Pop, mesh, Vx)}
        SizeMean[t] <- mean(Pop)
        ## Computation of expected growth, survival and recruitment
        GrMean <- rep(0, length(Pop))
        for (i in K_gr$Var){
            GrMean <- GrMean + Covs[, ..i][[1]] * K_gr[Var==i, 2][[1]]
        }
        SvMean <- rep(0, length(Pop))
        for (i in K_sv$Var){
            SvMean <- SvMean + Covs[, ..i][[1]] * K_sv[Var==i, 2][[1]]
        }
        RecMean <- 0
        for (i in K_rec$Var){
            RecMean <- RecMean + Covs[1, ..i][[1]] * K_rec[Var==i, 2][[1]]
        }
        P_sv <- (1 - fit_sgr$sv$family$linkinv(SvMean)) * (1-Harv)
        ## Sampling from Growth, Survival and Recruitment distribution
        Nrec <- rnbinom(1, mu=exp(RecMean) * SurfEch/0.03, size=Theta_rec)
        if(is.na(Nrec)){
            stop(paste0('Error Nrec NA at time', t, 'with mu=', exp(RecMean), 'Pop=', Pop, 'BA=', BA))
        }
        if (IsDeath==TRUE){
            Surv <- rbinom(length(P_sv), 1, P_sv)
        }else{
            Surv <- rep(1, length(P_sv))
        }
        NewRec[t] <- Nrec
        sizeT1 <- size + rlnorm(length(GrMean), meanlog=GrMean, sdlog=Sig_gr)
        Surv[sizeT1>(fit_sgr$maxDBH * 1.1)] <- 0
        PopT1 <- sizeT1[Surv==1]
        if (Ndelay >0){
            if (NrecMem[(t %% Ndelay + 1)] > 0){PopT1 <- c(rep(100, NrecMem[(t %% Ndelay + 1)]), PopT1)}
        }else{
            PopT1 <- c(rep(100, Nrec), PopT1)
        }
        Pop <- PopT1
        #:::::::::::::::::
        if (Ndelay > 0){
            NrecMem[(t %% Ndelay + 1)] <- Nrec
        }
        #:::::::::::::::::
    }
    return(list(BA=BAt, Nt=Nt, NHAt=Nt/SurfEch, NewRec=NewRec, SizeMean=SizeMean, Xe=Pop, V=Vt))
}

get_model_params <- function(s_res,
                             list_covs){
    require(data.table)
    params_sv <- s_res$params_m
    Iintercept <- which(grepl('ntercep', names(params_sv))==TRUE)
    K_intercept <- params_sv[Iintercept]
    params_sv <- params_sv[-Iintercept]

    DF_covs <- data.frame(name=names(list_covs), val=as.numeric(list_covs))
    tt <- function(N,split,k){strsplit(N, split=split)[[1]][k]}
    TabParams <- data.frame(V1=as.character(do.call(rbind,lapply(names(params_sv), tt, split=':', k=1))),
                            V2=as.character(do.call(rbind,lapply(names(params_sv), tt, split=':', k=2))),
                            K=as.numeric(params_sv)) %>% left_join(DF_covs, by=c('V1'='name')) %>%
        left_join(DF_covs, by=c('V2'='name'))
    TabParams <- as.data.table(TabParams)
    ListVarModel <- unique(c(as.character(TabParams$V1), as.character(TabParams$V2)))
    ListVar <- ListVarModel[!(ListVarModel %in% c(names(list_covs), NA))]
    K_Var <- data.frame(Var=ListVar, K=NA)
    for (i in 1:length(ListVar)){
        K_Var$K[i] <- sum(apply(TabParams[V1==ListVar[i]|V2==ListVar[i],3:5], 1, prod, na.rm=TRUE))
    }
    K_i <- apply(TabParams[!(V1 %in% c(ListVar)) & !(V2 %in% c(ListVar)),3:5], 1, prod, na.rm=TRUE)
    K_intercept <- K_intercept + sum(K_i)
    K <- rbind(K_Var, data.frame(Var='Intercept', K=K_intercept))
    return(as.data.table(K))
}

XtoPop <- function(X, mesh){ # transform X vector (stochastic) to Population vector
    Pop <- NULL
    for (k in 1:length(X)){
        if (X[k]>0){
            Pop <- c(Pop, rep(mesh[k], X[k]))
        }
    }
    return(Pop[Pop>0])
}

PoptoX <- function(Pop, mesh){ # transform Population vector to X vector
    X <- 0 * mesh
    Pop <- Pop[Pop>0]
    Int1 <- apply(abs(outer(Pop, mesh, '-')), 1, which.min)
    for (k in 1:length(Pop)){
        X[Int1[k]] <- X[Int1[k]] + 1
    }
    return(X)
}

PoptoV <- function(Pop, mesh, Vx){
    X <- PoptoX(Pop, mesh)
    V <- as.numeric(Vx %*% X)
    return(V)
}

# Function to retrieve specific outputs from Non Stochastic simulations (BAeq, Jac, DistrEq, BAALL, NtALL)
getSim <- function(spsel, type='BAeq', Harv=0.006, correction='', Type='NL', Ndelay=0, path=''){
    Sim <- readSimu(spsel=spsel, correction=correction, type='NonDem', Ndelay=Ndelay, path=path)
    Tlim <- max(Sim$BA$T)
    if (type=='BAeq'){
        Out <- filter(Sim$BAe, type==Type) %>% mutate(sp=spsel,Ndelay=Ndelay)
    }else if (type=='Jac'){
        Out <- filter(Sim$JacEq, type==Type) %>% mutate(sp=spsel,Ndelay=Ndelay)
    }else if (type=='DistrEq'){
        Df <- filter(Sim$Xe, type==Type)
        Nsim <- length(unique(Df$Sim))
        Nmesh <- dim(filter(Df, N==1, clim==1))[1]/Nsim
        Out <- group_by(Df, N, clim) %>%
            mutate(meshi=rep(1:Nmesh,Nsim)) %>% ungroup()
        Out <- group_by(Out, N, clim, meshi) %>%
            summarise(Xe=Xe[1], mesh=mesh[1]) %>% ungroup() %>%
            mutate(sp=spsel,Ndelay=Ndelay, SurfEch=Sim$SurfEch)
    }else if (type=='BAALL'){
        Out <- filter(Sim$BA,type==Type) %>% mutate(sp=spsel,Ndelay=Ndelay)
    }else if (type=='NtALL'){
        Out <- filter(Sim$Nt,type==Type) %>% mutate(sp=spsel,Ndelay=Ndelay)
    }else{
        Out <- NA
    }
    return(Out)
}

# Function to retrieve outputs from simulations (Dem, NonDem, Sigma and Perturb)
readSimu <- function(spsel,
                     Harv=0.006,
                     correction='',
                     SurfEch=300*1e-4,
                     path="",
                     Ndelay=0,
                     type='NonDem'){
    SuffCut <- ''
    SuffSurf <- ''
    SuffDelay <- ''
    SuffHarv <- gsub('[.]', '', as.character(Harv))
    if (correction=='cut'){SuffCut <- '_Corrcut'}
    if (type %in% c('Extinc', 'Dem')){
        if (SurfEch==0.03){
            SuffSurf <- '_Surf003'
        }else if (SurfEch==0.01){
            SuffSurf <- '_Surf001'
        }else if (SurfEch==0.1){
            SuffSurf <- '_Surf01'
        }else if (SurfEch==1){
            SuffSurf <- '_Surf1'
        }
    }
    if (Ndelay>0){SuffDelay <- paste0('_', Ndelay)}
    if (type=='Extinc'){
        DF <- NULL
        for (iC in (1:3)){
            filename <- file.path('outputSimu', path, paste0("Extinc_", spsel,'_',
                                                             iC, '_Harv', SuffHarv, SuffCut, SuffSurf, SuffDelay, '.Rds'))
            if (!(file.exists(filename))){
                print(paste0('Missing file ', filename))
                DF <- rbind(DF, data.frame(N=1:100, Nt=NA, Te=NA, BAe=NA, Npope=NA, NINI=NA, BAINI=NA, IsIni=FALSE, PrEx=NA, clim=iC, sp=spsel))
            }else{
                DF <- rbind(DF, mutate(readRDS(filename), clim=iC))
            }
        }
    }else if (type=='Sigma'){
        filename <- file.path('outputSimu', path, paste0('Sigma_', spsel, SuffCut, SuffDelay, '.Rds'))
        if (!(file.exists(filename))){
            print(paste0('Missing file ', filename))
            DF <- NULL
        }else{
            DF <- readRDS(filename)
        }
    }else if (type=='Perturb'){
        filename <- file.path('outputSimu', path , paste0('Perturb_', spsel, '.Rds'))
        if (!(file.exists(filename))){
            print(paste0('Missing file ', filename))
            return(NULL)
        }
        DF <- readRDS(filename)
    }else{
        filename <- file.path('outputSimu', path,
                              paste0("Simu", type, '_', spsel, '_Harv', SuffHarv, SuffCut, SuffSurf, SuffDelay, '.Rds'))
        if (!(file.exists(filename))){
            stop(paste0('No file found : ', filename))
        }
        DF <- readRDS(filename)
        if (type=='Dem'){
            names(DF)[1:2] <- c('DFT', 'XET')
            names(DF$DFT)[which(names(DF$DFT)=='Npop')] <- 'Nt'
        }
    }
    return(DF)
}

## Function to compute Left and reight eigenvectors of equilibrium state (normalized)
CalcVx <- function(Model, # require the mesh, LIPM and RecFun only
                   BAeq,
                   Harv=0.006,
                   correction='',
                   SurfEch=300*1e-4,
                   Ndelay=0){

    if (is.na(BAeq) | BAeq<2 | BAeq>199){return(list(Vx=rep(NA, length(Model$meshpts)+Ndelay), Ux=rep(NA, length(Model$meshpts)+Ndelay)))}
    mesh <- Model$meshpts
    ct <- Buildct(mesh,SurfEch)
    if (correction=='cut'){ct[length(mesh)] <- 0}
    LIPM <- Model$LIPM
    if (max(LIPM[[20]]) > 1e6){for (i in 1:length(LIPM)){LIPM[[i]] <- LIPM[[i]] * 1e-7}}
    if (correction=='cut'){for (i in 1:length(LIPM)){LIPM[[i]][,length(mesh)] <- 0;LIPM[[i]][length(mesh),] <- 0}}
    b <- c(rep(1/2,2),rep(0,length(mesh)-2))
    NIPM <- as.integer(floor(BAeq))
    NIPM <- min(NIPM,200)
    NIPM <- max(NIPM,1)
    RecMean <- exp(Model$RecFun(BAeq)) * SurfEch / (300 *1e-4)
    IPM0 <- as.matrix(LIPM[[NIPM]] * (1 - (BAeq-NIPM)) + LIPM[[NIPM+1]] * (BAeq-NIPM)) * (1-Harv)
    if (Ndelay>0){
        ct <- t(c(rep(0, Ndelay),ct))
        b <- c(b,rep(0, Ndelay))
        IPM0 <- cbind(matrix(0, nrow=length(mesh)+Ndelay, ncol=Ndelay),rbind(matrix(0, ncol=length(mesh), nrow=Ndelay), IPM0))
        for (i in 1:Ndelay){IPM0[i+1, i] <- 1}
        # TODO : equivalent a Add_Delay_Matrix
    }
    IPM <- IPM0 + RecMean * b %*% ct / BAeq
    Le <- Re(eigen(t(IPM))$vectors[,1])
    Re <- Re(eigen(IPM)$vectors[,1])
    Re <- Re / sum(Re)
    Le <- Le / as.numeric(Le %*% Re)
    return(list(Vx=Le,Ux=Re))
}
