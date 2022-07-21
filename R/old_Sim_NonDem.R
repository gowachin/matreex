JacobianIN <- function(Model,
                       BAeq,
                       Xeq,
                       Harv=0.006,
                       correction='',
                       SurfEch=300*1e-4,
                       Ndelay=10,
                       OUT='Metrics'){
    # TODO : library(Matrix) for %*%
    mesh <- Model$meshpts
    Nlength <- length(mesh)
    RecMean <- Model$RecFun
    ct <- Buildct(mesh,SurfEch)
    if (correction=='cut'){ct[Nlength] <- 0}
    LIPM <- Model$LIPM
    # if (max(LIPM[[20]]) > 1e6){for (i in 1:length(LIPM)){LIPM[[i]] <- LIPM[[i]] * 1e-7}} # TODO :why not do this prior to sim and jacobian ??
    if (correction=='cut'){for (i in 1:length(LIPM)){LIPM[[i]][,Nlength] <- 0;LIPM[[i]][Nlength,] <- 0}}
    b <- c(rep(1/2,2),rep(0,Nlength-2))
    NIPM <- as.integer(floor(BAeq))
    if (NIPM>190|NIPM<2){ # n'arrive pas car BAeq > 2 dans Run_Sim_NonDem_NL
        Eig <- data.frame(Period=NA, Damp=NA, Reactivity=NA, Reactivity2=NA, DampRatio=NA, Vec=NA)

        # # Maxime add to remove else part
        # if(OUT != "Metrics"){
        #     stop(sprintf("This is not possible for such BA value : %i", NIPM))
        # }
        # return(Eig)

    }else{
        xe <- Xeq
        Dr <- (exp(RecMean(BAeq+1))*SurfEch/(300*1e-4)-exp(RecMean(BAeq-1))*SurfEch/(300*1e-4))/2
        DIPM  <- (Matrix(LIPM[[NIPM+1]],sparse=FALSE) - Matrix(LIPM[[NIPM-1]],sparse=FALSE))/2
        DIPM <- DIPM * (1-Harv)
        IPM0 <- LIPM[[NIPM]] * (1-(BAeq-NIPM)) + LIPM[[NIPM+1]] * (BAeq-NIPM)
        IPM <- as.matrix(IPM0) * (1-Harv) + b %*% ct * exp(RecMean(BAeq)) * SurfEch / (300*1e-4) / BAeq
        if (Ndelay>0){
            DIPM <- rbind(matrix(0, ncol=Nlength+Ndelay,nrow=Ndelay), cbind(matrix(0, ncol=Ndelay, nrow=Nlength), DIPM))
            IPM <- rbind(matrix(0, ncol=Nlength+Ndelay, nrow=Ndelay), cbind(matrix(0, ncol=Ndelay, nrow=Nlength), IPM))
            for (i in 1:Ndelay){IPM[(i+1),i] <- 1}
            # TODO : this is Add_Delay_Matrix function
            # DIPM <- Add_Delay_Matrix(IPM, Nlength, Ndelay)
            # TODO : this is Add_Delay_Matrix function without the diagonal (add bool ?)
            ct <- t(c(rep(0, Ndelay), ct))
            mesh <- c(rep(0, Ndelay), mesh)
            b <- c(b, rep(0, Ndelay))
        }
        JFull <- Matrix(IPM,sparse=FALSE) - diag(Nlength+Ndelay) + Dr * b %*% ct + DIPM %*% xe %*% ct
        if (OUT=='Matrix'){return(as.matrix(JFull))}
        if (correction=='cut'){
            Nlength <- Nlength - 1
            EigFull <- eigen(JFull[1:(Nlength+Ndelay), 1:(Nlength+Ndelay)])
        }else{
            EigFull <- eigen(JFull)
        }
        imax <- which.max(Re(EigFull$values))
        Vec <- unique(abs(eigen(IPM)$values))
        Vec <- Vec[order(Vec,decreasing=TRUE)]
        DampRatio <- Vec[1] / Vec[2]
        B <- JFull + diag(dim(JFull)[1])
        Reactivity <- max((eigen(t(B)%*%B)$values))
        Reactivity2 <- max(abs(eigen(B)$values))
        Reactivity3 <- max(colSums(B))
        Eig <- data.frame(Period=abs(2*pi/Im(EigFull$values)),
                          Damp=log(0.5) / Re(EigFull$values), Reactivity=Reactivity,
                          Reactivity2=Reactivity2, DampRatio=DampRatio, Vec=EigFull$vectors[,imax])
    }

    if (OUT=='Metrics'){
        return(Eig)
    }else{
        # TODO : impossible exit if NIPM >190 | NIPM <2
        return(EigFull)
    }
}

