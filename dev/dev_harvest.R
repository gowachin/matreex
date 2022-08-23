
#' Compute the basal area to be cut to reach BAtarget
# x is the size distribution
# ct is the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
# BAtarget is the basal area to reach after maximum cut
# Pmax is the maximum proportion of BAcut / BA
# dBAmin is the minimum BA to perform cut
getBAcutTarget <- function(x,
                           ct,
                           BAtarget=20,
                           Pmax=0.25,
                           dBAmin=3){
    BA <- as.numeric(ct %*% x)
    if (BA==0){return(BAcut=0)}
    if ((BA - BAtarget) < dBAmin){return(BAcut=0)}
    Htarget <- (BA - BAtarget) / BA
    H <- min(Pmax, Htarget)
    return(BAcut=BA*H)
}

#' Compute the harvest curve for uneven stand (Pcut)
# x is the size distribution
# h is the simple harvest function (power=1)
# ct is the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
# hmax is the maximum harvest rate for a size class
# BAtarget is the basal area to reach after maximum cut
# Pmax is the maximum proportion of BAcut / BA
# dBAmin is the minimum BA to perform cut
## Two cases: BA above dha is enough to reach target and we only cut above
##              dbh such as we get target
##            BA above dha is not enough, we optimize power k of h
##              function to get target
getPcutUneven <- function(x,
                          h,
                          ct,
                          hmax=1,
                          BAtarget=20,
                          Pmax=0.25,
                          dBAmin=3,
                          Print=FALSE){
    BAcuttarget <- getBAcutTarget(x=x, ct=ct, BAtarget=BAtarget, Pmax=Pmax,
                                  dBAmin=dBAmin)
    getBAha <- function(x, h, ct, hmax=1){as.numeric(ct %*% (hmax*(h==1)*x))}
    getBAcut <- function(k, x, h, ct, hmax=1){as.numeric(ct %*% (h^k * x * hmax0))}
    getCostfunc <- function(k, x, h, ct, hmax=1, BAcuttarget){
        return((getBAcut(k, x, h, ct, hmax) - BAcuttarget)^2)
    }
    BAha <- as.numeric(getBAha(x=x, h=h, ct=ct, hmax=hmax))
    if (BAha >= BAcuttarget){
        Pcut <- (cumsum(ct[length(ct):1] * x[length(ct):1] * hmax) -
                     rep(BAcuttarget, length(ct)))<=0
        Pcut <- hmax * Pcut[length(ct):1]
    }else{
        kopt <-  optimize(f=getCostfunc, lower=0.1, upper=10, x=x, h=h,
                          ct=ct, hmax=hmax, BAcuttarget=BAcuttarget)$minimum
        Pcut <- h^kopt*hmax
    }
    if (Print==TRUE){
        print(paste0('BA initial : ', round(ct %*% x, digits=1),
                     ' / H target : ', round(BAcuttarget/(ct%*%x), digits=2)))
        print(paste0('BA cut : ', round(ct %*% (x*Pcut), digits=1),
                     ' / Target : ', round(BAcuttarget, digits=1)))
    }
    return(Pcut)
}


#' Perform the harvest Uneven
# x is the size distribution
# mesh is the meshpts associated with x
# ct is the vector to compute BA with x (ct=Buildct(mesh, SurfEch))
# dth is the minimum diameter at which we cut (same unit as x)

# dha is the harvest diameter
# hmax is the maximum harvest rate for a size class
# BAtarget is the basal area to reach after maximum cut
# Pmax is the maximum proportion of BAcut / BA
# dBAmin is the minimum BA to perform cut
HarvestIPMUneven <- function(x,
                             mesh,
                             ct,
                             dth=175,
                             dha=575,
                             hmax=1,
                             BAtarget=20,
                             Pmax=0.25,
                             dBAmin=3,
                             Print=FALSE){
    h <- (mesh-dth) / (dha - dth)
    h <- pmax(0, h)
    h <- pmin(1, h)
    # we only consider tree above dth in the algorithm
    x1 <- x
    x1[mesh<=dth] <- 0
    Pcut <- getPcutUneven(x=x1, h=h, ct=ct, hmax=hmax, BAtarget=BAtarget,
                          Pmax=Pmax, dBAmin=dBAmin, Print=Print)
    x <- x * (1 - Pcut)
    return(x)
}

