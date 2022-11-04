require(dplyr)
require(stringr)
require(tidyr)

# Idea: instead of performing the growth integration for each IPM conditions,
# We first create for a given species (defining the mesh and the growth variance)
# a matrix containing the integrated growth vectors for different values of growth (mu)
# mu is computed every stepMu, wich has to be defined (here default is 1e-3)
# With the matrix Tabgrowth containing growth integration vector for a large set of
# mu values, we can reconstruct the IPM using the vector of mu for each mesh point
#(at any BA, any climatic conditions, ...) without any additional computer cost
# One potential drawback : we do not integrate anymore on the size at the first time-step,
# can be fixed but need more time and be sure it's valuable
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


fit <- function(d_x1_x, mu, sig){
    out <- rep(0, length.out = length(d_x1_x))
    out[d_x1_x>0] <- dnorm(log(d_x1_x[d_x1_x>0]), mu[d_x1_x>0], sig) * 1/d_x1_x[d_x1_x>0]
    return(out)
}

# Integration with a part in GL, the other in midbin
getInt <- function(mu,
                   sig,
                   mesh,
                   level=120,
                   IntQuad=30,
                   IntTotal=300){
    meshQuad <- mesh[1:IntQuad]
    outQuad <- gaussQuadInt(-h/2, h/2, floor(level/3))
    meshQuad2 <- outer(meshQuad, outQuad$nodes, '+')
    tt <- structure(vapply(meshQuad2, FUN=fit, mu=mu, sig=sig, numeric(1)), dim=dim(meshQuad2))
    ValQuad <- tt %*% outQuad$weights
    meshMid <- mesh[(IntQuad+1):IntTotal]
    ValMid <- sapply(meshMid, fit, mu=mu, sig=sig) * h
    DF <- data.frame(mesh=mesh, Int=c(ValQuad, ValMid, rep(0, length(mesh)-length(ValQuad)-length(ValMid))),
                     method=c(rep('Quad', length(ValQuad)), rep('Mid', length(ValMid)), rep('None',
                                                                                            length(mesh)-length(ValQuad)-length(ValMid))))
    return(DF[, 2])
}

# Compute the growth integration for a vector of growth mu
# return a data.frame, each line is the growth integration for a given mu
# NB:we do not integrate the diagonal term: it is calculated as 1 - sum(Pi)
# Lot faster (level can be extremely lower), but dangerous as we cannot check properly
# integration, to be discussed
getTabGrowth <- function(minMu=-4,
                         maxMu=3,
                         stepMu=1e-3,
                         spsel='Picea abies',
                         level=120,
                         IntTotal=300){
    fi <- readRDS(paste0('~/Documents/IPM/full_ipm/output/New100/fit_sgr_all_', spsel, '.Rds'))
    Sig <- fi[[1]]$rec$sigma
    m_size <- 700
    L <- 90
    U <- 1.1 * fi[[1]]$maxDBH
    h <- (U - L) / m_size
    mesh <- seq(0, (IntTotal)*h, by=h)
    tabMu <- seq(minMu, maxMu, by=stepMu)
    TabGrowth <- t(do.call(rbind, lapply(tabMu, getInt, sig=Sig, mesh=mesh[-1],
                                         level=level, IntQuad=30, IntTotal=IntTotal)))
    TabGrowth <- rbind((1-apply(TT, 2, sum)), TT)
    return(t(TabGrowth))
}

# return the line in TabGrowth corresponding to the argument mu
LinkMuLine <- function(mu,
                       minMu=-4,
                       stepMu=1e-3){
    return(pmax(1, 1+(mu-minMu)/stepMu))
}

# Construct IPM using mumesh : vector of size mesh
# with the mu value calculated for each grid poinr
getIPM <- function(mumesh,
                   TabGrowth){
    Lines <- LinkMuLine(mumesh)
    IPM <- matrix(0, ncol=length(mumesh), nrow=length(mumesh))
    for (k in 1:length(mumesh)){
        rowMax <- min(length(mumesh), k+dim(IPMGr)[2]-1)
        IPM[k:rowMax, k] <- TabGrowth[Lines[k], 1:(rowMax-k+1)]
    }
    return(IPM)
}

#################################
####### Get the range of mu values
#################################
getMaxMinmu <- function(spsel='Picea abies',
                        Nmodel=1){
    print(spsel)
    meshSize <- seq(90, 900, by=10)
    fi <- readRDS(paste0('~/Documents/IPM/full_ipm/output/New100/fit_sgr_all_', spsel, '.Rds'))
    da <- readRDS(paste0('~/Documents/IPM/full_ipm/output/New100/data_plots_pred_', spsel, '.Rds'))
    for (Nc in 1:3){
        da1 <- pivot_longer(dplyr::select(filter(da, N==Nc), -N, -SDM, -PC1, -PC2), sgdd:sgdd2)
        sig <- fi[[Nmodel]]$gr$sigma
        params <- fi[[Nmodel]]$gr$params_m
        params <- data.frame(name=names(params), Val=as.numeric(params))
        Interact <- sum(grepl(':', params$name)) == 0
        minmu <- 100
        maxmu <- -100
        for (iBA in seq(0, 200, by=10)){
            da1i <- rbind(da1, data.frame(name='BATOTcomp', value=iBA))
            if (Interact==FALSE){
                par <- mutate(params, name1=str_split(name, ':', simplify=TRUE)[,1],
                              name2=str_split(name, ':', simplify=TRUE)[, 2]) %>%
                    left_join(da1i, by=c('name1'='name')) %>% left_join(da1i, by=c('name2'='name'))
            }else{
                par <- mutate(params, name1=str_split(name, ':', simplify=TRUE)[,1]) %>%
                    left_join(da1i, by=c('name1'='name'))
                names(par)[4] <- 'value.x'
                par <- mutate(par, value.y=NA)
            }
            par$value.x[is.na(par$value.x)] <- 1
            par$value.y[is.na(par$value.y)] <- 1
            Size <- filter(par, grepl('size', name)) %>% filter(!grepl('logsize', name))
            IntSize <- sum(Size$Val*Size$value.x*Size$value.y)
            LogSize <- filter(par, grepl('logsize', name))
            IntLogSize <- sum(LogSize$Val*LogSize$value.x*LogSize$value.y)
            Int <- filter(par, !grepl('size', name))
            IntInt <- sum(Int$Val*Int$value.x*Int$value.y)
            mu <- IntInt + IntSize*meshSize + IntLogSize*log(meshSize)
            minmu <- min(c(minmu, mu))
            maxmu <- max(c(maxmu, mu))
        }
    }
    return(data.frame(sp=spsel, minmu=minmu, maxmu=maxmu, sig=sig))
}

F <- function(spsel){
    cl <- makeForkCluster(5)
    H <- do.call(rbind, parLapply(cl, seq(1, 100, length.out=5), getMaxMinmu, spsel=spsel))
    stopCluster(cl)
    return(H)
}

G <- function(){
    Listsp <- str_split(list.files(path='~/Documents/IPM/full_ipm/output/New100', pattern='1.Rds'),
                        '_', simplify=TRUE)[, 2]
    H <- do.call(rbind, lapply(Listsp, F))
    return(H)
}

# Range of mu : -10 / 3







