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

    # TEMP dev
    mu <- 0
    sig <- 0.6
    h <- 2
    mesh <- seq(0, 16, by = h)
    level=120
    IntQuad=3
    IntTotal= 8
    Level <- 5
    N_int <- 5
    # TEMP dev

    meshQuad <- mesh[1:IntQuad]
    outQuad <- gaussQuadInt(-h/2, h/2, floor(level/3))
    meshQuad2 <- outer(meshQuad, outQuad$nodes, '+')
    tt <- structure(vapply(meshQuad2, FUN=fit, mu=mu, sig=sig, numeric(1)), dim=dim(meshQuad2))
    ValQuad <- tt %*% outQuad$weights


    meshMid <- mesh[(IntQuad+1):IntTotal]

    meshMid <- seq(mesh[(IntQuad+1)]- h/2, mesh[IntTotal] + h/2, by=h/Level)[-1]

    ca <- factor(rep(1:N_int, each = Level))
    ca <- .Internal(split(1:length(meshMid), ca))
    ValMid <- sapply(meshMid, fit, mu=mu, sig=sig) * h / Level
    ValMid <- unlist(lapply(ca, function(i) sum(ValMid[i])))
    # DF <- data.frame(mesh=mesh, Int=c(ValQuad, ValMid, rep(0, length(mesh)-length(ValQuad)-length(ValMid))),
    #                  method=rep(c("Quad", "Mid", "None"),
    #                             times = c(length(ValQuad), length(ValMid),
    #                                       length(mesh)-length(ValQuad)-length(ValMid)) ))

    res <- c(ValQuad, ValMid, rep(0, length(mesh)-length(ValQuad)-length(ValMid)))
    # return(DF[, 2])

    return(res)
}



# Compute the growth integration for a vector of growth mu
# return a data.frame, each line is the growth integration for a given mu
# NB:we do not integrate the diagonal term: it is calculated as 1 - sum(Pi)
# Lot faster (level can be extremely lower), but dangerous as we cannot check properly
# integration, to be discussed
getTabGrowth <- function(minMu=-4,
                         maxMu=3,
                         stepMu=1e-3,
                         species='Picea_abies',
                         level=120,
                         IntTotal=300){
    # TEMP dev
    minMu=-4
    maxMu=3
    stepMu=1e-3
    species='Picea_abies'
    level=120
    IntTotal=300
    # TEMP dev

    data(list = paste0("fit_", species))
    fit <- eval(parse(text=paste0("fit_", species)))
    sig_gr <- fit$gr$sig_grma
    m <- 700
    L <- 90
    U <- as.numeric(fit$info[["max_dbh"]]) * 1.1
    h <- (U - L) / m
    mesh <- seq(0, (IntTotal)*h, by=h)
    mu_tab <- seq(minMu, maxMu, by=stepMu)

    tic()
    tmp <- lapply(mu_tab, getInt, sig = sig_gr, mesh = mesh[-1],
                  level = level, IntQuad = 30, IntTotal = IntTotal)
    toc() # 32 sec
    tic()
    test <- map(mu_tab, ~ getInt(.x, sig = sig_gr, mesh = mesh[-1],
                                level = level, IntQuad = 30, IntTotal = IntTotal))
    toc()

    TabGrowth <- t(do.call(rbind, tmp))
    # TabGrowth <- rbind((1-apply(TT, 2, sum)), TT) # FIXME what is TT ??
    return(t(TabGrowth))
}


x <- getTabGrowth()

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
getRangemu <- function(species='Picea_abies',
                       fit, BA = 0:200, mesh_x = seq(90, 900, by = 10)){

    assertCharacter(species, len = 1)
    assertCount(BA, lower = 0, upper = 200)
    assertNumeric(mesh, lower = 0)

    climate_species <- subset(treeforce::climate_species,
                              sp == species, select = -sp)

    fres <- data.frame(min = 1:3, max = 1:3)
    for (Nc in 1:3){
        climate <- subset(climate_species, N == Nc, select = -N)
        climate <- drop(as.matrix(climate)) # we need it as a vector.
        list_covs <- c(climate, BATOTcomp = 0)
        # browser()
        res <- matrix(ncol = 2, nrow = length(BA),
                      dimnames = list(NULL, c("min", "max")))

        for (iBA in seq_along(BA)){

            list_covs["BATOTcomp"] <- BA[iBA]

            grFun <- exp_sizeFun(fit$gr$params_m, list_covs)
            mu <- grFun(mesh_x)

            res[iBA,] <- range(mu)
        }
        fres[Nc, ] <- c(min(res[,"min"]), max(res[,"max"]))

    }
    range <- c(min = min(fres$min), max= max(fres$max), sig = fit$gr$sigma)

    return(range)
}

load_all()
species <- "Picea_abies"
data(list = paste0("fit_", species))
fit <- eval(parse(text=paste0("fit_", species)))
tic()
(mu_range <- getRangemu(fit = fit, BA = seq(0, 200, by = 10), mesh = seq(90, 1900, by = 2)))
toc()

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



# Testing personnal function ####

load_all()
x <- make_mutrix(species = "Abies_alba", fit_Abies_alba, verbose = TRUE)
# Mu range done
# Launching mu computation loop
# GL integration occur on 25 cells
# midbin integration occur on 25 cells
# Loop done.
# Time difference of 2.45 mins

lobstr::obj_size(x) # 2.87 MB

# take as long as integration but is less dense.
x$mutrix %>%
    reshape2::melt() %>%
    ggplot(aes(x = Var1, y = Var2, fill = log(value))) +
    geom_tile() +
    scale_fill_viridis_c(na.value="transparent") +
    theme_dark() +
    NULL

x$mutrix[4000, 40] # why is it 0 ? So smoll values ?

