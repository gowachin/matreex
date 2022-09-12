rm(list = ls())
library(profvis)
library(tictoc)
source("dev/dev_old_mkipm.R")
#
# profvis({
#     res <- make_FullIPM_iClim_resample("Yggdrasil",
#                                        iClim = 1,
#                                         NbIPM = 1,
#                                        NbModel = 1,
#                                        m_size = 700)
# })


# in make_FullIPM_iClim
spsel <- "Yggdrasil"
iClim <- 1
NbIPM <- 2
NbModel <- 1
m_size <- 700
data_plots_pred <- readRDS(file.path('output', spsel, 'plots_pred.Rds'))
fit_sgr <- readRDS(file.path('output', spsel, 'fit_sgr_all.Rds'))[[1]]


level = 420
correction = "constant"
diag_tresh = 50
IsSurv=TRUE

data_plots_pred <- dplyr::filter(data_plots_pred, N==iClim)
# microbenchmark::microbenchmark( # plus rapide
#     filter = filter(data_plots_pred, N==iClim),
#     base = data_plots_pred[data_plots_pred$N == iClim,]
# )
data_plots_predBA <- data.frame(cbind(matrix(rep(as.matrix(data_plots_pred), each=NbIPM),
                                             nrow=NbIPM), 1:NbIPM))
names(data_plots_predBA) <- c(names(data_plots_pred), 'BATOTSP')
data_plots_predBA <- dplyr::mutate(data_plots_predBA, BATOTcomp=BATOTSP)
# # ligne alternative : (dt/8 in microsecond)
# data_plots_predBA <-cbind(data_plots_pred, data.frame(BATOTSP = 1:NbIPM,
#                                                       BATOTcomp = 1:NbIPM))

L <- 100*0.9
U <- 1.1* fit_sgr$maxDBH
h <- (U - L) / m_size
mesh_x <- seq(L+h/2,U-h/2,length.out=m_size)
N_int <- sum((mesh_x-min(mesh_x))<diag_tresh)
out2 <- gaussQuadInt(-h/2, h/2, floor(level/3))
WMat <- build_weight_matrix(out2$weights,N_int)

# in make_IPM_GL_2_i
data_plots_pred <- data_plots_predBA
i = 1

minSize <- 100*0.9
maxSize <- 1.1* fit_sgr$maxDBH
h <- (maxSize - minSize) / m_size
meshpts <- seq(minSize+h/2, maxSize-h/2, length.out = m_size)
list_m <- dplyr::filter(data_plots_pred, BATOTSP==i)
print("before mk_P_GL_2")
midPoint=TRUE
library(Matrix)
library(dplyr)

m = m_size
L = minSize
U = maxSize
g_res = fit_sgr$gr
s_res = fit_sgr$sv
list_covs = list_m

rm(data_plots_pred, data_plots_predBA,fit_sgr, list_m,
   out2, i, iClim, m_size, maxSize, minSize,
   mesh_x, meshpts, N_int, NbIPM, NbModel, spsel)
rm(make_FullIPM_iClim, make_FullIPM_iClim_resample,
   make_IPM_GL_2_i)

load_all()
source("dev/dev_old_mkipm.R")
tic()
set.seed(42)
old <- mk_P_GL_2(m, L, U, g_res, s_res, list_covs, diag_tresh= 50,
                 level=420, correction="none", WMat,
                 IsSurv=TRUE, midPoint=TRUE)
toc()
tic()
set.seed(42)
new <- mk_P_GL_2_stripe(m, L, U, g_res, s_res, list_covs, diag_tresh= 50,
                 level=420, correction="none", WMat,
                 IsSurv=TRUE, midPoint=TRUE)
toc()
all.equal(old, new)


microbenchmark::microbenchmark(
    old <- mk_P_GL_2(m, L, U, g_res, s_res, list_covs, diag_tresh= 50,
                     level=420, correction="none", WMat,
                     IsSurv=TRUE, midPoint=TRUE),
    new <- mk_P_GL_2_stripe(m, L, U, g_res, s_res, list_covs, diag_tresh= 50,
                            level=420, correction="none", WMat,
                            IsSurv=TRUE, midPoint=TRUE),
    times = 10
)
# beep(5)

# Unit: milliseconds
# expr
# old <- mk_P_GL_2(m, L, U, g_res, s_res, list_covs, diag_tresh = 50,      level = 420, correction = "none", WMat, IsSurv = TRUE, midPoint = TRUE)
# new <- mk_P_GL_2_stripe(m, L, U, g_res, s_res, list_covs, diag_tresh = 50,      level = 420, correction = "none", WMat, IsSurv = TRUE, midPoint = TRUE)
# min        lq     mean    median        uq       max neval
# 3037.1792 3084.0582 3242.552 3212.0862 3284.2053 4500.4443   100
# 621.1358  630.6961  694.538  642.0454  767.7873  940.2343   100

microbenchmark::microbenchmark(
    old = fuu(g, ca, h, Level),
    new = foo(g, ca, h, Level),
    best = fii(g, ca, h, Level),
    times = 100
)


# in fun_mid_point
rm(mk_P_GL_2)
h <- (U - L) / m
mesh_x <- seq(L+h/2, U-h/2, length.out=m)
N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
N_ini <- N_int+1
N_int <- 100
Level <- 100
gr <- exp_sizeFun(g_res$params_m, list_covs)
sv <- exp_sizeFun(s_res$params_m, list_covs)
sig_gr <- g_res$sigma
svlink<- s_res$family$linkinv

old <- unlist(lapply(ca, function(i) sum(g[i])))
new <- map_dbl(ca, ~ sum(g[.x]), g)


