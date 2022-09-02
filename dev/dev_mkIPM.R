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
   make_IPM_GL_2_i, mk_P_GL_2)


# in fun_mid_point
h <- (U - L) / m
mesh_x <- seq(L+h/2, U-h/2, length.out=m)
N_int <- sum((mesh_x - min(mesh_x)) < diag_tresh)
N_ini <- N_int+1
N_int <- 100
Level <- 100
gr <- g_res
sv <- s_res

microbenchmark::microbenchmark(
    init <- fun_mid_int(m, L, U, gr, sv, N_ini, N_int, list_covs, Level=100),
    edit <- fun_mid_int_stripe(m, L, U, gr, sv, N_ini, N_int, list_covs, Level=100),
    times =  10
)


source("dev/dev_old_mkipm.R")
tic()
set.seed(42)
init <- fun_mid_int(m, L, U, gr, sv, N_ini, N_int, list_covs, Level=100)
toc()
tic()
set.seed(42)
edit <- fun_mid_int_stripe(m, L, U, gr, sv, N_ini, N_int, list_covs, Level=100)
toc()
all.equal(init, edit)

gt <- c(3.11367626649636e-08, 2.82295872819631e-08, 2.56129586967654e-08,
        2.32559264224615e-08, 2.11310406647111e-08, 1.92139333796995e-08,
        1.74829522105495e-08, 1.59188402834669e-08)
cat <- c(0L, 0L, 1L, 1L, 1L, 1L, 2L, 2L)
h <- 2.29771428571429
Level <- 100
fuu <- function(gt, cat, h, Level){
    old <- data.frame(gt = gt, cat = cat) %>%
        group_by(cat) %>% summarise(P = sum(gt) * h / Level)
    old$P
}

foo <- function(gt, cat, h, Level){
    n <- unique(cat)
    res <- numeric(length = length(n))
    for(i in seq_along(n)){
        res[i] <- sum(gt[cat == n[i]]) * h /Level
    }
    res
}

fii <- function(gt, cat, h, Level){
    res <- split(gt, cat)
    res <-  unlist(lapply(res, sum)) * h / Level
}

microbenchmark::microbenchmark(
    old = fuu(g, ca, h, Level),
    new = foo(g, ca, h, Level),
    best = fii(g, ca, h, Level),
    times = 100
)



