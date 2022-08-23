document()
load_all()
library(profvis)
spe <- "Yggdrasil"
# Yggdrasil <- old_ipm2species(spe, path = "tests/testthat/testdata/", replicat = 1)
Yggdrasil <- old_ipm2species(spe, path = here(), replicat = 1)
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"

Forest <- forest(list(Yggdrasil))
# Forest$species$Yggdrasil$recruit_fun

load_all()
profvis({
set.seed(42)
res <- sim_deter_forest(Forest, tlim = 2000, equil_time = 1e4,
                 correction = "cut",
                 verbose = TRUE)
})



times <- as.numeric(sub("t", "", colnames(res)))
plot(times, res[grepl("BA",rownames(res)),], ylab = "Total BA", xlab = "time",
     cex = 0.1, type = "b", ylim = c(0, 100))
# plot(res[grepl("N",rownames(res)),])
text(700, res[grepl("BA",rownames(res)),ncol(res)],
     round(res[grepl("BA",rownames(res)),ncol(res)], 2), cex = .7, pos = 3)

set.seed(42)
Forest2 <- forest(list(Yggdrasil, Ents))
res2 <- sim_deter_forest(Forest2, tlim = 600, equil_time = 1e3,
                        correction = "cut",
                        verbose = TRUE)


value <- colSums(res2[grepl("BA",rownames(res2)),])
points(times, value, type = "b", cex = 0.1, col = "darkred")
text(700, value[ncol(res2)], round(value[ncol(res2)], 2), cex = .7, pos = 3)
points(times, res2[grepl("Ygg.*BA",rownames(res2)),],
       type = "b", cex = 0.1, col = "darkblue")
points(times, res2[grepl("Ent.*BA",rownames(res2)),],
       type = "b", cex = 0.1, col = "darkgreen")
legend("top", c("Yddgrasil solo", "Yggdrasil + Ents", "Yggdrasil", "Ents"),
       fill = c("black", "darkred", "darkblue", "darkgreen"))
text(700, res2[grepl("BA",rownames(res2)),ncol(res2)],
     round(res2[grepl("BA",rownames(res2)),ncol(res2)], 2), cex = .7, pos = 3)

# Benchmark to reduce time ####

get_ipm <- function(x, n){
    return(x$IPM$IPM[[n]])
}



lower_ba <- c(Yggdrasil  = 1)
higher_ba <- c(Yggdrasil  = 2)
t = 1
sim_BA <- c(1.111233, NA)

low_ba <- map2(Forest$species, lower_ba, get_ipm)
high_ba <- map2(Forest$species, higher_ba, get_ipm)

sim_ipm <- lapply(
    seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
        low_ba[[i]] * (1 - (floor(ba) - nipm[i])) + # TODO check if this is correct formula !
            high_ba[[i]] * ( floor(ba)  - nipm[i] )
    }, low_ba, high_ba, sim_BA[t], lower_ba
)

microbenchmark::microbenchmark(
    basic = lapply(
        seq_along(low_ba), function(i, low_ba, high_ba, ba, nipm){
            low_ba[[i]] * (1 - (floor(ba) - nipm[i])) +
                high_ba[[i]] * ( floor(ba)  - nipm[i] )
        }, low_ba, high_ba, sim_BA[t], lower_ba
    ),
    loop = {for (i in seq_along(low_ba)){
        sim_ipm[[i]] <- low_ba[[i]] * (1 - (floor(sim_BA[t]) - lower_ba[i])) +
            high_ba[[i]] * ( floor(sim_BA[t])  - lower_ba[i] )
    }}
)


## Test zone ####

expr(final <- exp(res) * SurfEch / 0.03 * distrib )

call2("<-", expr(final), call2("*", expr(exp(res)),
                               expr(SurfEch / 0.03 * distrib)))

