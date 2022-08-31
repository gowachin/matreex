document()
load_all()
library(profvis)

spe <- "Yggdrasil"
Yggdrasil <- old_ipm2species(spe, path = here(), replicat = 1,
                             harvest = def_harv)

load_all()
Yggdrasil$harvest_fun <- Uneven_harv
Forest <- forest(list(Yggdrasil),
                 harv_rules = c(Pmax = 0.7, dBAmin = 3,
                                freq = 30, alpha = 1))
# Forest$species$Yggdrasil$recruit_fun
targetBA <- 30

# profvis({
set.seed(42)
res <- sim_deter_forest(Forest, tlim = 600, equil_time = 600,
                 correction = "cut", targetBA = targetBA,
                 verbose = TRUE)
# })


times <- as.numeric(sub("t", "", colnames(res)))
plot(times, res[grepl("BA",rownames(res)),], ylab = "Total BA", xlab = "time",
     cex = 0.1, ylim = c(0, 100), type = "b")
text(700, res[grepl("BA",rownames(res)),ncol(res)],
     round(res[grepl("BA",rownames(res)),ncol(res)], 2), cex = .7, pos = 3)
abline(h = targetBA, lty = 3, col = "red")

# Two species
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"
Ents$init_pop <- function(mesh, SurfEch = 0.03) {
    ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
    ini <- exp(runif(1, -.005, .005) * mesh)
    alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    while(all(alea)){ # because god knows it's fucking possible.
        # and it will return NaN
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    }
    ini[alea] <- 0
    res <- as.numeric(ini / sum(ct * ini) )
    res <- res + 1e-10 # HACK to limit falling in floating point trap !
    res <- res * 80 # modif pour avoir un BA d'origine de 80
    return(res)
}
load_all()
set.seed(42)
Forest2 <- forest(list(Yggdrasil, Ents),
                  harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                 freq = 10, alpha = 1))
res2 <- sim_deter_forest(Forest2, tlim = 600, equil_time = 600,
                        correction = "cut", targetBA = targetBA,
                        verbose = TRUE)

library(ggplot2)
library(viridis)
tmp <- tree_format(res2)

tmp %>%
    filter(! var %in%  c("m", "h")) %>%
    filter(value > 0) %>% # For H
    ggplot(aes(x = time, y = value, color = species)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .5) +
    NULL

tmp %>%
    filter(var == "h") %>%
    dplyr::na_if(0) %>%
    ggplot(aes(x = time, y = mesh, fill = value)) +
    facet_wrap(~ species) +
    geom_tile() +
    scale_fill_viridis_c(na.value="transparent") +
    theme_dark() +
    NULL

tmp %>%
    filter(var == "m") %>%
    dplyr::na_if(0) %>%
    ggplot(aes(x = time, y = mesh, fill = value)) +
    facet_wrap(~ species) +
    geom_tile() +
    scale_fill_viridis_c(na.value="transparent") +
    theme_dark() +
    NULL


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

