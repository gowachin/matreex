document()
load_all()
library(profvis)

spe <- "Yggdrasil"
Yggdrasil <- old_ipm2species(spe, path = here(), replicat = 1,
                             harvest = def_harv)


load_all()
Yggdrasil$harvest_fun <- Uneven_harv
body(Yggdrasil$init_pop)[[length(body(Yggdrasil$init_pop)) - 1]] <- expr(res <- res * 20)
body(Yggdrasil$recruit_fun)[[length(body(Yggdrasil$recruit_fun)) - 2]] <- expr(distrib <- c(1, numeric(mesh - 1)))
Forest <- forest(
    list(
        delay(Yggdrasil, delay = 0)
    ),
    harv_rules = c(Pmax = 0.7, dBAmin = 3,
                   freq = 30, alpha = 1))
# Forest$species$Yggdrasil$recruit_fun
targetBA <- 30

# profvis({
set.seed(42)
res <- sim_deter_forest(Forest, tlim = 1000, equil_time = 1000,
                 correction = "cut", targetBA = targetBA,
                 verbose = TRUE)
# })


n <- 15
memor <- vector("list", n+1)
memor[[1]] <- tree_format(res)


for(d in 1:n){
    if(d %% 10 == 0){
        cat("\n========", d, "========\n")
    }
    set.seed(42)
    Forest <- forest(
        list(
            delay(Yggdrasil, delay = d)
        ),
        harv_rules = c(Pmax = 0.7, dBAmin = 3,
                       freq = 30, alpha = 1))
    res <- sim_deter_forest(Forest, tlim = 1000, equil_time = 1000,
                            correction = "cut", targetBA = targetBA,
                            verbose = TRUE)
    memor[[d+1]] <- tree_format(res)
}

memor <- map2(memor, 1:(n+1), ~ mutate(.x, delay = .y))
memor <- do.call("rbind", memor)


memor %>%
    filter(var == "m", equil, delay %in% c(1, 5, 10, 15)) %>%
    mutate(mesh = mesh - delay + 1) %>%
    ggplot(aes(x = mesh, y = value, color = factor(delay))) +
    geom_line(size = .4) +
    geom_vline(xintercept = 1, color = "red") +
    annotate("text", x = 0, y = .1, label = "Minimal mesh in BA", size = 3,
             color = "red", angle = 90) +
    NULL

memor %>%
    filter(var == "BAsp", ! equil, value != 0) %>%
    # dplyr::na_if(0) %>%
    ggplot(aes(x = time, y = value, color = factor(delay))) +
    geom_line(size = .4) +
    NULL

memor %>%
    filter(var == "N", equil) %>%
    ggplot(aes(x = delay, y = value, color = factor(delay))) +
    geom_point(size = .4) +
    # geom_vline(xintercept = 1, color = "red") +
    # annotate("text", x = 0, y = .1, label = "Minimal mesh in BA", size = 3,
    # color = "red", angle = 90) +
    NULL

memor %>%
    filter(var == "m", delay %in% c(1, 5, 15)) %>%
    dplyr::na_if(0) %>%
    ggplot(aes(x = time, y = mesh, fill = value)) +
    facet_wrap(~ delay) +
    geom_tile() +
    scale_fill_viridis_c(na.value="transparent") +
    theme_dark() +
    NULL

# Two species
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"
# body(Ents$init_pop)[[length(body(Ents$init_pop)) - 1]] <- expr(res <- res * 80)
Ents <- delay(Ents, 10)
load_all()
set.seed(42)
Forest2 <- forest(list(Yggdrasil, Ents),
                  harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                 freq = 10, alpha = 1))
res2 <- sim_deter_forest(Forest2, tlim = 1000, equil_time = 1000,
                        correction = "cut", targetBA = targetBA,
                        verbose = TRUE)

library(ggplot2)
library(viridis)
tmp <- tree_format(res2)

tmp %>%
    filter(! var %in%  c("m", "h")) %>%
    filter(time < 50) %>%
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

tmp %>%
    filter(var == "m", equil == T) %>%
    filter(value != 0) %>%
    ggplot(aes(x = value)) +
    geom_density() +
    facet_wrap(~ species) +
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
