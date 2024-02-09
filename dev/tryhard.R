document()
load_all()
library(profvis)
library(ggplot2)
library(here)

# Benchmark to reduce time ####

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
