Running multi species deterministic simulations
================

``` r
library(matreex)
library(ggplot2)
library(dplyr)
library(magrittr)
```

# Introduction

Simulation run on the Forest class object, that contains all elements.
It’s composed of *n* species. Each species is defined by multiples
functions :

-   init_pop is the function to initiate the population at time 1.
-   recruit_fun is the function that depend on the total basal area.
-   harvest_fun is the function that return the harvested population
    that should be subtracted to the population at time t.
-   IPM is the list of integrated matrices for several basal area.

# Loading previous work IPM

Functions are defined to load the previous IPM saved as *.Rds*. The
files must located in specific directory format which is :

**`path/output/species/IPM_Clim_X.Rds`**

Once the files are in place, the code below allow to create the species
of interest. One last argument is *replicat* which select one of the 100
simulations of previous work. **A mean IPM should replace this later**.

Default functions of initialization and harvest are used. You can write
your own functions as long as they use the same arguments. Try to write
function that don’t initialize negative populations.

Example data set is a cropped IPM with a mesh of length 30 with BA from
1 to 10.

``` r
spe <- "Yggdrasil" # Example dataset given with the package
# Warning, loading an old IPM takes around 10s
Yggdrasil <- old_ipm2species(
    spe, climatic = 1, replicat = 1,
    path = system.file("extdata", package = "matreex", mustWork = TRUE),
)
def_init
#> function(mesh, SurfEch = 0.03) {
#>     ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
#>     ini <- exp(runif(1, -.005, .005) * mesh)
#>     alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
#>     while(all(alea)){ # because god knows it's fucking possible that alea is
#>                       # all FALSE and it will return NaN
#>         alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
#>     }
#>     ini[alea] <- 0
#>     res <- as.numeric(ini / sum(ct * ini) )
#>     res <- res + 1e-4 # HACK to limit falling in floating point trap !
#>                       # also line to add BA later if needed
#>     return(res)
#> }
#> <environment: namespace:matreex>
```

A forest is just a list of the species. **Later it will also require
harvest rules with target BA and species**

``` r
Forest <- forest(list(Yggdrasil))
```

# Simulation

Simulations are done with the function `sim_deter_forest`. Such
simulations are deterministic once the initial population is defined. As
the species *Yggdrasil* we have defined use random generation for the
initial population, we will use a seed before each simulation.

## Single species

A simulation is mostly defined by the time limitation (*tlim*) and an
equilibrium time (*equil_time*).

``` r
set.seed(42)
sim1sp <- sim_deter_forest(Forest, tlim = 60, equil_time = 1e3,
                        correction = "cut", equil_dist = 50,
                        verbose = TRUE)
#> apply a IPM cut correction
#> Starting while loop. Maximum t = 1000
#> Simulation ended after time 63
#> BA stabilized at 2.38 with diff of 0.98 at time 63
#> Time difference of 0.241 secs
```

The output is a single table with time in column and different variables
in rows. The variables are given at each time point plus equilibrium
state. They are defined below :

-   Distribution in the mesh of the species
-   Basal Area
-   Sum of the distribution
-   Distribution of harvest
-   Sum of harvest

*Only few rows and columns are displayed below*

``` r
m <- length(Forest$species$Yggdrasil$IPM$mesh)
sim1sp[c(1:2, m:(m+3), (2 *m + 2):(2 *m + 3)), c(1:3,30:31)]
#>                          t1           t2          t3          t30          t31
#> Yggdrasil.m1    0.003333333 7.747297e+00  8.84561167 7.248014e+00 7.269409e+00
#> Yggdrasil.m2    0.003333333 7.749060e+00 12.94675813 1.190738e+01 1.194154e+01
#> Yggdrasil.m30   0.003333333 0.000000e+00  0.00000000 0.000000e+00 0.000000e+00
#> Yggdrasil.BAsp  1.001247781 1.130797e+00  1.18682383 2.467116e+00 2.444072e+00
#> Yggdrasil.N    73.344115088 8.685144e+01 95.88504142 2.111872e+02 2.097068e+02
#> Yggdrasil.h1    0.000000000 3.520878e-06  0.00821321 7.956546e-03 7.975125e-03
#> Yggdrasil.h30   0.000000000 0.000000e+00  0.00000000 0.000000e+00 0.000000e+00
#> Yggdrasil.H     0.000000000 4.307325e-01  0.48842130 1.203184e+00 1.194027e+00
times <- as.numeric(sub(".*t", "", colnames(sim1sp)))

tree_format(sim1sp) %>%
    filter(var == "BAsp") %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) + ylab("BA") +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/single_print-1.png)<!-- -->

## Multiple species

Below, I create a mock second species to illustrate how to combine
species in simulation and recreate the same plot for total basal area.

We expect the same equilibrium result (minus the initial population
random effect), since the species dynamics are similar.

``` r
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"

Forest2 <- forest(species = list(Yggdrasil = Yggdrasil,
                                 Ents = Ents))

set.seed(42)
sim2sp <- sim_deter_forest(Forest2, tlim = 30, equil_time = 1e3,
                        correction = "cut", equil_dist = 50,
                        verbose = TRUE)
#> apply a IPM cut correction
#> Starting while loop. Maximum t = 1000
#> Simulation ended after time 66
#> BA stabilized at 4.65 with diff of 0.87 at time 66
#> Time difference of 0.376 secs
```

The result is the same table with each species tables grouped one under
the other. Since the species names is pasted with the variable name in
row names, one can extract variables of interest.

*Only few rows and columns are displayed below*

    #>                          t1           t2          t30        eqt66
    #> Yggdrasil.m1    0.003333333 7.548911e+00 6.950370e+00 7.073081e+00
    #> Yggdrasil.m2    0.003333333 7.550683e+00 1.152393e+01 1.173678e+01
    #> Yggdrasil.m30   0.003333333 0.000000e+00 0.000000e+00 0.000000e+00
    #> Yggdrasil.BAsp  1.001247781 1.128276e+00 2.425893e+00 2.322848e+00
    #> Yggdrasil.N    73.344115088 8.650693e+01 2.072597e+02 2.019230e+02
    #> Yggdrasil.h1    0.000000000 3.620717e-06 8.127348e-03 8.257006e-03
    #> Yggdrasil.h30   0.000000000 0.000000e+00 0.000000e+00 0.000000e+00
    #> Yggdrasil.H     0.000000000 4.310481e-01 1.183411e+00 1.149976e+00
    #> Ents.m1         0.003333333 7.548911e+00 7.002202e+00 7.073474e+00
    #> Ents.m2        13.470737078 9.865542e+00 1.161425e+01 1.173750e+01
    #> Ents.m30        0.003333333 0.000000e+00 0.000000e+00 0.000000e+00
    #> Ents.BAsp       1.001247781 1.014490e+00 2.365301e+00 2.322459e+00
    #> Ents.N         91.963362618 9.921262e+01 2.037742e+02 2.019034e+02
    #> Ents.h1         0.000000000 3.620717e-06 8.194788e-03 8.257563e-03
    #> Ents.h30        0.000000000 0.000000e+00 0.000000e+00 0.000000e+00
    #> Ents.H          0.000000000 5.077425e-01 1.161881e+00 1.149854e+00

``` r
tree_format(sim2sp) %>%
    filter(var == "BAsp") %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) + ylab("BA") +
    stat_summary(fun = sum, color = 'black', 
                 geom ='line', linetype = "dashed", size = .3) +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/two_plot-1.png)<!-- -->

## Changing harvesting rules.

Once a species is defined, you can modify it and say use a different
harvest function. For example, the default function apply a constant
harvesting rate of 0.6 percent per year.

``` r
Yggdrasil$harvest_fun
#> function(x, species, ...){
#> 
#>     dots <- list(...)
#>     ct <- dots$ct
#> 
#>     rate <- 0.006 * (ct > 0)
#>     return(x * rate)
#> }
#> <environment: namespace:matreex>
```

We can use another function that set an Uneven harvest a specified
intervals. To modify this interval, we need to edit the harvest rules
(*harv_rules*) when creating a Forest

``` r
Yggdrasil$harvest_fun <- Uneven_harv

# Because this example is a reduced mesh, I need to modify the harv_lim of the species
Yggdrasil$harv_lim <- c(dth = 95, dha = 110, hmax = 1)
Forest_harv <- forest(list(Yggdrasil),
                      harv_rules = c(Pmax = 1, dBAmin = 0.2, freq = 10, alpha = 1))
targetBA <- 2
set.seed(42)
sim1harv <- sim_deter_forest(Forest_harv, tlim = 60, 
                        equil_time = 60, equil_dist = 5,
                        harvest = "Uneven", targetBA = targetBA,
                        correction = "cut", verbose = TRUE)
#> apply a IPM cut correction
#> Starting while loop. Maximum t = 60
#> Simulation ended after time 60
#> BA stabilized at 2.19 with diff of 0.32 at time 60
#> Time difference of 0.219 secs
```

    #>                          t1        t2        t3        t30        t31
    #> Yggdrasil.m1    0.003333333  7.747301  8.842460   7.131691   7.510876
    #> Yggdrasil.m2    0.003333333  7.749075 12.968488  11.756685  12.129711
    #> Yggdrasil.m30   0.003333333  0.000000  0.000000   0.000000   0.000000
    #> Yggdrasil.BAsp  1.001247781  1.136997  1.199171   2.137105   2.329908
    #> Yggdrasil.N    73.344115088 87.282173 96.749186 195.133715 206.155726
    #> Yggdrasil.h1    0.000000000  0.000000  0.000000   0.000000   0.000000
    #> Yggdrasil.h30   0.000000000  0.000000  0.000000   0.000000   0.000000
    #> Yggdrasil.H     0.000000000  0.000000  0.000000  26.438140   0.000000

``` r
tree_format(sim1harv) %>%
    filter(var %in% c("BAsp", "H"),  !equil) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) +
    facet_wrap(~ var, scales = "free_y") +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/harv_plot-1.png)<!-- -->

<!-- Here we observe that harvest is not enough to reach targetBA. This could be explained by a restriction in tree size target. -->

The uneven harvest is possible with multiple species. Multiple
parameters will scale the effect of harvest :

-   Pmax, the maximum percentage of harvest possible
-   freq, the frequency of harvest
-   alpha, which will modulate the distribution of the cut on the
    species according to their abundances. When alpha = 1, the species
    are cut according to their abundance, which tends to balance the
    ratio of species. Alpha \> 1 will overcut the most present species,
    which can lead to unstable equilibrium. Alpha \< 1 will undercut the
    abundant species as well.

``` r
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"

Forest_harv2 <- forest(species = list(Yggdrasil = Yggdrasil,
                                      Ents = Ents),
                       harv_rules = c(Pmax = 1, dBAmin = 0.2, freq = 5, alpha = 1))

set.seed(42)
sim2harv <- sim_deter_forest(Forest_harv2, tlim = 60, 
                        equil_time = 60, equil_dist = 5,
                        harvest = "Uneven", targetBA = targetBA,
                        correction = "cut", verbose = TRUE)
#> apply a IPM cut correction
#> Starting while loop. Maximum t = 60
#> Simulation ended after time 60
#> BA stabilized at 2.36 with diff of 1.36 at time 60
#> Time difference of 0.331 secs
```

``` r
tree_format(sim2harv) %>%
    filter(var == "BAsp", !equil) %>% 
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) + ylab("BA") +
    stat_summary(fun = sum, color = 'black', na.rm = TRUE,
                 geom ='line', linetype = "dashed", size = .3) +
    geom_hline(yintercept = targetBA, linetype = "dotted")+
    NULL
```

![](multisp_deter_sim_files/figure-gfm/harv_two_plot-1.png)<!-- -->

## Recruitment delay

We can modify a species to add more delay for recruitment of new
individuals.

``` r
n_delay <- 5

Yggdrasil$harvest_fun <- def_harv
Yggdrasil_d5 <- delay(Yggdrasil, n_delay)
Yggdrasil_d5$info["species"] <- "Yggdrasil_d5"
Forest_delay <- forest(list(Yggdrasil_d5))
```

Simulation doesn’t change anything, the delay is only defined at the IPM
level.

``` r
set.seed(42)
sim5d <- sim_deter_forest(Forest_delay, tlim = 60, equil_time = 1e3,
                        correction = "cut", equil_dist = 50,
                        verbose = TRUE)
#> apply a IPM cut correction
#> Starting while loop. Maximum t = 1000
#> Simulation ended after time 69
#> BA stabilized at 2.39 with diff of 0.93 at time 69
#> Time difference of 0.236 secs
```

Equilibrium BA should be really close ($\Delta_{BA} < 1$). N is expected
to increase with delay since delayed mesh cell with seeds are counted
in.

``` r
tree_format(sim5d) %>%
    rbind(tree_format(sim1sp)) %>%
    filter(var %in% c("BAsp", "N")) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) + ylab("BA") +
    facet_wrap(~ var, scales = "free_y") +
    NULL
#> Warning in data.frame(species = sub(pattern, "\\1", var, perl = TRUE), var =
#> sub(pattern, : NAs introduits lors de la conversion automatique
```

![](multisp_deter_sim_files/figure-gfm/delay_plot-1.png)<!-- -->

Despite a really close BA, the size distribution is different at
equilibrium.

*Note : Values below the redline does not count in BA computation.*

``` r
tree_format(sim5d) %>%
    mutate(mesh = mesh - n_delay) %>%
    rbind(tree_format(sim1sp)) %>%
    filter(var == "m", equil) %>%
    ggplot(aes(x = mesh, y = value, fill = species, group = desc(species))) +
    geom_area(color = "black", alpha = .3, size = .4, position = "identity") +
    geom_vline(xintercept = 1, color = "red") +
    annotate("text", x = 0, y = .1, label = "Minimal mesh in BA", size = 3,
             color = "red", angle = 90) +
    NULL
#> Warning in data.frame(species = sub(pattern, "\\1", var, perl = TRUE), var =
#> sub(pattern, : NAs introduits lors de la conversion automatique
```

![](multisp_deter_sim_files/figure-gfm/delay_dist_plot-1.png)<!-- -->
