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
#> Time difference of 0.191 secs
```

The output is a single table with a row for each observation and a
column per variable. Variables are given at each time point plus
equilibrium state. They are defined below :

-   Distribution in the mesh of the species
-   Basal Area
-   Sum of the distribution
-   Distribution of harvest
-   Sum of harvest

``` r
sim1sp %>% mutate_if(is.character, as.factor) %>% summary()
#>       species       var            time            mesh           size          equil        
#>  Yggdrasil:3843   BAsp:  61   Min.   : 1.00   Min.   : 1.0   Min.   : 91.15   Mode :logical  
#>                   h   :1830   1st Qu.:16.00   1st Qu.: 8.0   1st Qu.:107.23   FALSE:3780     
#>                   H   :  61   Median :31.00   Median :15.5   Median :124.47   TRUE :63       
#>                   m   :1830   Mean   :31.03   Mean   :15.5   Mean   :124.47                  
#>                   N   :  61   3rd Qu.:46.00   3rd Qu.:23.0   3rd Qu.:141.70                  
#>                               Max.   :63.00   Max.   :30.0   Max.   :157.78                  
#>                                               NA's   :183    NA's   :183                     
#>      value          
#>  Min.   :  0.00000  
#>  1st Qu.:  0.03758  
#>  Median :  0.06730  
#>  Mean   :  6.10113  
#>  3rd Qu.:  6.72356  
#>  Max.   :223.80286  
#> 

sim1sp %>%
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
#> Time difference of 0.34 secs
```

The result is the same table.

``` r
sim2sp %>% mutate_if(is.character, as.factor) %>% summary()
#>       species       var            time            mesh           size          equil        
#>  Ents     :1953   BAsp:  62   Min.   : 1.00   Min.   : 1.0   Min.   : 91.15   Mode :logical  
#>  Yggdrasil:1953   h   :1860   1st Qu.: 8.00   1st Qu.: 8.0   1st Qu.:107.23   FALSE:3780     
#>                   H   :  62   Median :16.00   Median :15.5   Median :124.47   TRUE :126      
#>                   m   :1860   Mean   :17.13   Mean   :15.5   Mean   :124.47                  
#>                   N   :  62   3rd Qu.:24.00   3rd Qu.:23.0   3rd Qu.:141.70                  
#>                               Max.   :66.00   Max.   :30.0   Max.   :157.78                  
#>                                               NA's   :186    NA's   :186                     
#>      value          
#>  Min.   :  0.00000  
#>  1st Qu.:  0.03538  
#>  Median :  0.06254  
#>  Mean   :  5.64523  
#>  3rd Qu.:  6.53854  
#>  Max.   :218.81133  
#> 

sim2sp %>%
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
#> Time difference of 0.18 secs
```

``` r
sim1harv %>% mutate_if(is.character, as.factor) %>% summary()
#>       species       var            time            mesh           size          equil        
#>  Yggdrasil:3843   BAsp:  61   Min.   : 1.00   Min.   : 1.0   Min.   : 91.15   Mode :logical  
#>                   h   :1830   1st Qu.:16.00   1st Qu.: 8.0   1st Qu.:107.23   FALSE:3780     
#>                   H   :  61   Median :31.00   Median :15.5   Median :124.47   TRUE :63       
#>                   m   :1830   Mean   :30.98   Mean   :15.5   Mean   :124.47                  
#>                   N   :  61   3rd Qu.:46.00   3rd Qu.:23.0   3rd Qu.:141.70                  
#>                               Max.   :60.00   Max.   :30.0   Max.   :157.78                  
#>                                               NA's   :183    NA's   :183                     
#>      value        
#>  Min.   :  0.000  
#>  1st Qu.:  0.000  
#>  Median :  0.000  
#>  Mean   :  6.292  
#>  3rd Qu.:  7.163  
#>  Max.   :233.016  
#> 

sim1harv %>%
    filter(var %in% c("BAsp", "H"),  !equil) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) +
    facet_wrap(~ var, scales = "free_y") +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/harv_plot-1.png)<!-- -->

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
#> Time difference of 0.279 secs
```

``` r
sim2harv %>%
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
Yggdrasil_d5$info["species"] <- "Yggdrasil_d"
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
#> Time difference of 0.232 secs
```

Equilibrium BA should be really close ($\Delta_{BA} < 1$). N is expected
to increase with delay since delayed mesh cell with seeds are counted
in.

``` r
sim5d %>%
    rbind(sim1sp) %>%
    filter(var %in% c("BAsp", "N")) %>%
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linetype = "dashed", size = .3) +
    geom_point(size = .7) + ylab("BA") +
    facet_wrap(~ var, scales = "free_y") +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/delay_plot-1.png)<!-- -->

Despite a really close BA, the size distribution is different at
equilibrium.

*Note : Values below the redline does not count in BA computation.*

``` r
sim5d %>%
    mutate(mesh = mesh - n_delay) %>%
    rbind(sim1sp) %>%
    filter(var == "m", equil) %>%
    ggplot(aes(x = mesh, y = value, fill = species, group = desc(species))) +
    geom_area(color = "black", alpha = .3, size = .4, position = "identity") +
    geom_vline(xintercept = 1, color = "red") +
    annotate("text", x = 0, y = .1, label = "Minimal mesh in BA", size = 3,
             color = "red", angle = 90) +
    NULL
```

![](multisp_deter_sim_files/figure-gfm/delay_dist_plot-1.png)<!-- -->
