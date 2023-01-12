Harvesting models
================

This vignette illustrate how to use harvesting scenarii with `{matreex}`
package. An harvesting scenario is a simulation modification that will
trigger some additional death in defined size distribution. This
distribution is differenciated from base mortality and is saved and
exported for each simulation.

This document will skip the usage of basic function that are shown in [a
previous introduction vignette](matreex.html).

# Simulations input

## Define a species

The first step of IPM integration. This part is common with basic usage
of the package so nothing is very important here.

**Please keep in mind this computation is intensive and may take few
minutes !**

``` r
library(matreex)
library(dplyr)
library(ggplot2)

# Load fitted model for a species
# fit_species # list of all species in dataset
data("fit_Picea_abies")

# Load associated climate
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)
# see ?climate_species to understand the filtering of N.
climate
#>        sgdd       wai        sgddb      waib      wai2   sgdd2      PC1        PC2 N       SDM
#> 62 1444.667 0.4519387 0.0006922012 0.6887343 0.2042486 2087062 1.671498 0.02602064 2 0.6760556

Picea_ipm <- make_IPM(
    species = "Picea_abies", 
    climate = climate, 
    fit = fit_Picea_abies,
    clim_lab = "optimum clim", 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:50, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)
#> Launching integration loop
#> GL integration occur on 32 cells
#> midbin integration occur on 25 cells
#> Loop done.
#> Time difference of 25.8 secs
```

Harvesting rules have different scales : some rules are defined for each
species, others at the forest level and lastly input parameters when
launching simulation. These rules can always be defined but will only be
used when the correct module is triggered by simulations. The table
below summarise the required parameters for each scenario. Each scenario
have an example in the following document.

| Harvesting scenario       | Species function                     | Species parameters | Forest parameters                | Simulation parameters       |
|---------------------------|--------------------------------------|--------------------|----------------------------------|-----------------------------|
| default                   | `def_harv()` <br> or custom by user. | $\emptyset$        | `harv_rules["freq"]`             | $\emptyset$                 |
| Uneven                    | `Uneven_harv()`                      | `harv_lim`         | `harv_rules`                     | `targetBA`                  |
| Even <br> *in active dev* | `Even_harv()`                        | `rdi_coef`         | `harv_rules` <!--`FinalHarvT`--> | `targetRDI` <br> `targetKg` |

Short example : When using an Uneven scenario, each species needs to use
the `Uneven_harv()` function with `harv_lim` parameters. The forest will
need to have `harv_rules` parameters. The `sim_deter_forest()` require
`targetBA` argument. If `targetRDI` is defined, it will not be used.

In multi-specific simulation, it is possible to combine Uneven or Even
managed species with default species.

| species 1 | species 2 | Harvesting scenario to use |
|-----------|-----------|----------------------------|
| default   | default   | default                    |
| Uneven    | default   | Uneven                     |
| Even      | default   | Even                       |
| Uneven    | Even      | **Not possible !**         |

# Default scenario

## Presentation

The default scenario is based on Kunstler et al.
([2021](#ref-kunstler2021)). This mean there is a constant harvest rate
triggered each year. This harvest rate is uniform on the distribution.
This is the scenario used when no modification is given.

This rate is coded in `def_harv()` function as shown below and the
frequency (each year) is given in `harv_rules["freq"]`. Note that the
rate only impact all size classes except classes used in delay (with
`* (ct > 0)`).

``` r
def_harv
#> function(x, species, ...){
#> 
#>     dots <- list(...)
#>     ct <- dots$ct
#> 
#>     rate <- 0.006 * (ct > 0)
#>     return(x * rate)
#> }
#> <environment: namespace:matreex>
Picea_sp <- species(IPM = Picea_ipm, init_pop = def_initBA(30))
Picea_for <- forest(species = list(Picea = Picea_sp), 
                    harv_rules = c(Pmax = 0.25, dBAmin = 3, 
                                   freq = 1, alpha = 1))
```

With this scenario, the simulation can be launched without any
additional parameter.

``` r
set.seed(42) # The seed is here for initial population random functions.
Picea_sim <- sim_deter_forest(
    Picea_for,
    tlim = 200,
    equil_time = 200, equil_dist = 10, equil_diff = 1,
    harvest = "default", # this is the default value but we write it.
    SurfEch = 0.03,
    verbose = TRUE
)
#> Starting while loop. Maximum t = 200
#> Simulation ended after time 200
#> BA stabilized at 45.16 with diff of 0.69 at time 200
#> Time difference of 0.821 secs
```

Once the simulation is done, we can extract the basal area and the
number of individual at each step of the simulation. The new variable we
can look now is $H$, the sum of harvest distribution at each step. This
distribution is also exported with an harvest value for each mesh cell
$h_i$.

In this case, $H$ is correlated with $N$ since it’s a constant
percentage not linked with a size distribution. The first step being the
initialization step, it’s normal to have no harvest.

``` r
Picea_sim  %>%
    dplyr::filter(var %in% c("BAsp", "N", "H"), ! equil) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .2) + geom_point(size = 0.4) 
```

![](Harvesting_files/figure-gfm/sp1plot-1.png)<!-- -->

## Modulation

This section will just illustrate variation of the default scenario.
First we modify the frequency of the harvest. When the harvest is not
triggered, the value returned is 0.

``` r
set.seed(42) # The seed is here for initial population random functions.
Picea_sim_f20 <- sim_deter_forest(
    forest(species = list(Picea = Picea_sp), 
                      harv_rules = c(Pmax = 0.25, dBAmin = 3, 
                                     freq = 20, alpha = 1)),
    tlim = 50,
    equil_time = 50, equil_dist = 10, equil_diff = 1,
    harvest = "default", 
    SurfEch = 0.03,
    verbose = TRUE
)
#> Starting while loop. Maximum t = 50
#> Simulation ended after time 50
#> BA stabilized at 30.14 with diff of 0.02 at time 50
#> Time difference of 0.27 secs
Picea_sim_f20  %>%
    dplyr::filter(var %in% c("BAsp", "N", "H"), ! equil) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .2) + geom_point(size = 0.4) 
```

![](Harvesting_files/figure-gfm/def_freq-1.png)<!-- -->

A more advanced step consist in modification of the harvest function for
more custom effects. Obviously, this kind of modification is more prone
to error so don’t hesitate to contact `{matreex}` maintainer. For
example, this is function where we multiply a constant rate with the
mesh, meaning that the larger the tree get, the more we its size class.

``` r
Picea_harv <- Picea_sp
Picea_harv$harvest_fun <- function(x, species, ...){

    dots <- list(...)
    ct <- dots$ct

    rate <- 6e-4 * (ct > 0) * species$IPM$mesh 
    return(x * rate)
}

set.seed(42) # The seed is here for initial population random functions.
Picea_sim_f20 <- sim_deter_forest(
    forest(species = list(Picea = Picea_harv), 
                      harv_rules = c(Pmax = 0.25, dBAmin = 3, 
                                     freq = 20, alpha = 1)),
    tlim = 250,
    equil_time = 250, equil_dist = 10, equil_diff = 1,
    harvest = "default",
    SurfEch = 0.03,
    verbose = TRUE
)
#> Starting while loop. Maximum t = 250
#> Simulation ended after time 250
#> BA stabilized at 25.86 with diff of 5.25 at time 250
#> Time difference of 1.07 secs
Picea_sim_f20  %>%
    dplyr::filter(var %in% c("BAsp", "N", "H"), ! equil) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .2) + geom_point(size = 0.4) 
```

![](Harvesting_files/figure-gfm/edit_def_harv-1.png)<!-- -->

# Uneven scenario

## Theory

Uneven harvest scenario consist in harvesting tree depending on the
basal area of the stand and the distribution of the tree sizes. This
should lead to stands with uneven size distribution.

### Monospecific case

#### Harvest proportion

We note $P_{cut}$ the global harvest proportion which will determine the
amount of basal area to be cut and harvested.

$$
P_{cut} = \left\{
 \begin{array}{ll}
    0 & if (BA_{stand} - BA_{target}) < \Delta BA_{min}\\
    min(\frac{BA_{stand}-BA_{target}}{BA_{stand}}, P_{max}) & if (BA_{stand} - BA_{target}) \geq \Delta BA_{min}
 \end{array}
 \right\}
$$

For example, $\Delta BA_{min} = 3 m^2ha^{-1}$, $BA_{target} = 20, 25$ or
$30 m^2ha^{-1}$ (depending on species) and $P_{max} = 0.25$.

Note that stand basal area $BA_{stand}$ is computed only considering
trees with a dbh above $d_{th}$ (see below).

#### Harvest curve

Each tree harvest probability only depends on its diameter ($d$). There
is a minimum diameter of harvest ($d_{th}$), harvest probability then
increases with diameter until $d_{ha}$, the harvesting diameter after
which harvest probability is high and constant.

We therefore considered the harvesting function (which associates a dbh
to an harvesting probability) $$
h(d) = \left\lbrace
 \begin{array}{ll}
    0 & \text{if } d < d_{th} \\
    h_{max} (\frac{d - d_{th}}{d_{ha} - d_{th}})^{k} & \text{if } d_{th}\leq d < d_{ha} \\
    h_{max} & \text{if } d \geq d_{ha} 
 \end{array}
 \right\rbrace
$$

The parameter $h_{max}$ can be tuned so that the probability for a large
tree to be harvested approaches 1 after several harvesting operations:
$h_{max} = 1 - \sqrt[n]{1-p}$ with $n$ the number of harvesting
operations. Parameter $k$ defines how much the harvest preferentially
selects large trees.

![Harvest curve example, $d_{th} = 17.5cm$, $d_{ha}=57.5cm$,
$h_{max}=0.8$,
$k=2$.](Harvesting_files/figure-gfm/harvest_curve_plot-1.png)

#### Harvesting algorithm

The algorithm only target tree contributing to $BA_{stand}$, that is the
trees above $d_{th}$

Given $\phi(x)$ the density function of diameters, the basal area
harvested is

$$ 
\begin{array}{ll}
BA_{harv} & = \pi/4\int_{d_{th}}^{d_{max}}x^2 h(x) \phi(x)dx \\
          & = \pi/4\int_{d_{th}}^{d_{ha}}x^2 h(x) \phi(x)dx + h_{max} \pi/4\int_{d_{ha}}^{d_{max}}x^2 \phi(x)dx\\
      & = BA_{th} + BA_{ha}
\end{array}
$$

$\bullet$ If $BA_{ha} >= P_{cut} \times BA_{stand}$, there is enough
large trees (diameter above $d_{ha}$) so that the harvest will only
concern large trees: $BA_{th} = 0$.

We then have two different strategies : we can either adapt $h_{max}$ or
$d_{ha}$ so that:
$$  BA_{ha} = h_{max} \pi/4 \int_{d_{ha}}^{d_{max}}x^2 \phi(x)dx = P_{cut} \times BA_{stand}$$
The `{matreex}` package will cut the larger trees until
$P_{cut} \times BA_{ha} - targetBA <= 0$, since $h_{max}$ and $d_{ha}$
are values set by the user and won’t be modified. This is equivalent to
find an harvest diameter $d_t$ such as $d_{ha} < d_t < d_{max}$.

$\bullet$ If $BA_{ha} < P_{cut} \times BA_{stand}$, we first harvest
$BA_{ha}$ and then compute $k$ such as:
$BA_{th} = P_{cut} \times BA_{stand} - BA_{ha}$.

### Multispecific case

<!-- #### Abundance-based preference -->

As in the monospecific case, we define the global harvest rate
$P_{cut} = \frac{BA_{harv}}{BA_{stand}}$.

Here, $BA_{stand}$ is divided between $s$ species:
$BA_{stand} = \sum_{i=1}^{s} BA_{stand, i}$ and
$BA_{harv} = \sum_{i=1}^{s} BA_{harv, i}$

We note $p_i = \frac{BA_{stand, i}}{BA_{stand}}$ and
$P_{cut, i} = \frac{BA_{stand, i} - BA_{harv, i}}{BA_{stand, i}} = f(p_i) * P_{cut}$

We suppose that harvesting rate increases with abundance (we harvest
preferentially trees with the highest proportion), which means $f$ is an
increasing function.

By definition,

$$
\begin{array}{ll}
BA_{harv} & = \sum_{i=1}^{s} BA_{harv, i} = \sum_{i=1}^{s} BA_{stand, i} * (1 - P_{cut, i}) \\
          & = BA_{stand} \sum_{i=1}^{s} p_i * (1 - f(P_i) P_{cut}) \\
      & = BA_{stand} (1 - P_{cut}) \\
\end{array}
$$

So that we have the constraint on $f$:
$\sum_{i=1}^{s} p_i (1-f(p_i)P_{cut}) = 1 - P_{cut} \sum_{i=1}^{s} p_i f(p_i) = 1-P_{cut}$

which is equivalent to $\sum_{i=1}^{s} p_i f(p_i) = 1$

The case $f(p_i) = 1$ works, which leads to $P_{cut,i} = P_{cut}$. In
that case the harvest rate is the same for every species $i$.

More broadly, \$ \> 0\$

$$f(p_i) = \frac{p_i^{\alpha - 1}}{\sum_{i=1}^{s} p_i ^{\alpha}}$$

For $\alpha = 2$, we for example have

$$f(p_i) =\frac{p_i}{\sum_{i=1}^{s} p_i ^2} $$ To sum up $\alpha = 1$,
an abundant species will be more harvested, in effort to balance the
species. For $\alpha < 1$, the abundant species will be less harvested.
Finally, $\alpha > 1$ the most abundant species will be cut more, in a
greater proportion than that allowing to tend towards an equilibrium.

<!--
#### Favoured species

In some case, we may want to favour some species compared to others.
We note $Q$ the $q$ species we want to favour, and $P_Q=\sum_{i=1}^{q} p_i$
We first compute the harvest rate $H_Q$ and $H_{1-Q}$ for respectively the favoured/other species.

If $P_Q \geq 0.5$, we take $H_Q = H_{1-Q} = H$ (the species to be favoured are already dominant).


If $P_Q < 0.5$, we compute $H_Q=f(P_Q) H$ and $H_{1-Q}=f(1-P_Q) H$ with $\alpha > 1$.
By definition, we will get $H_Q \leq H$.

We then apply for each species $i$ the harvest rate $H_i = H_Q$ or $H_i=H_{1-Q}$ depending on which group it belongs to.
-->

## Examples

All the parameters described above are input either in the `species()`,
`forest()` or `sim_deter_forest()` functions. Additional parameter
`dBAmin` is the difference between `BA` and `targetBA` under which an
harvesting will not be triggered.

``` r
Picea_Uneven <- species(IPM = Picea_ipm, init_pop = def_initBA(30), 
                        harvest_fun = Uneven_harv,
                        harv_lim = c(dth = 175, dha = 575, hmax = 1))
Picea_for_Uneven <- forest(species = list(Picea = Picea_Uneven), 
                      harv_rules = c(Pmax = 0.25, dBAmin = 3, 
                                     freq = 5, alpha = 1))
```

``` r
set.seed(42) # The seed is here for initial population random functions.
Picea_sim_f20 <- sim_deter_forest(
    Picea_for_Uneven,
    tlim = 200,
    equil_time = 200, equil_dist = 10, equil_diff = 1,
    harvest = "Uneven", targetBA = 30, # We change the harvest and set targetBA.
    SurfEch = 0.03,
    verbose = TRUE
)
#> Starting while loop. Maximum t = 200
#> Simulation ended after time 200
#> BA stabilized at 31.86 with diff of 2.71 at time 200
#> Time difference of 0.948 secs
Picea_sim_f20  %>%
    dplyr::filter(var %in% c("BAsp", "N", "H"), ! equil) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .2) + geom_point(size = 0.4) 
```

![](Harvesting_files/figure-gfm/Uneven_sim-1.png)<!-- -->

We notice that the basal area obtained by the simulation is higher than
the targeted one. This can be explained by the fact that the cutting
calculation is done on $BA_{stand}$, which does not take into account
individuals smaller than $d_{th}$.

<!-- # Even scenario -->
<!-- $$ -->
<!-- dg = \sqrt{\frac{\sum_{i = 0}^n d_i^2 x_i}{ \sum_{i = 0}^n x_i }} \\ -->
<!-- dg_{cut} = \sqrt{\frac{\sum_{i = 0}^n d_i^2 x_i Pc_i }{ \sum_{i = 0}^n x_i Pc_i }} -->
<!-- $$ -->
<!-- mesh in cm. -->
<!-- $$ -->
<!-- RDI = \frac{ \sum_{i = 0}^n x_i }{ e^{ RDI_{int} + RDI_{slope} \times \frac{ \log( \frac{ \sum_{i = 0}^n mesh_i^2  x_i }{ \sum_{i = 0}^n x_i} ) }{2}   } } -->
<!-- $$ -->

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-kunstler2021" class="csl-entry">

Kunstler, Georges, Arnaud Guyennon, Sophia Ratcliffe, Nadja Rüger,
Paloma Ruiz-Benito, Dylan Z. Childs, Jonas Dahlgren, et al. 2021.
“Demographic Performance of European Tree Species at Their Hot and Cold
Climatic Edges.” *Journal of Ecology* 109 (2): 1041–54.
<https://doi.org/10.1111/1365-2745.13533>.

</div>

</div>