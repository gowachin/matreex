# matreex (development version)

* `targetBA` argument of `sim_deter_forest()` is now required per hectare and not for the corresponding `SurEch`.

* `sim_deter_forest()` output table changed the name for variable `m_i` (distribution of density by mesh i) for `n_i`. #8

* Correction of vignettes 

# matreex 0.2.0

## Breaking changes

* `tree_format()` will get deprecated and fully integrated in `sim_deter_forest()`.

## Documentation

* New functions `summary()` for `ipm` and `species` class.

* New vignette `Harvesting.Rmd`.

* New vignette `matreex.Rmd`.

## User interface

* `make_IPM()` is more flexible for climate input to simplify scripts. 
Single rows matrix and data.frame are now accepted. 
This replace `climate <- drop(as.matrix(subset(climate_species, ...)`.

* `species()` class constructor has now default values and should be preferred to `new_species()`.

# matreex 0.1.0

* Added a `NEWS.md` file to track changes to the package.

* Changed package name to `{matreex}`.

* Add disturbance model from Julien.

* Class `ipm` was modified and *invalidate* previous script output (if saved as
.Rds/.Rdata).

* Idea from Arnaud to precompute mu values to speed up in-simulation integration
and modify climate. This kind of simulations takes longer but is very flexible. 
New class `mu_gr`.

* Refactorisation of harvest process for Uneven and Even process. 
Uneven works with multi species.

* Refactorisation of integration process to gain time and allow a more general
function (added input choices).

* Multi species simulations works with 1 to n species in a general setup. 
Also created a generic output format usable in `{tidyverse}`.

* Add functions to load old IPMs and modify them in new class (`ipm`, `species`
, `forest`).

* Packaging functions from [Kunstler *et al* (2020)](https://doi.org/10.1111/1365-2745.13533) under initial name `{treeforce}`.
