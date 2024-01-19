# matreex (development version)

# matreex 0.4.0 dev

## News

* Mutlispecific even harvesting functionnality. #24

* `sim_rdikg()` function to compute RDI and Kg after simulations. #24

* Disturbance use a mixture effect in case of biotic disturbance. #17

* Favoured species for Uneven harvest model.

* Spatial effects are under development. #15

## Breaking changes

* `mu_sgr` class now use the delay as normal IPMs. Previous object may not work.

* `Even_harv()` function is not used, `getPcutEven()` is called directly.

## Documentation

* New vignette `Disturbance.Rmd` after completion of the model by Julien and Jasper. #17

* Completions of `Harvesting.Rmd` after completion of the model by Jasper. #24

* New files for metadat and HAL upload in prevision.

## Datasets

* Update of the disturbance coefficient by Julien.

# matreex 0.3.0

## News

* Disturbance is now fully implemented and happens before harvesting. An 
occurance of disturbance cancel harvesting and the size distribution of
disturbed trees is saved as harvest distribution output.

* Added BAstand in the output of the model

## Breaking changes

* `def_init_k()` require the distribution of size per hectare and don't need to multiply this distribution by `SurfEch` value anymore. See #6

* `sim_deter_forest()` output table changed the name for variable `m_i` (distribution of density by mesh i) for `n_i`. #8.

* `delay` now accept negative values, to remove a delay to any object.

* IPM build with `make_IPM()` now use the default delay set for each species (see #10). Previous simulation used `delay = 0`.

* Default parameter `correction` in `make_IPM()` and `make_mu_gr()` is now `"none"`. Previous one was `"constant"`. #11

## Documentation

* Add `{pkgdown}` usage to document the package on the www.

* Correction of vignettes

* New warning if `equil_diff` value is too low in `sim_deter_forest()`.

## Datasets

* Default delay set for each species (see #10)

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
