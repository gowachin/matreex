---
title: 'matreex: Simulating European forest dynamics with IPM.'
tags:
  - R
  - forestry
  - dynamics
  - ecology
authors:
  - name: Maxime Jaunatre
    orcid: 0009-0002-2816-1677
    corresponding: true
    affiliation: "1"
  - name: Julien Barrere
    orcid: 0000-0002-6686-726X
    affiliation: "1,2"
  - name: Björn Reineking
    orcid: 0000-0001-5277-9181
    affiliation: "1"
  - name: Arnaud Guyennon
    orcid: 0000-0003-2178-3801
    affiliation: "3"
  - name: Thomas Cordonnier
    orcid: 0000-0003-3684-4662
    affiliation: "4"
  - name: Georges Kunstler
    orcid: 0000-0002-2544-1940
    corresponding: true
    affiliation: "1"
affiliations:
 - name: Univ. Grenoble Alpes, INRAE, LESSEM, St-Martin-d'Hères, France
   index: 1
   ror: 01a0ez112
 - name: INRAE, Aix Marseille Université, UMR RECOVER, Aix-en-Provence, France
   index: 2
 - name: Independent Researcher, France
   index: 3
 - name: ONF, Département Recherche Développement et Innovation, 21 rue du Muguet, 39100 Dole.
   index: 4
date: 11 February 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

<!--docker run --rm --volume $PWD/inst/paper:/data  --user $(id -u):$(id -g) --env JOURNAL=joss openjournals/inara -->

# Summary

Integral projection models (IPMs) are powerful tools for studying the temporal dynamics of populations structured by continuous traits, allowing for predictions of changes in trait distributions over time [@ellner2016].
Unlike individual-based or cohort-based models, which represent populations as discrete populations, IPMs describe populations as continuous populations, integrating over the uncertainty of demographic processes.
This removes demographic stochasticity and results in fully deterministic simulations, which are complementary to individual-based models (IBMs).

Here, we introduce `matreex`, an R package specifically designed to build IPMs for European forest tree species.
Our package includes pre-fitted species-specific growth, survival and recruitment functions that account for the effect of climate and competition, and functions to efficiently integrate IPMs and run temporal simulations of single-species or multispecies forest communities until equilibrium.
This package complements existing R packages for IPMs, such as `ipmr` [@Metcalf2013] and `IPMpack` [@Levin2021], which are not specifically designed for forest ecosystems.

# Statement of need

<!--
Statement of need: A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work.
-->

IPMs are rarely applied to forest ecosystems due to the complexity of tree growth kernels, which are challenging to integrate, making the construction of forest IPMs particularly difficult.
`matreex` is an R package specifically designed to build IPMs and run simulations for European forests.
These IPMs can be used by the forest ecology research community, to study population dynamic and forest communities trajectories.
In `matreex` IPM simulations, it is also possible to include temporally variable climatic conditions, natural disturbances, harvesting scenarios and regional dispersal affecting population dynaminc depending on tree species sensitivity and stand structure.

# State of the field

Other, more generalist, packages exist to build IPMs for various types of organisms (`IPMpack` [@Metcalf2013], `ipmr` [@Levin2021]).
However they focus on analysis on the matrix itself whereas `matreex` use the integrated matrices to simulate population dynamics.
In addition `matreex` integrates numerous functions to run simulations for multispecies communities under various harvesting and disturbance scenarios.

# Software design

A key feature of the `matreex` package is the development of IPM integration functions designed for complex forest tree growth kernels.
As trees have a small annual growth rate, most of the integration effort is concentrated near the diagonal of the matrix.
We achieve this by combining different integration methods (Gauss-Legendre and Mid-bin) at different distances of the diagonal, which helps speed up the integration.
This speed is crucial as we integrate a matrix per competition level (based on species basal area) to account for density dependence.
More details about the integration are provided in the [`matreex` webpage](https://lessem.pages-forge.inrae.fr/rpackages/matreex/articles/building_ipm.html)

![Figure 1: Combination of different integration methods. Dashed line is the identity where $z_t=z_t+1$ and dark blue distribution is an expected distribution of the growth kernel.\label{fig:band_matrix](fig/figures_files/figure-html/band_matrix-1.png)

A crucial development objective was to simplify the workflow for ecological researchers, who may work on multispecies models with climate variation, disturbances, and harvesting.
To facilitate this, `matreex` provides fitted vital models for European tree species directly in the package [@kunstler2021; @guyennon2023; @barrere2024], although other new vital models can be used.
The object-oriented architecture limits code complexity, allowing users to focus on designing large simulation experiments to tackle their ecological questions.

```
library(matreex)
library(dplyr)
library(ggplot2)

# select a climate to run in
data("climate_species")
# N is the climate number to user, 2 being the optimum climate for the species
climate <- subset(climate_species, N == 2 & sp == "Fagus_sylvatica", select = -sp)

# integrate ipms and store them in species object
Picea <- species(IPM = make_IPM(
    species = "Picea_abies", fit = fit_Picea_abies,
    climate = climate, clim_lab = "optimum clim", 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:60
), init_pop = def_initBA(1))
Betula <- species(IPM = make_IPM(
    species = "Betula", fit = fit_Betula,
    climate = climate, clim_lab = "optimum clim", 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Betula) * 1.1),
    BA = 0:60
), init_pop = def_initBA(1))
Fagus <- species(IPM = make_IPM(
    species = "Fagus_sylvatica", fit = fit_Fagus_sylvatica,
    climate = climate, clim_lab = "optimum clim", 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
    BA = 0:60
), init_pop = def_initBA(1))
Abies <- species(IPM = make_IPM(
    species = "Abies_alba", fit = fit_Abies_alba,
    climate = climate, clim_lab = "optimum clim", 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Abies_alba) * 1.1),
    BA = 0:60
), init_pop = def_initBA(1))

# assemble species in a forest object
Forest <- forest(species = list(Picea = Picea, Abies = Abies, 
                                   Fagus = Fagus, Betula = Betula))
set.seed(42) # The seed is here for initial population random functions.

# Run simulation and plot it
Sim <- sim_deter_forest(
    Forest, 
    tlim = 1000, 
    equil_time = 1000, equil_dist = 50, equil_diff = 1,
    SurfEch = 0.03,
    verbose = TRUE
)
Sim  %>%
    dplyr::filter(var == "BAsp", ! equil) |>
    ggplot(aes(x = time, y = value, color = species)) +
    geom_line(linewidth = .4) + 
    ylab("Basal Area (m2)") + xlab("Simulation time (years)") 
```

![Figure 2: Simulation output for 4 species. \label{fig:simulation}](fig/simulation.png)

## Climatic temporal variability

Modelling forest dynamics under fluctuating climatic conditions can be computationally expensive because the IPMs growth kernel must be integrated for every climatic condition.
To avoid this high computation cost, we implement a new method involving the pre-integration of IPM growth matrix blocks for different mean growth rates.
The IPM for each climatic condition is then reassembled from these mean growth rate IPM matrix blocks (mu matrix).
This allows to speed-up the simulations.
This is described in the [matreex climate variation vignette](https://lessem.pages-forge.inrae.fr/rpackages/matreex/articles/mu_simulation.html).

## Disturbance

One key originality of `matreex` is the possibility to apply storm, fire, biotic of snow disturbances of varying intensity (ranging from 0 to 1) at any time of the IPM simulations.
When a disturbance strikes the stand a given year of the simulation, the survival function is replaced by the species-specific equations from @barrere2024.
These equations quantify the annual mortality probability of a tree in a disturbed stand as a function of its species, diameter at breast height, stand structure, nature and intensity of the disturbance.
A disturbance striking two different stands with the same intensity will thus result in different mortality rates, depending notably on the sensitivity of the tree species present in the stand.
Disturbances in `matreex` are described in details in @barrere2024 and in @barrereprep.
[Figure 3](@fig:disturbance) shows an example of multispecies simulations with storm disturbance.

![Figure 3: Simulation output for 3 species with a disturbance at $time = 2600$.](fig/disturbance.png){#fig:disturbance}

## Regional dispersal

Most stand-scale forest dynamics models simulate closed systems, where only the tree species already present in the stand contribute to the recruitment of new trees.
This limitation perclude the simulation of immigration from external species, which is a key process of forest dynamics, particularily under climate change.
To overcome this limitation, we included the possibility to split the recruitment function in two component : (i) within-plot dispersal that depends on the sum of basal area of fecund tree species in the plot, and (ii) external dispersal, that depends on the regional pool.
This regional pool approach is extensively presented in @barrereprep.

## Harvesting

Since most temperate forests are managed, it is crucial to incorporate silvicultural interventions into the simulations.
We implemented three management strategies:

- First, we implemented a simple constant annual harvesting rates accounting for the effect of the harvesting rates observed in NFI data used for the model calibration [@kunstler2021].

- Second, we implement an even-aged management.
The objective is to apply typical even-aged harvesting, based on a single cohort.
Trees are harvested with successive thinning during stand development until the final harvest.
Thinning harvest are based on the distance to a self-thinning boundary, based on @Aussenac2021.
This is easily connected with management guidelines.

- Third, we implement an unven-aged harvesting.
The uneven-aged harvest scenario consists in selective harvesting across all size classes with the objective to reach a stable size structure with continuous replacement of large mature trees.
This scenario depends on the basal area of the stand and the size distribution of the trees (building on @Lafond2014).
These three managements are described in the [matreex harvesting vignette](https://lessem.pages-forge.inrae.fr/rpackages/matreex/articles/Harvesting.html).

# Research impact statement

`matreex` was designed to be easily adapted to various research ideas and is in continuous development.
It has already been used in several scientific publications tackling diverse scientific questions [@kunstler2021; @guyennon2023; @barrere2024; @barrereprep; @barangerprep].
The ability to simulate forests easily with a dedicated R package will help ecologists to analyse the effect of climate change, shifting disturbance regimes, and their interplay with forest management across European forests.

`matreex` is an open-source package made available under the MIT license.
Installation and usage instructions can be found at the [website](https://lessem.pages-forge.inrae.fr/rpackages/matreex/)

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.

# Acknowledgements

JB, MJ, BR and GK are funded through the BiodivClim ERA-Net Cofund (joint BiodivERsA Call on “Biodiversity and Climate Change”, 2019-2020) with national co–funding through ANR (France, project ANR-20-EBI5-0005-03).
MJ and GK were funded by the ANR DECLIC (grant ANR-1520-CE32-0005-01) and REGE-ADAPT PEPR FORESTT France 2030 (ANR-24-PEFO-0006).
JB, MJ, BR and GK are funded by the RESONATE H2020 project (grant 101000574).
GK, BR and AG received support from the REFORCE – EU FP7ERA-NET Sumforest 2016 through the call ‘Sustainable forests for the society of the future’, with the ANR as national funding agency (grant ANR-16-SUMF-0002).

# References
