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
    equal-contrib: true
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Julien Barrere
    orcid: 0000-0002-6686-726X
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Bjoern Reiniking
    orcid: 0000-0001-5277-9181
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Laura Touzot
    orcid: 0000-0003-0445-554X
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Arnaud Guyennon
    orcid: 0000-0003-2178-3801
    affiliation: "2"
  - name: Georges Kunstler
    orcid: 0000-0002-2544-1940
    equal-contrib: true
    corresponding: true
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Univ. Grenoble Alpes, INRAE, LESSEM, St-Martin-d'Hères, France
   index: 1
   ror: 00hx57361
 - name: Independent Researcher, France
   index: 2
date: 11 February 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Introduction

Integrated projection models (IPMs) are powerful tools for studying the temporal dynamics of populations structured by continuous traits, allowing predictions of changes in trait distributions over time (@ellner2016). Unlike individual-based or cohort-based models, which represent populations as finite (discrete) populations, IPMs describe populations as infinite (continuous) populations, integrating over the uncertainty of demographic processes. This removes demographic stochasticity and results in fully deterministic simulations which is complementary to IBM models. IPMs are rarely applied to forest ecosystems due to the complexity of tree growth kernels, which are challenging to integrate, making the construction of forest IPMs particularly difficult.

Here, we introduce an R package specifically designed to build IPMs for European forest tree species. Our package includes fitted species-specific functions of growth, survival and recruitment accounting for the effect of climate and competition, and functions to efficiently integrate IPMs and run temporal simulations of single-species or multispecies forest communities until equilibrium. We also included the possibility to simulate natural disturbances (storm, fire, biotic and snow) affecting population survival depending on tree species sensitivity and stand structure. This package complements existing R packages for IPMs, such as ipmr and IPMpack, which are not specifically designed for forest ecosystems.

# Statement of need

`matreex` is an R package specifically designed to build IPMs and run simulations for European forests. Other, more generalist, packages exist to build IPMs for various types of organisms (`IPMpack`  for @Metcalf2013, `ipmr` for @Levin2021). In addition `matreex` integrates numerous functions to run simulations for multispecies communities with harvesting and disturbance scenarios.

A specificity of `matreex` package is the development of IPM integration functions focused on trees. As trees have a small annual growth rate, most of the integration effort is put near the diagonal of the matrice. We achieve this by combining different integration methods (Gauss-Legendre and Mid-bin) at different distances of the diagonal, which helps speed up the integration (REF TO FIGUE). This speed is crucial as we integrated a matrix per competition level (based on species basal area) to account for density dependence. More details about the integration are provided in the `matreex` webpage REF.

<!--
figure sur l'intégration ? GK= Juste la 1er fig pas celle donnant le detail de la band matrix
-->

A crucial development effort was also to simplify usage for ecological researchers, who may work on multi-specific models with climate variation, disturbances, and harvesting. This is made easier by providing fitted vital models for European tree species directly in the package RAJOUTER REF PAPIER, although other new models can still be used. The object-oriented method limits the complexity of the code so that users can concentrate on setting up large simulation experiments to tackle their ecological questions. The figure XX REF shows an example of multispecies simulations with storm disturbance from Barrere et al. 2024 ADD REF. ADD FIGURE

 <!-- but as forest are commonly managed, it was important to add different management algoritm as well as easily exploitable outputs for foresters. -->

<!--
## Climatic temporal variability

While the first IPM model was designed for a static climate, we developed new methods to approximate the IPM in a temporally variable climate in an efficient way in the simulations.
-->

## Disturbance

One key originality of matreex is the possibility to apply storm, fire, biotic of snow disturbances of different intensity (ranging from 0 to 1) at any time of the  IPM simulations. When a disturbance strikes a given year of the simulation, the survival function is replaced by the species-specific equations from Barrere et al. (2023), which quantify the annual mortality probability of a tree in a disturbed stand as a function of its species, diameter at breast height, stand structure, and the nature and intensity of the disturbance. A disturbance striking two different stands with the same intensity will thus result in different mortality rates, depending notably on the sensitivity of the tree species present in each plot to that specific disturbance.

## Regional dispersal

Most stand-scale forest dynamics models tend to simulate closed systems, where only the tree species already present in the stand contribute to the recruitment of new trees. This limitation prevents the possibility to simulate immigration from external species, which is a key process of forest dynamics, particularily in a context of climate change. To overcome this limitation, we included the possibility to split the recruitment function in two component : (i) within-plot dispersal that depends on the summed basal area of fecund tree species in the plot, and (ii) external dispersal, that depends on a regional pool associated with the plot simulated. This regional pool approach is extensively presented in Barrere et al. (in revision).

<!--
## Harvesting
TODO
-->

## Usage and availability

`matreex` was designed to be expanded to fit new research ideas and is in continuous development and have been used in different scientific publications (@kunstler2021, @guyennon2023, @barrere2024; Barrere et al. in prep, Baranger et al in prep**). The ability to simulate forests easily with a designed R package will help ecologists analyse the effect of climate change, change in disturbance regimes, and the interplay with forest management of European forests.

`matreex` is an open-source package made available under the MIT license. Installation and usage instructions can be found at the  website [TODO mettre le site en ligne avec une vrai url](https://forgemia.inra.fr/lessem/matreex)

<!--
# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

-->

# Acknowledgements

**TODO : Projets européens ? Arnaud ? RESONATE FUNPOTENTIAL**

JB, MJ, BR and GK are funded through the BiodivClim ERA-Net Cofund,(joint BiodivERsA Call on “Biodiversity and Climate Change”, 2019-2020) with national co–funding through ANR (France, project ANR-20-EBI5-0005-03).
GK and LT were funded by the ANR DECLIC (grant ANR-1520-CE32-0005-01) and REGE-ADAPT PEPR FORESTT France 2030 (ANR-24-PEFO-0006).
JB, MJ, BR and GK are funded by the RESONATE H2020 project (grant 101000574).
G.K. and A.G. received support from the REFORCE – EU FP7ERA-NET Sumforest 2016 through the call ‘Sustainable forests for the society of the future’, with the ANR as national funding agency (grant ANR-16-SUMF-0002).

# References
