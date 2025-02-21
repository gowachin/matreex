---
title: 'matreex: Simulation IPM for forest dynamic.'
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

Integrated projection models (IPMs) are important tools to study population dynamics toward equilibrium, by predicting continuous trait distributions along time (@ellner2016). By contrast with individual based model, the continuous distribution remove demographic stochasticity and lead to identical simulations given the same input. This point is crucial in studying difference between vital rates, climatic conditions, community composition or perturbations in ecology.

# Statement of need

`matreex` is an R package speciallised on forestry dynamic. Other packages exist to use IPMs for various types of organisms (`pack`  for @Metcalf2013, `ipmr` for @Levin2021), however they focus on analysis on the matrix itselft whereas `matreex` package use the integrated matrices to simulate population dynamics.

A specifity of `matreex` package is the focus on trees. As theses species have a small annual growth rate, most of the integration effort is grouped near the diagonal of matrices. We resolved this issue using different integrations methods (Gauss-Legendre and Mid-bin) at different distances of the diagonal, which help speed up the integration. This speed is crucial as we integrated a matrix per competition level (based on species basal area) for density dependance.

<!--
figure sur l'intégration ?
-->

A crucial development effort was also given in order to simplify usage for researcher in ecology, who may work on multi-specific models with climate variation. This is made easier by adding fitted vital model for european tree species directly in the package, altough other models can still be used. The object-oriented method also makes it possible to limit the complexity of the code so that we can concentrate on a larger number of simulations to explore different starting conditions.


 <!-- but as forest are commonly managed, it was important to add different management algoritm as well as easily exploitable outputs for foresters. -->

<!--
## Climatic variabily

While the first IPM model was designed for a static climatic point, we exploited computation approximation to compute integrated matrices during simulations.
-->


## Usage and availability

`matreex` was designed to be expanded to fit new research ideas and is still in development, after being used in differents scientific publications (@kunstler2021, @guyennon2023, **TODO Barrere, Baranger ?**). The ability to simulate forest easily will help learning about management methods, climatic perturbation and species cohabitation in future european forest, by producing clear outputs that ecologist explore in R language.

`matreex` is an open source package made available under the MIT license. Installation and usage instructions can be found at the  website [TODO mettre le site en ligne avec une vrai url](https://forgemia.inra.fr/lessem/matreex)

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

**TODO : Projets européens ? Arnaud ?**

# References