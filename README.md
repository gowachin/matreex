# {matreex} <img src="man/figures/logo.png" align="right" height="138" />
<!-- badges: start -->
[![R build status](https://github.com/gowachin/matreex/workflows/R-CMD-check/badge.svg)](https://github.com/gowachin/matreex/actions?workflow=R-CMD-check)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://img.shields.io/badge/devel%20version-0.4.0-blue.svg)](https://github.com/gowachin/matreex)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/gowachin/matreex/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/gowachin/matreex)
<!-- badges: end -->

The goal of this package is to run integrated projection models of tree species in single or multi-specific density dependence context. The simulations return the size distribution dynamics along time. These models can be completed with different harvest and disturbance models and be runned untile equilibrium.

Main methods have been developped for [Kunstler *et al* (2020)](https://doi.org/10.1111/1365-2745.13533) and [Guyennon *et al* 2023](https://onlinelibrary.wiley.com/doi/10.1111/geb.13640) as well as european treee species growth/survival/recruitment models.

## Installation

### Dependencies

This package relies on very few packages listed below, that you can install with the following code.

```
deps <- c('checkmate', 'cli', 'dplyr', 'here', 'Matrix', 
          'purrr', 'rlang', 'statmod', 'stats', 'tidyr')
for (i in deps ){
  if(!require(i,character.only = TRUE))
    install.packages(i)
}
```

### Stable version

<!-- 
Be aware that anyone who installs directly from GitHub will need to explicitly request vignettes, e.g. with devtools::install_github(dependencies = TRUE, build_vignettes = TRUE).
-->

You can install the `{matreex}` package from [github](https://github.com/gowachin/matreex) with :

```
# install.packages("remotes")
remotes::install_github("gowachin/matreex")
```

### Development version

You can install the development version of {matreex} from github with :

```
# install.packages("remotes")
remotes::install_github("gowachin/matreex", ref = "dev")
```


**If you install the development version of this package, please check the NEWS.md or [changelog](https://gowachin.github.io/matreex/news/index.html) regularly. The package is in active development and some changes may break your code !**

## Usage Guide

[Getting started vignette](https://gowachin.github.io/matreex/articles/matreex.html) is now available to run basic simulations with `{matreex}` package.

A second [vignette about harvesting](https://gowachin.github.io/matreex/articles/Harvesting.html) module is also available. 

A third [vignette about disturbance](https://gowachin.github.io/matreex/articles/Disturbance.html) module is also available. 

## Support

Issues are centralized on [the github project.](https://github.com/gowachin/matreex/issues). 

## Roadmap

Future dev is listed on the github repository project to keep tracks of what we are working on.
Link is [here](https://github.com/gowachin/matreex/projects/1)

## License

Project is under MIT Licence.

<!--
## Suggestions for a good README
Every project is different, so consider which of these sections apply to yours. The sections used in the template are suggestions for most open source projects. Also keep in mind that while a README can be too long and detailed, too long is better than too short. If you think your README is too long, consider utilizing another form of documentation rather than cutting out information.

## Name
Choose a self-explaining name for your project.

## Description
Let people know what your project can do specifically. Provide context and add a link to any reference visitors might be unfamiliar with. A list of Features or a Background subsection can also be added here. If there are alternatives to your project, this is a good place to list differentiating factors.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
-->

