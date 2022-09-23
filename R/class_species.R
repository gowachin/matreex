#' Constructor for species class
#'
#' Only used in the treeforce package
#'
#' @param IPM ipm class object from the treeforce package.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most useful when creating random
#' population.
#' @param harvest_fun Function to impact the population with harvest rule.
#' Argument must be \code{x}, \code{species},  \code{harv_rule},
#' \code{BAtarget}, \code{ct} and \code{t}.
#' Should return a population state as it's take it in input, with less
#' population than before. Unless you want zombie trees. It represent the
#' distribution of the population to harvest
#' @param harv_lim Limits of harvest for a population size distribution.
#' \describe{
#'   \item{dth}{minimum diameter at which we cut the size distribution}
#'   \item{dha}{harvest diameter}
#'   \item{hmax}{maximum harvest rate for a size class}
#' }
#' @param recruit_fun Function to recruit new individual into the population.
#' Argument must be \code{BATOTSP}. # TODO rename argument !
#' Should return a single numeric value.
#'
#' @keywords internal
#' @export
new_species <- function(IPM, init_pop,
                        harvest_fun,
                        harv_lim = c(dth = 175, dha = 575, hmax = 1),
                        recruit_fun){

    species <- list(
        IPM = IPM, init_pop = init_pop,
        harvest_fun = harvest_fun, harv_lim = harv_lim,
        recruit_fun = recruit_fun,
        info = c(species = sp_name(IPM), climatic = climatic(IPM))
    )

    class(species) <- "species"

    return(species)
}

#' validator for IPM class.
#'
#' @param x IPM class object
#'
#' @import checkmate
#'
#' @noRd
validate_species <- function(x){

    assertClass(x, "species")
    values <- unclass(x)
    names <- attr(x, "names")

    # check names of the object ####
    assertCharacter(names)
    if(any(names != c("IPM", "init_pop", "harvest_fun", "harv_lim",
                      "recruit_fun", "info"))){
        stop(paste0("IPM class must be composed of elements IPM, init_pop,",
                    " harvest_fun, recruit_fun and info"))
    }

    # check all values ####
    validate_ipm(values$IPM)
    assertFunction(values$init_pop, args = c("mesh", "SurfEch"))
    assertFunction(values$harvest_fun,
                   args = c("x", "species", "targetBAcut", "ct"))
    # assertNumeric(values$harv_lim[1:2], lower = 0)
    # assertNumber(values$harv_lim[3], lower = 0, upper = 1)
    # TODO : check that X return >= 0 values of same length
    assertFunction(values$recruit_fun, args = c("BATOTSP", "BATOTNonSP",
                                                "mesh", "SurfEch"))
    # check infos ####
    assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("species", "climatic"))){
        stop("IPM class must have info of elements species and climatic")
    }

    x
}

#' Create a new species for simulation
#'
#' Only used in the treeforce package
#'
#' @inheritParams new_species
#'
#' @export
species <- function(IPM, init_pop, harvest_fun,
                    harv_lim = c(dth = 175, dha = 575, hmax = 1),
                    recruit_fun){

    res <- validate_species(new_species(
        IPM = IPM, init_pop = init_pop, harvest_fun = harvest_fun,
        recruit_fun = recruit_fun
    ))

    return(res)
}

#' Dev function to read old ipm
#'
#' Read old ipm and put them in the new class format.
#'
#' @param species Name of the species to run simulation on. Single char.
#' @param climatic Vector of climatic situations to run on. IPM must exist for
#' each one or else this climatic value will be skipped. int.
#' @param path Place to save the resulting file. Single Char.
#' @param replicat Numeric for the simulation to select. By default, the 42th.
#' @param harvest Function to impact the population with harvest rule.
#' Argument must be \code{pop}.
#' Should return a population state as it's take it in input, with less
#' population than before. Unless you want zombie trees.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most usefull when creating random
#' population.
#' @param delay Number of year delay between the recruitment of an individual
#' and it's inclusion in the IPM. This will enlarge the IPM and add sub diagonal
#' values of 1. # TODO see code{link{treeforce}{delay.ipm}}.
#'
#' @import checkmate
#' @import here
#'
#' @export
old_ipm2species <- function(species, climatic = 1, path = here(), replicat = 42,
                            harvest = def_harv, init_pop = def_init, delay = 0){

    assertCharacter(species, len = 1)
    assertCharacter(path, len = 1)
    assertCount(climatic)
    assertCount(replicat)
    assertCount(delay)

    fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
    raw_IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
    assertNumber(replicat, lower = 1, upper = length(raw_IPM))
    raw_IPM <- raw_IPM[[replicat]]

    res_ipm <- new_ipm(
        IPM = raw_IPM$LIPM, BA = 1:length(raw_IPM$LIPM), mesh = raw_IPM$meshpts,
        species = species, climatic = climatic, compress = TRUE, delay = 0
    )

    if(delay > 0){
        res_ipm <- delay(res_ipm, delay)
    }

    res <- species(
        IPM = res_ipm, init_pop = init_pop, harvest_fun = harvest,
        recruit_fun = exp_recFun(params = raw_IPM$rec$params_m,
                                 list_covs = raw_IPM$list_m)
    )

    return(res)
}


#' Default population initiation
#'
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param SurfEch Value of plot size surface in \eqn{m^2}
#'
#' @export
def_init <- function(mesh, SurfEch = 0.03) {
    ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
    ini <- exp(runif(1, -.005, .005) * mesh)
    alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    while(all(alea)){ # because god knows it's fucking possible.
                      # and it will return NaN
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    }
    ini[alea] <- 0
    res <- as.numeric(ini / sum(ct * ini) )
    res <- res + 1e-4 # HACK to limit falling in floating point trap !

    return(res)
}

#' Default population harvest
#'
#' @param x population state at time t
#' @param species ignored
#' @param targetBAcut ignored
#' @param ct define this # TODO
#'
#' @return
#' Distribution of population to harvest.
#' Values are between 0 (null harvest) and Xi.
#'
#' @export
def_harv <- function(x, species, targetBAcut, ct){
    rate <- 0.006 * (ct > 0)
    return(x * rate)
}

