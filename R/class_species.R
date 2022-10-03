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
        info = c(species = sp_name(IPM), clim_lab = climatic(IPM))
    )

    class(species) <- "species"

    return(species)
}

#' validator for species class.
#'
#' @param x species class object
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
        stop(paste0("species class must be composed of elements IPM, init_pop,",
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
    if(any(names(values$info) != c("species", "clim_lab"))){
        stop("species class must have info of elements species and clim_lab")
    }

    x
}

#' Create a new species for simulation
#'
#' Species are defined by an IPM which is a transition matrix from size between
#' t and t+1, recruitment and harvest functions. Each species has these items
#' defined for a given climate.
#' An additionnal vector of harvest parameers is required with minimal size to
#' harvest (dth), size above wich harvest is constant (dha).
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
#' Read old ipm and put them in the new species class format.
#'
#' @param species Name of the species to load. Single char.
#' @param climatic Numeric that coded old climatic state. int.
#' @param path Place to load the previous IPM file. Single Char.
#' @param replicat Numeric for the model to select. By default, the 42th.
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
    assertCount(delay)

    fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
    raw_IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
    assertNumber(replicat, lower = 1, upper = length(raw_IPM))
    raw_IPM <- raw_IPM[[replicat]]

    res_ipm <- new_ipm(
        IPM = raw_IPM$LIPM, BA = 1:length(raw_IPM$LIPM), mesh = raw_IPM$meshpts,
        species = species, climatic = drop(as.matrix(raw_IPM$list_m)),
        clim_lab = climatic, compress = TRUE, delay = 0
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
#' The population will initiate with a random distribution to match a basal area
#' of 1.
#'
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param SurfEch Value of plot size surface in \eqn{m^2}
#'
#' @importFrom stats runif rbinom
#'
#' @export
def_init <- function(mesh, SurfEch = 0.03) {
    ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
    ini <- exp(runif(1, -.005, .005) * mesh)
    alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    while(all(alea)){ # because god knows it's fucking possible that alea is
                      # all FALSE and it will return NaN
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    }
    ini[alea] <- 0
    res <- as.numeric(ini / sum(ct * ini) )
    res <- res + 1e-4 # HACK to limit falling in floating point trap !

    return(res)
}


#' Init population at BA
#'
#' This function modify the def_init function to start at a given BA with
#' the same process of random distribution.
#'
#' @param BA Basal area targeted. This single value must be above 0 but can be
#' very close (minimal accepted value is 1e-10)
#'
#' @return
#' Function similar to def_init
#'
#' @import checkmate
#'
#' @export
def_initBA <- function(BA = 1){

    assertNumber(BA, lower = 1e-10)
    fun <- def_init
    force(BA)
    res <- NULL # HACK because res must be binded
    body(fun)[[8]] <- call2("<-", expr(res), call2("*", expr(res), BA))

    return(fun)
}


#' Init population with precise distribution
#'
#' @param x Distribution to draw systematically. This distribution should
#' be composed of values in \code{[0, Inf]} values with a sum superior to 0.
#'
#' @details
#' The resulting function will check if the provided vector is the same length
#' as mesh.
#'
#' @return
#' Function similar def_init but with no random effect anymore.
#'
#' @import checkmate
#' @export
def_init_k <- function(x){

    assertNumeric(x, lower = 0, any.missing = FALSE)
    assertTRUE(sum(x) > 0)

    force(x)
    fun <- function(mesh, SurfEch = 0.03) {
        if(length(x) != length(mesh)){
            stop(paste0("A species initiate with a define distribution with ",
                        "different length that it's mesh. Check sp$init_pop ",
                        "functions using def_init_k !"))
        }
        return(x)
    }

    return(fun)
}


#' Default population harvest
#'
#' Constant rate harvest of 0.06 percent per year
#' (check if harvest frequence is 1 in forest object).
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

