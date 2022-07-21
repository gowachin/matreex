#' Constructor for forest class
#'
#' Only used in the treeforce package
#'
#' @param IPM ipm class object from the treeforce package.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most useful when creating random
#' population.
#' @param harvest_fun Function to impact the population with harvest rule.
#' Argument must be \code{pop}.
#' Should return a population state as it's take it in input, with less
#' population than before. Unless you want zombie trees.
#' @param recruit_fun Function to recruit new individual into the population.
#' Argument must be \code{BATOTSP}. # TODO rename argument !
#' Should return a single numeric value.
#'
#' @noRD
new_forest <- function(species = list(), ){

    forest <- list(
        IPM = IPM, init_pop = init_pop,
        harvest_fun = harvest_fun, recruit_fun = recruit_fun,
        info = c(forest = forest(IPM), climatic = climatic(IPM))
    )

    class(forest) <- "forest"

    return(forest)
}

#' validator for IPM class.
#'
#' @param x IPM class object
#'
#' @import checkmate
#'
#' @noRd
validate_forest <- function(x){

    values <- unclass(x)
    names <- attr(x, "names")

    # check names of the object ####
    assertCharacter(names)
    if(any(names != c("IPM", "init_pop", "harvest_fun",
                      "recruit_fun", "info"))){
        stop(paste0("IPM class must be composed of elements IPM, init_pop",
                    " harvest_fun, recruit_fun and info"))
    }

    # check all values ####
    validate_forest(values$IPM)
    assertFunction(value$init_pop, args = c("mesh", "SurfEch}"))
    assertFunction(values$harvest_fun, args = c("x"))
    assertFunction(values$recruit_fun, args = c("BATOTSP")) # FIXME args must be changed.
    # check infos ####
    assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("forest", "climatic"))){
        stop("IPM class must have info of elements forest and climatic")
    }

    x
}

#' Create a new forest for simulation
#'
#' Only used in the treeforce package
#'
#' @param IPM ipm class object from the treeforce package.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most useful when creating random
#' population.
#' @param harvest_fun Function to impact the population with harvest rule.
#' Argument must be \code{pop}.
#' Should return a population state as it's take it in input, with less
#' population than before. Unless you want zombie trees.
#' @param recruit_fun Function to recruit new individual into the population.
#' Argument must be \code{BATOTSP}. # TODO rename argument !
#' Should return a single numeric value.
#'
#' @export
forest <- function(IPM, init_pop, harvest_fun, recruit_fun){

    res <- validate_forest(new_forest(
        IPM = IPM, init_pop = init_pop, harvest_fun = harvest_fun,
        recruit_fun = recruit_fun
    ))

    return(res)
}

#' Dev function to read old ipm
#'
#' Read old ipm and put them in the new class format.
#'
#' @param forest Name of the forest to run simulation on. Single char.
#' @param climatic Vector of climatic situations to run on. IPM must exist for
#' each one or else this climatic value will be skipped. int.
#' @param path Place to save the resulting file. Single Char.
#' @param replicat Numeric for the simulation to select. By default, the 42th.
#' @param harvest_fun Function to impact the population with harvest rule.
#' Argument must be \code{pop}.
#' Should return a population state as it's take it in input, with less
#' population than before. Unless you want zombie trees.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most usefull when creating random
#' population.
#'
#' @import checkmate
#'
#' @export
old_ipm2forest <- function(forest, climatic = 1, path = here(), replicat = 42,
                            harvest = def_harv, init_pop = def_init){

    assertCharacter(forest, len = 1)
    assertCharacter(path, len = 1)
    assertCount(climatic)
    assertCount(replicat)

    fIPM <- here(path, "output", forest, paste0("IPM_Clim_", climatic, ".Rds"))
    IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
    assertNumber(replicat, lower = 1, upper = length(IPM))
    IPM <- IPM[[42]]

    res_ipm <- new_ipm(
        IPM = IPM$LIPM, BA = 1:length(IPM$LIPM), mesh = IPM$meshpts,
        forest = forest, climatic = climatic, compress = TRUE
    )

    res <- forest(
        IPM = res_ipm, init_pop = init_pop, harvest_fun = harvest,
        recruit_fun = IPM$RecFun
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

    return(res)
}

#' Default population harvest
#'
#' @param x population state at time t
#'
#' @export
def_harv <- function(x){
    return(x * (1 - 0.006))
}

