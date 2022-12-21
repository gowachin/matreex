#' Constructor for forest class
#'
#' Only used in the matreex package
#'
#' @param species List of species created with matreex package.
#' @param harv_rules Vector for harvest rules at the scale of the forest.
#' \describe{
#'   \item{Pmax}{maximum proportion of BAcut / BA}
#'   \item{dBAmin}{the minimum BA to perform cut}
#'   \item{freq}{Frequence at which the harvest will be executed.}
#' }
#' @importFrom purrr map_chr
#'
#' @keywords internal
#' @export
new_forest <- function(species = list(),
                       harv_rules = c(Pmax = 0.25, dBAmin = 3,
                                      freq = 1, alpha = 1)
                       ){

    sp <- names(species) <- map_chr(species, sp_name)
    forest <- list(
        species = species, harv_rules = harv_rules,
        info = list(species = sp,
                 clim_lab = map_chr(species, climatic))
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

    map(values$species, validate_species)
    # TODO check forest harv rules

    clim_lab <- values$info$clim_lab
    if(length(unique(clim_lab)) > 1){
        clim_ipm <- clim_lab[clim_lab != "mu_gr"]
        if(length(clim_ipm) > 1){ # D & F
            stop(paste0("Some ipm species are not defined with the same climatic name.",
                        "Check it with : map_chr(species, climatic)"))
        }
    }


    invisible(x)
}

#' Create a new forest for simulation
#'
#' A forest is a group of one of multiple species to silumate along time using
#' the IPM defined for each species and harvest rules.
#'
#' @inheritParams new_forest
#'
#' @export
forest <- function(species = list(),
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1)
                   ){

    res <- validate_forest(new_forest(
        species = species,
        harv_rules = harv_rules
    ))

    return(res)
}

#' Dev function to read old ipm
#'
#' Read old ipm and put them in the new class format.
#'
#' @param sp_name Name of the species to add in a single forest object. char.
#' @param climatic Vector of climatic situations to run on. IPM must exist for
#' each one or else this climatic value will be skipped. int.
#' @param path Place to save the resulting file. Single Char.
#' @param replicat Numeric for the simulation to select. By default, the 42th.
#'
#' @import checkmate
#'
#' @export
old_ipm2forest <- function(sp_name, climatic = 1, path = here(), replicat = 42){

    assertCharacter(sp_name)

    res <- forest(
        species = lapply(sp_name, old_ipm2species, climatic, path, replicat)
    )

    return(res)
}

