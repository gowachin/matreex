#' Constructor for forest class
#'
#' Only used in the treeforce package
#'
#' @param species List of species created with treeforce package.
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
                                      freq = 10)){

    sp <- names(species) <- map_chr(species, sp_name)
    forest <- list(
        species = species, harv_rules = harv_rules,
        info = list(species = sp,
                 climatic = map_dbl(species, climatic))
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
    if(length(unique(values$info$climatic)) > 1){
        stop("All species are not defined for the same climatic.")
    }

    x
}

#' Create a new forest for simulation
#'
#' Only used in the treeforce package
#'
#' @inheritParams new_forest
#'
#' @export
forest <- function(species = list(),
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 10)
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

