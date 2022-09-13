#' Constructor for IPM class
#'
#' @param IPM list of the IPM for different basal area for the species of
#' interest. List of dtCMatrix.
#' @param BA values of BA defined for the IPM, with the same names as IPM. num.
#' @param mesh Vector of all values between L and U of the IPM
#' (L being the smallest size and U the largest.). Length of this vector is
#' the number of size class in the IPM. num.
#' @param species Name of the species to run simulation on. Single char.
#' @param climatic Vector of climatic situations to run on. IPM must exist for
#' each one or else this climatic value will be skipped. int.
#' @param compress Is the IPM matrix compressed as integer (via \code{x * 1e7}).
#' Help to limit  size when saved on disc. FALSE by default. lgl.
#' @param delay Number of year delay between the recruitment of an individual
#' and it's inclusion in the IPM. This will enlarge the IPM and add sub diagonal
#' values of 1. # TODO see code{link{treeforce}{delay.ipm}}.
#'
#' @export
new_ipm <- function(IPM, BA, mesh, species, climatic,
                    compress = FALSE, delay = 0){

    IPM <- list(IPM = IPM, BA = BA, mesh = mesh,
                info = c(species = species, climatic = climatic,
                         compress = compress, delay = delay))
    class(IPM) <- "ipm"

    return(IPM)
}

#' validator for IPM class.
#'
#' @param x IPM class object
#'
#' @import checkmate
#'
#' @noRd
validate_ipm <- function(x){

    assertClass(x, "ipm")
    values <- unclass(x)
    names <- attr(x, "names")

    # check names of the object ####
    assertCharacter(names)
    if(any(names != c("IPM", "BA", "mesh", "info"))){
        stop("IPM class must be composed of elements IPM, BA, mesh and info")
    }

    # check the IPM part ####
    assertList(values$IPM, types = c("dtCMatrix", "dgCMatrix"),
               any.missing = FALSE, min.len = 1)
    # check other values ####
    assertNumeric(values$BA, lower = 0, any.missing = FALSE, min.len = 1)
    assertNumeric(values$mesh, lower = 0, any.missing = FALSE,
                  len = dim(values$IPM[[1]])[1])
    # check infos ####
    assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("species", "climatic", "compress", "delay"))){
        stop(paste0("IPM class must have info of elements species,",
                    " climatic, compress and delay"))
    }
    x
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
#' @param delay Number of year delay between the recruitment of an individual
#' and it's inclusion in the IPM. This will enlarge the IPM and add sub diagonal
#' values of 1. # TODO see code{link{treeforce}{delay.ipm}}.
#'
#' @import checkmate
#' @import here
#'
#' @noRd
old_ipm2ipm <- function(species, climatic = 1, delay = 0, path = here(),
                        replicat = 42){

    assertCharacter(species, len = 1)
    assertCharacter(path, len = 1)
    assertCount(climatic)
    assertCount(delay)
    assertCount(replicat)

    fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
    IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
    assertNumber(replicat, lower = 1, upper = length(IPM))
    IPM <- IPM[[replicat]]

    res <- validate_ipm(new_ipm(
        IPM = IPM$LIPM, BA = 1:length(IPM$LIPM), mesh = IPM$meshpts,
        species = species, climatic = climatic, delay = 0, compress = TRUE
    ))

    if(delay > 0){
        res <- delay(res, delay)
    }

    return(res)
}


