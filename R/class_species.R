#' Constructor for species class
#'
#' Only used in the matreex package
#'
#' @param IPM ipm class object from the matreex package.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most useful when creating random
#' population.
#' @param harvest_fun Function to impact the population with harvest rule.
#' Argument must be \code{x}, \code{species},  \code{...},
#' @param disturb_fun Function to impact the population with possibles
#' disturbances. Extra care is needed to give this function all needed
#' parameters. Default is \code{def_disturb}.
#' @param disturb_coef Species coefficient  for disturbance reaction. These
#' values and names are highly dependent on the disturbance function.
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
#' @param rdi_coef Coefficient for RDI curve used in even harvest.
#' The model require the intercept and slope.
#' @param type Type of the tree, choosing between "Broadleaf" and "Coniferous".
#' This value is only used during biotic disturbance with a specific disturb_fun.
#' This is experimental
#'
#' @keywords internal
#' @export
new_species <- function(IPM, init_pop,
                        harvest_fun, disturb_fun,
                        harv_lim = c(dth = 175, dha = 575, hmax = 1),
                        rdi_coef = NULL,
                        disturb_coef = NULL,
                        type = c("Undefined", "Broadleaf", "Coniferous")
                        ){

    type <- match.arg(type)

    if(inherits(IPM, "ipm")){
        rec <- exp_recFun(params = IPM$fit$rec$params_m,
                          list_covs = IPM$climatic)
    } else if(inherits(IPM, "mu_gr")){
        rec <- "to define"
    } else {
        stop("IPM must either be an ipm or mu_gr object.")
    }

    species <- list(
        IPM = IPM, init_pop = init_pop,
        harvest_fun = harvest_fun, harv_lim = harv_lim,
        disturb_fun = disturb_fun,
        rdi_coef = rdi_coef, disturb_coef = disturb_coef,
        recruit_fun = rec,
        info = c(species = sp_name(IPM), clim_lab = climatic(IPM), type = type)
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
                      "disturb_fun", "rdi_coef", "disturb_coef",
                      "recruit_fun", "info"))){
        stop(paste0("species class must be composed of elements IPM, init_pop,",
                    " harvest_fun, harv_lim, disturb_fun, rdi_coef, disturb_coef, ",
                    "recruit_fun and info"))
    }

    # check all values ####
    if(inherits(values$IPM, "ipm")){
        validate_ipm(values$IPM)
    } else if(inherits(values$IPM, "mu_gr")){
        validate_mu_gr(values$IPM)
    } else {
        stop("IPM must either be an ipm or mu_gr object.")
    }

    assertFunction(values$init_pop, args = c("mesh", "SurfEch"))
    assertFunction(values$harvest_fun,
                   args = c("x", "species", "..."))
    assertFunction(values$disturb_fun,
                   args = c("x", "species", "disturb", "..."))
    # assertNumeric(values$harv_lim[1:2], lower = 0)
    # assertNumber(values$harv_lim[3], lower = 0, upper = 1)
    # TODO : check that X return >= 0 values of same length
    if(inherits(values$IPM, "ipm")){
        assertFunction(values$recruit_fun, args = c("BATOTSP", "BATOTNonSP",
                                                    "mesh", "SurfEch"))
    } else if(inherits(values$IPM, "mu_gr")){
        assertString(values$recruit_fun)
    }
    # check infos ####
    assertCharacter(values$info, any.missing = FALSE)
    if(any(names(values$info) != c("species", "clim_lab", "type"))){
        stop("species class must have info of elements species, clim_lab and type")
    }

    invisible(x)
}

#' Create a new species for simulation
#'
#' Species are defined by an IPM which is a transition matrix from size between
#' t and t+1, recruitment and harvest functions (see Details). Each species has these items
#' defined for a given climate.
#' An additionnal vector of harvest parameters is required with minimal size to
#' harvest (dth), size above wich harvest is constant (dha).
#'
#' @inheritParams new_species
#'
#' @family functions for initiating species population during simulation
#' @family functions that defines harvest rules for a species.
#'
#' @details
#' A species is defined by an IPM that is an integrated prediction matrix
#' for growth and survival functions of the species. Since the species has other
#' functions defined, they are accessible and editable.
#'  \describe{
#'   \item{\code{init_pop}}{Function to initiate a new population. Default is
#'   \code{\link[matreex]{def_init}}.
#'   }
#'   \item{\code{recruit_fun}}{Function that give a distribution for recruits.
#'   The default is defined from models associated with the IPM
#'   (\code{x$IPM$rec}) but it's possible to replace it. For example you can
#'   nullify the recruitment to simulate extinction.
#'   }
#'   \item{\code{harvest_fun}}{Function that give harvest density distribution
#'   when an harvest event occurs (this frequence is set at the forest scale.).
#'   The default function is \code{\link[matreex]{def_harv}} with a constant
#'   harvest rate of 0.6 percent. Other functions are
#'   \code{\link[matreex]{Uneven_harv}} and \code{\link[matreex]{Even_harv}}.
#'   }
#' }
#'
#' @aliases harvest_fun init_pop recruit_fun
#'
#' @export
species <- function(IPM, init_pop = def_init, harvest_fun = def_harv,
                    disturb_fun = def_disturb,
                    harv_lim = c(dth = 175, dha = 575, hmax = 1),
                    rdi_coef = NULL, disturb_coef = NULL,
                    type = c("Undefined", "Broadleaf", "Coniferous")){

    type <- match.arg(type)

    res <- validate_species(new_species(
        IPM = IPM, init_pop = init_pop, harvest_fun = harvest_fun,
        disturb_fun = disturb_fun,
        harv_lim = harv_lim, rdi_coef = rdi_coef, disturb_coef = disturb_coef,
        type = type
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
#' @param disturb Function to impact the population with possibles
#' disturbances. Extra care is needed to give this function all needed parameters/
#' Default is \code{def_disturb}.
#' @param init_pop Function to initiate the population at simulation start.
#' Arguments must be \code{mesh} and \code{SurfEch}.
#' Using the arguments is not mandatory, it's most usefull when creating random
#' population.
#' @param delay Number of year delay between the recruitment of an individual
#' and it's inclusion in the IPM. This will enlarge the IPM and add sub diagonal
#' values of 1. See \code{\link[matreex]{delay}}.
#'
#' @import checkmate
#' @import here
#'
#' @export
old_ipm2species <- function(species, climatic = 1,
                            path = here(), replicat = 42,
                            harvest = def_harv,
                            disturb = def_disturb,
                            init_pop = def_init,
                            delay = 0){

    assertCharacter(species, len = 1)
    assertCharacter(path, len = 1)
    assertCount(climatic)
    assertCount(delay)

    fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
    raw_IPM <- readRDS(assertFileExists(fIPM)) # note 10" to load...
    assertNumber(replicat, lower = 1, upper = length(raw_IPM))
    raw_IPM <- raw_IPM[[replicat]]

    res_ipm <- new_ipm(
        IPM = raw_IPM$LIPM, BA = 1:length(raw_IPM$LIPM), mesh = raw_IPM$meshpts,
        species = species, correction = "constant",
        climatic = drop(as.matrix(raw_IPM$list_m)),
        fit = old_fit2fit(species, path = path, replicat = replicat, mean = FALSE),
        clim_lab = climatic, compress = TRUE, survival = TRUE, delay = 0
    )

    if(delay > 0){
        res_ipm <- delay(res_ipm, delay)
    }

    rdi <- matreex::rdi_coef
    rdi <- drop(as.matrix(rdi[rdi$species == species,c("intercept", "slope")]))

    disturb_c <- matreex::disturb_coef
    disturb_c <- disturb_c[disturb_c$species == species,]

    type <- matreex::tree_type
    type <- type[type$species == species, "type"]

    res <- species(
        IPM = res_ipm, init_pop = init_pop, harvest_fun = harvest,
        disturb_fun  = disturb, rdi_coef = rdi, disturb_coef = disturb_c,
        type = type
    )

    return(res)
}


#' Default population initialization
#'
#' The population will initiate with a random distribution to match a basal area
#' of 1.
#'
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param SurfEch Value of plot size surface in \eqn{m^2}
#'
#' @details
#' This function is the default function for \code{init_pop} function of a species
#' and always takes arguments \code{mesh} and \code{SurfEch}.
#'
#' @importFrom stats runif rbinom
#' @family functions for initiating species population during simulation
#'
#' @export
def_init <- function(mesh, SurfEch = 0.03) {
    ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
    ini <- exp(runif(1, -.005, -1e-4) * mesh)
    alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    while(all(alea)){ # because god knows it's fucking possible that alea is
                      # all FALSE and it will return NaN
        alea <- rbinom(length(mesh), 1, runif(1, .6, .9)) == 1
    }
    ini[alea] <- 0
    res <- as.numeric(ini / sum(ct * ini) )
    res <- res + 1e-4 # HACK to limit falling in floating point trap !
                      # also line to add BA later if needed
    return(res)
}


#' Default even population initialization
#'
#' The population will initiate with only individual in the 5 first cells of the
#' mesh.
#'
#' @param mesh all possible states of a population, based on an IPM.
#' Minimal and maximal values are respectively U and L, for a total number of
#' m states.
#' @param SurfEch Value of plot size surface in \eqn{m^2}
#'
#' @importFrom stats runif
#' @family functions for initiating species population during simulation
#'
#' @export
def_init_even <- function(mesh, SurfEch = 0.03) {
  ct <- drop(Buildct(mesh = mesh, SurfEch = SurfEch))
  x <- exp(runif(1, -.005, .005) * mesh)
  if(any(mesh == 0)){
      sel <- max(which(mesh == 0)) + 1 # mesh starting point once delay is applied
  } else {
      sel <- 1
  }
  x[-(sel:(sel+4))] <- 0
  res <- as.numeric(x / sum(ct * x))
  res <- res + 1e-15 # FIXME for pending point...while BA = 0 is not integrated
  res <- res # line to add BA later if needed
  return(res)
}


#' Init population at BA
#'
#' This function modify the def_init function to start at a given BA with
#' the same process of random distribution.
#'
#' @param BA Basal area targeted. This single value must be above 0 but can be
#' very close (minimal accepted value is 1e-10)
#' @param fun Function to modify, single chr in choices.
#'
#' @return
#' Function similar to def_init
#'
#' @import checkmate
#' @family functions for initiating species population during simulation
#'
#' @export
def_initBA <- function(BA = 1, fun = c("def_init", "def_init_even")){

    assertNumber(BA, lower = 1e-10)
    fun <- match.arg(fun)

    fun <- switch(fun, def_init = def_init, def_init_even = def_init_even)
    force(BA)
    res <- NULL # hack because res must be binded
    l <- length(as.list(body(fun))) - 1 # last - 1 line to edit.
    body(fun)[[l]] <- call2("<-", expr(res), call2("*", expr(res), BA))

    return(fun)
}



#' Init population with precise distribution
#'
#' @param x Distribution to draw systematically. This distribution should
#' be composed of values in \code{[0, Inf]} values with a sum superior to 0.
#' This is the distribution per hectare and not for the sampled plot.
#'
#' @details
#' The resulting function will check if the provided vector is the same length
#' as mesh.
#'
#' @return
#' Function similar def_init but with no random effect anymore.
#'
#' @import checkmate
#' @family functions for initiating species population during simulation
#'
#' @export
def_init_k <- function(x){

    assertNumeric(x, lower = 0, any.missing = FALSE)
    assertTRUE(sum(x) >= 0)
    if(sum(x) == 0){
        warning(paste0("sum(x) is equal to 0, the species will not be present",
                       " in the forest. Be sure this is intentional."))
    }
    # assertTRUE(sum(x) > 0)

    force(x)
    fun <- function(mesh, SurfEch = 0.03) {
        if(length(x) != length(mesh)){
            stop(paste0("A species initiate with a define distribution with ",
                        "different length that it's mesh. Check sp$init_pop ",
                        "functions using def_init_k !"))
        }
        return(x * SurfEch)
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
#' @param ... Variables used in this case of Uneven harvest
#' \describe{
#'  \item{ct}{is the vector to compute BA with x (ct = Buildct(mesh, SurfEch))}
#' }
#'
#' @details
#' This function is the default function for \code{harvest_fun} function of
#' a species and always takes arguments \code{x}, \code{species} plus specific
#' argument from different harvest models..
#'
#' @return
#' Distribution of population to harvest.
#' Values are between 0 (null harvest) and Xi.
#'
#' @family functions that defines harvest rules for a species.
#'
#' @export
def_harv <- function(x, species, ...){

    dots <- list(...)
    ct <- dots$ct

    rate <- 0.006 * (ct > 0)
    return(x * rate)
}


#' Default disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Default disturbance function does not require
#'
#' @export
def_disturb <- function(x, species, disturb = NULL, ...){

    if(! is.null(disturb)){
        warning(paste0("default disturbance function does not impact populations",
                       ". Please add your own disturbance function."))
    }
    Pkill <- numeric(length(x))

    return(x* Pkill)
}
