#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#'  \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#'
#' @details
#' Delayed mesh cells takes the value of the minimal mesh size.
#'
#' @export
disturb_fun <- function(x, species, disturb = NULL, ...){

    dots <- list(...)
    qmd <- dots$qmd
    size <- species$IPM$mesh
    coef <- species$disturb_coef
    if(any(disturb$type %in% coef$disturbance)){
        coef <- subset(coef, disturbance == disturb$type)
    } else {
        stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                     sp_name(species), disturb$type))
    }

    # edits for delay
    size[size == 0] <- min(size[size !=0])

    logratio <-  log(size / qmd)
    dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
    logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
    Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                        coef$b * disturb$intensity ^(coef$c * dbh.scaled))

    return(x* Pkill) # always return the mortality distribution
}


#' Disturbance function
#'
#' @inheritParams disturb_fun
#'
#' @details
#' Delayed mesh cells takes the value of the minimal mesh size.
#'
#' @export
disturb_fun_mixt <- function(x, species, disturb = NULL, ...){

    dots <- list(...)
    qmd <- dots$qmd
    coni_perc <- dots$perc_coni # New
    size <- species$IPM$mesh
    coef <- species$disturb_coef
    if(any(disturb$type %in% coef$disturbance)){
        coef <- subset(coef, disturbance == disturb$type)
    } else {
        stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                     sp_name(species), disturb$type))
    }

    # edits for delay
    size[size == 0] <- min(size[size !=0])

    logratio <-  log(size / qmd)
    dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
    logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
    Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                        coef$b * disturb$intensity ^(coef$c * dbh.scaled))

    # Jasper equation #*********************************************************
    if(species$info["type"] == "Coniferous" && disturb$type == "biotic"){
        # just a message to check for me
        message(sprintf("sp : %s | coni : %.3f | modif : %.2f",
                        species$info["species"], coni_perc,
                        (0.4767 + 0.6015 * coni_perc)))
        Pkill <- Pkill * (0.4767 + 0.6015 * coni_perc)
    }
    # End of Edit #*************************************************************


    Pkill <- pmin(Pkill, 1)

    return(x* Pkill) # always return the mortality distribution
}
