#' @param x a named values for a climatic variable.
#' @param inv_null If TRUE, the inverse variable is \code{varb = 1/(var+1)}, in
#' case of var can take 0 for value.
#' Example : \code{x = c(sgdd = 2000)}
#'
#' @noRd
value2sqb <- function(x, inv_null = FALSE){

    nms <- names(x)
    # res <- numeric(3)
    res <- as.data.frame(matrix(ncol = 3, nrow = nrow(x)))
    names(res) <- c(nms, paste0(nms, c("2", "b")))
    res[, nms] <- x[,nms]
    res[, paste0(nms, "2")] <- x[,nms]^2
    res[, paste0(nms, "b")] <- 1/(x[,nms] + inv_null)

    res
}

#' @param climate Named vector of climatic variables.
#' @param inv_null If TRUE, the inverse variable is \code{varb = 1/(var+1)}, in
#' case of var can take 0 for value.
#' Example : \code {climate = c(sgdd = 2000, wai = 0.16) ; inv_null = c(sgdd = FALSE, wai = TRUE)}
#'
#' @noRd
expand_clim <- function(climate, inv_null){

    # assertNumeric(climate, any.missing = FALSE)
    assertLogical(inv_null, any.missing = FALSE)

    nms <- names(climate)
    res <- vector("list", length(climate))
    for(x in seq_along(nms)){
        res[[x]] <- value2sqb(climate[nms[x]], inv_null[nms[x]])
    }
    res <- do.call("cbind", res)

    return(res)
}

#' Build climatic gradient
#'
#' @param clim_table climate table from data climate_species, filtered on a species
#' @param start_clim strating climate in opti, hot or cold
#' @param end_clim same as start_clim
#' @param times duration of the simulation without extrapolation
#' @param noise function to draw noise from. Default is rnorm. Shoud have a
#' n value to draw and an sd parameter.
#' @param sigma standard deviation value for the noise, for the two climatic variable
#' sgdd and wai.
#' @param seed for the random effect of the noise function
#' @param extrapolate extrapolation of the climatic tendancy in percent of the
#' times. Example, for 100 year and an extrapolation of 50% will add 50 years of
#' simulations beyond the last climate.
#'
#' @importFrom stats rnorm
#' @details
#' The climate will be linear between the two climate, using a seq() function.
#' clim_table must encode the climate with numeric 1, 2 and 3 respectively
#' values for hot, opti and cold.
#'
#' @export
clim_gradient <- function(clim_table,
                          start_clim = c("opti", "hot", "cold"),
                          end_clim = c("opti", "hot", "cold"), times = 100,
                          noise = rnorm, sigma = c(1, 0.001), seed = 42,
                          extrapolate = 0){

    match.arg(start_clim)
    match.arg(end_clim)

    N <- NULL

    assertNumber(extrapolate, lower = 0)
    extrapolate <- extrapolate / 100

    start_clim <- switch(start_clim, opti = 2, hot = 1, cold = 3)
    end_clim <- switch(end_clim, opti = 2, hot = 1, cold = 3)

    beg <- subset(clim_table, N == start_clim, select = c("sgdd", "wai"))
    end <- subset(clim_table, N == end_clim, select = c("sgdd", "wai"))

    if(extrapolate > 0){
        end <- end + ((end - beg) * extrapolate)
        times <- times + times * extrapolate
    }

    # add the sigma
    pre_table <- rbind(beg, end, sigma)

    set.seed(seed = seed)
    basic_clim <- apply(
        pre_table, 2,
        function(x, t, foo){
            seq(from = x[[1]], to = x[[2]], length.out = t) + foo(t, sd = x[[3]])
        },
        t = times, foo = noise)
    basic_clim <- as.data.frame(basic_clim)

    # Now expand the climate
    res <- expand_clim(basic_clim, inv_null = c(sgdd = FALSE, wai = TRUE))

    return(res)
}
