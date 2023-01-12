
devtools::load_all()
ipm <- make_IPM("Abies_alba",
                climate = drop(as.matrix(
                    subset(climate_species, sp == "Abies_alba" & N == 2, select = -c(N, sp))
                    )),
                clim_lab = "optim", fit = fit_Abies_alba,
                mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Abies_alba) * 1.1),
                BA = 0:10, verbose = TRUE)

sp <- species(ipm, init_pop = def_init)


object <- sp

print.ipm <- function(x, ...){
    summary(x)
}

#' IPM summary
#'
#' @param object IPM object.
#' @param ...
#' \describe{
#'  \item{fline}{Print first line. TRUE by default.}
#' }
#'
#' @keywords Internal
#' @export
summary.ipm <- function(object, ...){

    info <- object$info
    int <- object$int
    dim <- dim(object$IPM[[1]])
    BA <- object$BA
    mesh <- object$mesh

    dots <- list(...)
    fline <- dots$fline
    fline <- ifelse(is.null(fline), TRUE, fline)

    res <- sprintf(
        paste0(
            ifelse(fline, "%s%s",
                   "IPM object for species %s at climate '%s' \n\n"
                   ),
            "Integation was done on a mesh from %.2f to %.2f with %.0f cells for BA between %.0f and %.0f. \n",
            "Gauss-Legendre was used on %.0f cells with %.0f x %.0f levels (size at t * t+1) \n",
            "Midbin was used on %.0f cells with %.0f levels. \n",
            "The correction was %s, the IPM does %scontain survival and the recruitment delay is %.0f. \n"
        ),
        ifelse(fline, "", info["species"]),
        ifelse(fline, "", info["clim_lab"]),
        mesh[1], tail(mesh, 1), length(mesh), BA[1], tail(BA, 1),
        int["gl_tresh"], int["gl1"], int["gl2"],
        int["mb_tresh"], int["mid_level"],
        info["correction"], ifelse(as.logical(info["surv"]), "", "not "),
        as.numeric(info["delay"])
    )

    message(res)

    invisible(object)
}

summary(object, sum_ipm = FALSE)

summary.mu_gr <- function(object, ...){
    NULL
}

#' species summary
#'
#' @param object species object.
#' @param ...
#' \describe{
#'  \item{sum_ipm}{Print IPM summary. TRUE by default.}
#' }
#'
#' @keywords Internal
#' @export
summary.species <- function(object, ...){

    info <- object$info
    harv_lim <- object$harv_lim
    rdi_coef <- object$rdi_coef
    disturb_coef <- object$disturb_coef

    dots <- list(...)
    sum_ipm <- dots$sum_ipm

    res <- sprintf(
        paste0(
            "IPM object for species %s at climate '%s' \n\n",
            "For Uneven harvest, dth = %.1f, dha = %.1f and hmax = %.2f.\n",
            "rdi_coef are %sdefined for this species. Even harvest %spossible\n",
            "disturb_coef are %sdefined for this species. Disturbance models %spossible\n"
        ),
        info["species"],  info["clim_lab"],
        harv_lim["dth"], harv_lim["dha"], harv_lim["hmax"],
        ifelse(is.null(rdi_coef), "not ", ""),
        ifelse(is.null(rdi_coef), "not ", ""),
        ifelse(is.null(disturb_coef), "not ", ""),
        ifelse(is.null(disturb_coef), "not ", "")
    )

    message(res)

    if(is.null(sum_ipm) || sum_ipm){
        summary(object$IPM, fline = TRUE)
    }
    invisible(object)
}

summary.forest <- function(object, ...){
    NULL
}
