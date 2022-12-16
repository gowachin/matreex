
devtools::load_all()
ipm <- make_IPM("Abies_alba",
                climate = drop(as.matrix(
                    subset(climate_species, sp == "Abies_alba" & N == 2, select = -c(N, sp))
                    )),
                clim_lab = "optim", fit = fit_Abies_alba,
                mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Abies_alba) * 1.1),
                BA = 0:10, verbose = TRUE)

object <- ipm

print.ipm <- function(x, ...){
    summary(x)
}

summary.ipm <- function(object, ...){

    info <- object$info
    int <- object$int
    dim <- dim(object$IPM[[1]])
    BA <- object$BA
    mesh <- object$mesh

    res <- sprintf(
        paste0(
            "IPM object for species %s at climate '%s' \n",
            "Integation was done on a mesh from %.2f to %.2f with %.0f cells between %.0f and %.0f of BA. \n",
            "Gauss-Legendre was used on %.0f cells with %.0f x %.0f levels (size at t * t+1) \n",
            "Midbin was used on %.0f cells with level %.0f. \n",
            "The correction was %s ,the IPM does %s contain survival and the recruitment delay is %.0f. \n"
        ),
        info["species"],  info["clim_lab"],
        mesh[1], tail(mesh, 1), length(mesh), BA[1], tail(BA, 1),
        int["gl_tresh"], int["gl1"], int["gl2"],
        int["mb_tresh"], int["mid_level"],
        info["correction"], ifelse(as.logical(info["surv"]), "", "not"),
        as.numeric(info["delay"])
    )

    message(res)

    invisible(object)
}

summary.mu_gr <- function(object, ...){
    NULL
}

summary.species <- function(object, ...){
    NULL
}

summary.forest <- function(object, ...){
    NULL
}
