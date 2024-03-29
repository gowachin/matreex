% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_mu_gr.R
\name{make_mu_gr}
\alias{make_mu_gr}
\title{Build IPM for a given species and climate}
\usage{
make_mu_gr(
  species,
  fit,
  climate = NULL,
  mesh = c(m = 700, L = 90, U = 1500),
  BA = 0:200,
  correction = c("none", "constant", "ceiling", "sizeExtremes"),
  stepMu = 0.001,
  level = c(3, 140),
  diag_tresh = 50,
  midbin_tresh = 25,
  mid_level = 5,
  year_delta = 1,
  IsSurv = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{species}{The species names to be registered in the object}

\item{fit}{Fitted model for growth and survival of the species and climate
given. Functions will depend on size and basal area.}

\item{climate}{Climate table for the species.
Optionnal if the species is defined in the package.
The format is climatic variable
in column and different climate by row. An example is in the package with
\code{cliamte_species}.}

\item{mesh}{vector of mesh variables. m is the number of bins, L is the
minimum size and U the maximum size. h will be defined in the function as
\eqn{h <- (U - L) / m}.}

\item{BA}{Vector of basal area to integrate on. Integrating on 0 is important
so use it. Integrating above 200 is absurd.}

\item{correction}{Correction to apply to the IPM matrix for eviction. Choices
constant (default), ceiling, sizeExtremes and none.}

\item{stepMu}{Step between each mu in the species range. This value has effect
simulation. Default is 1e-3. Values below 1e-5 takes longer than classical
matrix integration.}

\item{level}{Number of point to use for integration in a cell during
Gauss-Legendre integration. This value will be divided by 3 since size t is
integrated at level = 3 and size t+1 at level = level/3. single int
(default 420).}

\item{diag_tresh}{Threshold for Gauss-Legendre integration, which a distance
to the diagonal. Number of cell integrated is the number of cell for which
size t+1 - size t is inferior to this threshold. single dbl (default 50).}

\item{midbin_tresh}{Number of cells external to the GL integration to
integrate with the mid bin method.}

\item{mid_level}{Number of point to use for integration in a cell during
mid bin integration.}

\item{year_delta}{Number of year between 2 obersavtion when using this model.
default 1, single int. NOTE : value for dev usage only !}

\item{IsSurv}{Adding survival to the IPM. Set to FALSE is useful to test for
eviction of the model. TRUE by default.}

\item{verbose}{Print message. FALSE by defaul}
}
\description{
Integrate IPM for growth and survival function at a specific climate for a
species on a basal area variation.
}
\details{
The check between climate variables and fitted variable will assert if all
variables in the model are provided expect variables derived from "size"
(size, size2, logsize), "intercept" and "BATOTcomp". An error will be
triggered if the climate variable is missing.

One can desactivate each kind of integration with some treshold values.
A negative value in diag_tresh (ex: -1) will cancel the Gauss-Legendre
integration and a midbin_tresh null value (ex: 0) will cancel the midbin
integration.

This is a working function to test faster integration but it integrate a
ba value
}
