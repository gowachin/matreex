% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/disturbance.R
\name{disturb_fun_mixt}
\alias{disturb_fun_mixt}
\title{Disturbance function}
\usage{
disturb_fun_mixt(x, species, disturb = NULL, ...)
}
\arguments{
\item{x}{population state distribution at time t}

\item{species}{The species class object of interest to get mesh and RDIcoef
values from. RDIcoef is a one line dataframe with RDI coefficient for one
species.}

\item{disturb}{Disturbance parameters. Highly depend on the disturbance
impact parameters given to the species.}

\item{...}{Not used in this case.
\describe{
\item{qmd}{Forest Quadratic Mean Diameter}
}}
}
\description{
Disturbance function
}
\details{
Delayed mesh cells takes the value of the minimal mesh size.
}
