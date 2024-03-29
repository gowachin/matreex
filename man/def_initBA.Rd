% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_species.R
\name{def_initBA}
\alias{def_initBA}
\title{Init population at BA}
\usage{
def_initBA(BA = 1, fun = c("def_init", "def_init_even"))
}
\arguments{
\item{BA}{Basal area targeted. This single value must be above 0 but can be
very close (minimal accepted value is 1e-10)}

\item{fun}{Function to modify, single chr in choices.}
}
\value{
Function similar to def_init
}
\description{
This function modify the def_init function to start at a given BA with
the same process of random distribution.
}
\seealso{
Other functions for initiating species population during simulation: 
\code{\link{def_init_even}()},
\code{\link{def_init_k}()},
\code{\link{def_init}()},
\code{\link{species}()}
}
\concept{functions for initiating species population during simulation}
