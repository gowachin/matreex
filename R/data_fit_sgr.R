#' Species fited models dataset.
#'
#' Reformated fitted models from Kunstler et al. 2021. Each species
#' is fitted for growth, survival and recruitment along climatic conditions
#' based on sgdd and wai. All objects were modified with the oldfit2fit function
#' from this package.
#'
#' The fit_species dataset contains the vector of species stored with the package.
#' The lag_species dataset contains the vector of default delay for each species.
#'
#' @details
#'
#' The lag was computed by Kunstler and missing values for Betula, 
#' Juniperus_thurifera, Prunus_padus, Quercus_faginea and Quercus_pyrenaica
#' are set with the mean of other species. More details is given in 
#' https://github.com/gowachin/matreex/issues/10
#'
#' @name fit_data
#' @aliases fit_species
#'
#' @source https://doi.org/10.1111/1365-2745.13533
#'
"fit_species"

#' @rdname fit_data
"lag_species"

#' @rdname fit_data
"fit_Abies_alba"

#' @rdname fit_data
"fit_Acer_campestre"

#' @rdname fit_data
"fit_Acer_pseudoplatanus"

#' @rdname fit_data
"fit_Alnus_glutinosa"

#' @rdname fit_data
"fit_Betula"

#' @rdname fit_data
"fit_Carpinus_betulus"

#' @rdname fit_data
"fit_Fagus_sylvatica"

#' @rdname fit_data
"fit_Fraxinus_excelsior"

#' @rdname fit_data
"fit_Juniperus_thurifera"

#' @rdname fit_data
"fit_Larix_decidua"

#' @rdname fit_data
"fit_Picea_abies"

#' @rdname fit_data
"fit_Pinus_halepensis"

#' @rdname fit_data
"fit_Pinus_nigra"

#' @rdname fit_data
"fit_Pinus_pinaster"

#' @rdname fit_data
"fit_Pinus_pinea"

#' @rdname fit_data
"fit_Pinus_sylvestris"

#' @rdname fit_data
"fit_Pinus_uncinata"

#' @rdname fit_data
"fit_Populus_tremula"

#' @rdname fit_data
"fit_Prunus_padus"

#' @rdname fit_data
"fit_Quercus_faginea"

#' @rdname fit_data
"fit_Quercus_ilex"

#' @rdname fit_data
"fit_Quercus_petraea"

#' @rdname fit_data
"fit_Quercus_pubescens"

#' @rdname fit_data
"fit_Quercus_pyrenaica"

#' @rdname fit_data
"fit_Quercus_robur"

#' @rdname fit_data
"fit_Quercus_suber"

#' @rdname fit_data
"fit_Salix_caprea"
