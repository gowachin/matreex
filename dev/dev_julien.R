#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#' @name disturbance.R
#' @description R script to disturb a population simulated by IPM
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Et du coup, j'en profites pour t'envoyer le rdata contenant les
# paramètres pour l'équation de sensibilité aux perturbations, et
# un script qui montre l'utilisation du fichier pour appliquer
# une perturbation à une population donnée. J'ai codé une fonction
# à qui tu indiques la nature, l'intensité et la durée de la
# perturbation et qui l'applique à la dernière itération d'une simulation
# (dans un contexte monospécifique). 
#
#
# Les autres fonctions, pas la peine de t'en préoccuper, c'est juste
# le début du taf que j'ai commencé pour faire tourner IPM et simulations
# d'un seul coup pour beaucoup d'espèces en passant par des listes.
# Pour toi, les fonctions intéressantes sont: 
#
# - load_and_format_dist.param: pour formatter les paramètres
# de perturbations pour avoir une valeur moyenne par espèce et par
# perturbation (au lieu de toutes les itérations posterior)
#
# - disturb.population pour perturber une population à la dernière
# itération de la simulation. 
#
# En tout cas encore merci et du coup bon courage pour cette
# histoire de découplage survie / mortalité extrême !


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Clear environment and load packages and functions -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rm(list = ls())
packages <- c("dplyr", "ggplot2", "matreex", "tidyr")
lapply(packages, require, character.only = TRUE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Functions ---------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to load and format disturbance parameters per species
#' @param disturbance_parameters_file file containing parameters value
#' @param calc.type character indicating whether to use mean parameter value ("mean")
#'                  or all posterior iterations ("all")
load_and_format_dist.param <- function(disturbance_parameters_file, calc.type = "mean"){

    # Load the rdata file
    load(disturbance_parameters_file)

    # Remove unused disturbance types and species
    out <- param_per_iteration %>%
        filter(!(species %in% c("Other broadleaf", "Other conifer"))) %>%
        filter(disturbance %in% c("storm", "fire", "biotic"))

    # If the calculation type is mean, average parameters over all mcmc iteration
    if(calc.type == "mean"){
        out <- out %>%
            gather(key = "parameter", value = "value", "a0", "a1", "b", "c", "dbh.intercept",
                   "dbh.slope", "logratio.intercept", "logratio.slope") %>%
            group_by(disturbance, species, parameter) %>%
            summarize(value.mean = mean(value)) %>%
            spread(key = "parameter", value = "value.mean")
    }

    # Return output
    return(out)
}


#' #' Function to load demographic parameters for several species
#' #' @param species.names character vector of species ("Genus_species" format)
#' load_param_demo <- function(species.names){
#'
#'     # Initialize list of fits (one element per species)
#'     fit.list <- list()
#'
#'     # Loop on all species to load and store demographic parameters
#'     for(i in 1:length(species.names)){
#'         eval(parse(text=paste0("fit.list$", species.names[i], " <- fit_", species.names[i])))
#'     }
#'
#'     # Return the list
#'     return(fit.list)
#'
#' }



#' #' Function to load species climate data and format it as input for the IPM
#' #' @param species.name Name of the species foor which to load climate (Genus_species format)
#' #' @param N.ref Integer specifying whether to keep species cold margin (1), optimum (2) or hot margin (3)
#' load_and_format_climate_ipm <- function(species.name, N.ref = 2){
#'
#'     # Load the data
#'     data("climate_species")
#'
#'     # Format to fit IPM
#'     out  <- drop(as.matrix(
#'         climate_species %>%
#'             filter(N == N.ref & sp == species.name) %>%
#'             dplyr::select(-sp)
#'     ))
#'
#'     # Return formatted vector
#'     return(out)
#' }



#' #' Fit IPM for several species and one climate
#' #' @param fit.list List containing the demographic parameters of each species
#' #' @param climate.ref Climate to use for each fit of the IPM
#' #' @param clim_lab.in Name of the reference climate
#' #' @param mesh.m numeric indicating the mesh size
#' #' @param mesh.L numeric indicating the lower bound of size
#' #' @param BA.max maximum basal area for the integration
#' make_IPM_multispecies <- function(fit.list, climate.ref, clim_lab.in,
#'                                   mesh.m = 700, mesh.L = 100, BA.max = 200){
#'
#'     # Identify the species names
#'     species.names <- names(fit.list)
#'
#'     # Initialize the list of IPM
#'     IPM.list <- list()
#'
#'     # Loop on all species
#'     for(i in 1:length(species.names)){
#'         # Print the species
#'         print(paste0("fit IPM for species ", i, "/", length(species.names), " : ", gsub("\\_", "\\ ", species.names[i])))
#'         # Make the IPM
#'         ipm.i <- make_IPM(
#'             species = species.names[i], climate = climate.ref, fit =  fit.list[[i]],
#'             clim_lab = clim_lab.in,
#'             mesh = c(m = mesh.m, L = mesh.L, U = as.numeric(fit.list[[i]]$info[["max_dbh"]]) * 1.1),
#'             BA = 0:BA.max, verbose = TRUE
#'         )
#'         # Add to the list
#'         eval(parse(text=paste0("IPM.list$", species.names[i], " <- ipm.i")))
#'     }
#'
#'     # Return the list generated
#'     return(IPM.list)
#' }
#'
#'
#'
#' #' Generate a list of species object from the list of IPM
#' #' @param IPM.list list of fitted IPM
#' #' @param f.init Function to initialize basal area
#' generate_species.list <- function(IPM.list, f.init = def_initBA(20)){
#'
#'     # Names of the species
#'     species.names <- names(IPM.list)
#'
#'     # Initialize the list of species
#'     species.list <- list()
#'
#'     # Loop on all IPM (and thus on all species)
#'     for(i in 1:length(species.names)){
#'         # Generate species object for species i
#'         species.i <- new_species(IPM.list[[i]], init_pop = f.init,
#'                                  harvest_fun = Uneven_harv)
#'         # Add it to the list
#'         eval(parse(text=paste0("species.list$", species.names[i], " <- species.i")))
#'     }
#'
#'     # Return the output list
#'     return(species.list)
#' }



#' Function to disturb the population at the last iteration of a simulation
#' @param sim simulation object
#' @param disturbance.in which disturbance to apply
#' @param intensity.in intensity (between 0 and 1) of the disturbance
#' @param duration.disturbance how long should the disturbance kill trees
#' @param disturbance_parameters Parameters averaged over all iterations
disturb.population <- function(sim, disturbance.in, intensity.in,
                               duration.disturbance, disturbance_parameters){

    # Identify the last year of simulation
    disturbance.time <- dim(sim)[2] - 1

    # Identify the species in the simulation
    species.in <- unique(gsub("\\..+", "", rownames(sim)))

    # Extract vector to disturb
    pop.to.disturb <- tree_format(sim) %>%
        # Restrict to last year post disturbance and to tree number per size class
        filter(time == disturbance.time & var == "m" & species %in% species.in) %>%
        dplyr::select(-equil) %>%
        distinct()

    # Calculate quadratic diameter at time of disturbance
    dqm.disturbance.time <- (pop.to.disturb %>%
                                 # group_by(size) %>%
                                 # summarize(value = sum(value)) %>%
                                 # ungroup() %>%
                                 summarize(dqm = sqrt(sum(size*size*value)/sum(value))))$dqm

    # Disturb the population with the mean of parameter values
    disturbed.proba <- pop.to.disturb %>%
        # Add disturbance parameters
        mutate(species = gsub("\\_", "\\ ", species),
               disturbance = disturbance.in) %>%
        left_join(disturbance_parameters, by = c("species", "disturbance")) %>%
        # Calculate variables used to calculate probability
        mutate(dqm = dqm.disturbance.time,
               logratio = log(size/dqm),
               dbh.scaled = dbh.intercept + size*dbh.slope,
               logratio.scaled = logratio.intercept + logratio*logratio.slope,
               I = intensity.in,
               t = duration.disturbance) %>%
        # Calculate probability of survival from disturbance per size class
        mutate(p = (1 - plogis(a0 + a1*logratio.scaled + b*I^(c*dbh.scaled)))^t) %>%
        dplyr::select(size, value, p) %>%
        mutate(distribution.disturbed = value*p)

    # Length of the mesh
    mesh <- max(as.numeric(gsub(".+\\.m", "", rownames(sim))), na.rm = TRUE)

    # Species distribution disturbed
    distr.before.disturbance <- sim[grep(paste0(species.in,".m"), rownames(sim)),
                                    grep("eqt", colnames(sim))]
    distr.disturbed <- distr.before.disturbance
    distr.disturbed[c(1:mesh)] <- disturbed.proba$distribution.disturbed

    # Return disturbed population
    return(distr.disturbed)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Load parameters, conduct simulation and disturb population ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## - Disturbance parameters

# Name of the file containing disturbance parameters
disturbance_parameters_file <- "data/parameters_disturbance.Rdata"
# Formatted parameters for disturbance
disturbance_parameters <- load_and_format_dist.param(disturbance_parameters_file)




## - Demographic parameters

# Species for which to conduct simulations
species.names <- "Fagus_sylvatica"
# Get demographic parameters
fit.list <- load_param_demo(species.names)



## - IPM

# Get the optimum climate as reference climate
climate.ref <- load_and_format_climate_ipm(species.names[1])
# Fit IPM for all species in species.names
IPM.list <- make_IPM_multispecies(fit.list, climate.ref, clim_lab.in = "climate_Fag.syl")


## - Make simulations

# Create a list of species object (one element per species in species.names)
species.list <- generate_species.list(IPM.list)
# Make simulation with the first species
sim <- sim_deter_forest(
    species.list[[species.names[1]]], tlim = 1500, equil_time = 1500, equil_dist = 1,
    harvest = "default", targetBA = 40, SurfEch = 0.03, correction = "cut", verbose = TRUE
)




## - Disturb the population at the last iteration of the simulation with
#    a storm of intensity 0.7 over 5 years
test <- disturb.population(
    sim, disturbance.in = "storm", intensity.in = 0.7,
    duration.disturbance = 5, disturbance_parameters)



