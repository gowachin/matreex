# load("dev/evforest.Rds")
# load("dev/evsim.Rds")

library(purrr)
library(dplyr)
library(ggplot2)
# load_all()
# rm(list = ls())
sim_rdikg <- function(sim, rdi_c = NULL){

    # # DEV ####
    # forest = Jasper_for_Even_dev
    # sim = Jasper_sim_f20_dev
    # rdi_c = NULL
    #********

    sp <- unique(sim$species)

    # TODO if rdi_coef
    if(is.null(rdi_c)){
        rdi_c <- matreex::rdi_coef
        if(! all(sp %in% rdi_c$species)){
            stop("A simulated species is not present in the mackage rdi coefficient. Please provide the rdi_coef argument.")
        }
    }
    rdi_c <- filter(rdi_c, species %in% sp)

    sim <- sim %>%
        filter(size > 0, ! equil) %>%
        select(- equil, -mesh) %>%
        mutate(size = (size / 10)^2)
    # Compute rdi ####
    rdi_sp <- sim %>%
        filter(var == "n") %>%
        left_join(rdi_c, by = "species") %>%
        group_by(species, time) %>%
        group_modify(~ data.frame(rdi = RDI(
            x = .x$value,
            RDI_int = unique(.x$intercept), RDI_slo = unique(.x$slope),
            meshcm2 = .x$size)))# %>%

    rdi_val <- rdi_sp %>%
        group_by(time) %>%
        summarise(rdi = sum(rdi), species = "All") %>%
        select(species, time, rdi) %>%
        bind_rows(rdi_sp) %>%
        pivot_longer(rdi, names_to = "var")

    # Compute Kg
    tmp <- sim %>%
        filter(var %in% c("n", "h")) %>%
        tidyr::pivot_wider(names_from = var, values_from = value)
    kg_val <- tmp %>%
        group_by(time, species) %>%
        mutate(X = n + h) %>%
        summarise(
            surfx = drop(size %*% X),
            tx = sum(X),
            surfcut = drop(size %*% h),
            tcut = sum(h),
            .groups = "drop") %>%
        group_by(time) %>%
        summarise(
            Dg2 = sum(surfx) / sum(tx),
            Dgcut2 = sum(surfcut) / sum(tcut),
            Kg = Dgcut2 / Dg2)  %>%
        replace_na(list(Kg = 0)) %>%
        mutate(species = "All") %>%
        select(species, time, Dg2, Dgcut2, Kg) %>%
        pivot_longer(cols = c("Dg2", "Dgcut2", "Kg"), names_to = "var")

    # output all
    res <- bind_rows(rdi_val, kg_val)

    return(res)
}


# x <- sim_rdikg(sim = Jasper_sim_f20_dev,
#          rdi_c = NULL)
#
# Jasper_sim_f20_dev %>%
#     filter(var %in% c("N", "BAsp", "H")) %>%
#     select(species, time, var, value) %>%
#     bind_rows(x) %>%
#     filter(var != "Dgcut2") %>%
#     ggplot(aes(x = time, y = value, color = species)) +
#     geom_line() + #geom_point() +
#     facet_wrap(~ var, scales = "free_y") +
#     geom_hline(yintercept = 0.9, color = "red") +
#     geom_hline(yintercept = 0.9, color = "green") +
#      NULL
#
#
# t <- Jasper_sim_f20_dev %>%
#     filter(var == "H", value > 0) %>%
#     pull(time) %>% unique()
#
# x <- x %>%
#     filter(time %in% t, species == "All", var %in% c("Kg", "rdi")) %>%
#     pivot_wider(names_from = var, values_from = value)
#
# for(i in 1:nrow(x)){
#     cat(sprintf(
#         "Kg : %.4f | RDI : %.4f \n",
#         x$Kg[i], x$rdi[i])
#     )
# }

#' yeay ça marche mais du coup question, est-ce que c'est valable pour la version d'origine ?
# N <- length(mesh)
# Dg1 <- mesh^2 %*% (x*(1-Pc)) / N
# RDI <- sum(x*(1-Pc)) / exp(RDIcoef$rqIntercept + RDIcoef$rqSlope/2*log(Dg1))


# load_all()
# set.seed(42) # The seed is here for initial population random functions.
# Picea_sim_f20 <- sim_deter_forest(
#     Picea_for_Even,
#     tlim = 100,
#     equil_time = 100, equil_dist = 10, equil_diff = 1,
#     harvest = "Even", targetRDI = 0.6, targetKg = 0.1,
#     final_harv = 150,
#     SurfEch = 0.03,
#     verbose = TRUE
# )
# # Picea_sim_f20  %>%
#     # dplyr::filter(var %in% c("BAsp", "N", "H"), ! equil) %>%
#     # ggplot(aes(x = time, y = value)) +
#     # facet_wrap(~ var, scales = "free_y") +
#     # geom_line(linewidth = .2) + geom_point(size = 0.4)
#
# x <- sim_rdikg(sim = Picea_sim_f20,
#                 rdi_c = NULL)
