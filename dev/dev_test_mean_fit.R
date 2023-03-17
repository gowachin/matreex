# Test difference between mean_fit models ####

if(!interactive()){
    stop(paste("This script copy files and depend on external relative paths.",
               "Please run it interactively."))
}


## Library ####
library(here)
library(dplyr)
library(tibble)
library(matreex)
library(purrr)
library(ggplot2)

## Functions ####
compt_clim_resp <- function(
        mean, models, climate, mesh, BATOTSP, BATOTNonSP, minmax = FALSE){

    inv <- models[[1]]$sv$family$linkinv

    # Compute mean answer
    mean_gr <- matreex:::exp_sizeFun(mean$gr$params_m, climate) %>%
        do.call(list(size = mesh))
    mean_sv <- matreex:::exp_sizeFun(mean$sv$params_m, climate) %>%
        do.call(list(size = mesh)) %>%
        inv() %>% `-`() %>% `+`(1)
    l_gr <- purrr::map(models, ~ matreex:::exp_sizeFun(.x$gr$params_m, climate)) %>%
        purrr::map(~ purrr::exec(.x, mesh)) %>%
        `names<-`(paste0("gr_m",1:length(models))) %>%
        do.call(what = "cbind")

    l_sv <- purrr::map(models, ~ matreex:::exp_sizeFun(.x$sv$params_m, climate)) %>%
        purrr::map(~ 1 - inv(purrr::exec(.x, mesh))) %>%
        `names<-`(paste0("sv_m",1:length(models))) %>%
        do.call(what = "cbind")

    if(minmax){
        l_gr <- t(apply(l_gr, 1, range)) %>%
            `colnames<-`(c("gr_min", "gr_max"))
        l_sv <- t(apply(l_sv, 1, range)) %>%
            `colnames<-`(c("sv_min", "sv_max"))
    }

    res <- dplyr::bind_cols(
        clim = switch (climate["N"], "hot", "opti", "cold"),
        species = mean$info["species"], mesh = mesh,
        l_gr, l_sv, gr_mean = mean_gr, sv_mean = mean_sv) %>%
        tidyr::pivot_longer(cols = -c("clim", "species", "mesh"),
                            names_to = "var", values_to = "value") %>%
        tidyr::separate(var, sep = "_", into = c("vital", "models"))


    return(res)
}

compt_clim_rec <- function(
        mean, models, climate, mesh, BATOTSP, BATOTNonSP, minmax = FALSE){

    n <- length(models)

    mean_rec <- list(matreex:::exp_recFun(mean$rec$params_m, climate)) %>%
        purrr::cross2(BATOTSP) %>%
        `names<-`(paste0("rec_mean_BA", BATOTSP)) %>%
        map(~ do.call(what = purrr::exec, c(.x, list(BATOTNonSP, mesh)))) %>%
        map_dfc(~ sum(.x))

    l_rec <- purrr::map(models, ~ matreex:::exp_recFun(.x$rec$params_m, climate)) %>%
        purrr::cross2(BATOTSP) %>%
        `names<-`(
            paste0("rec_m", 1:n, "_BA", rep(BATOTSP, each = n))
        ) %>%
        map(~ do.call(what = purrr::exec, c(.x, list(BATOTNonSP, mesh)))) %>%
        map_dfc(~ sum(.x))

    res <- dplyr::bind_cols(
        clim = switch (climate["N"], "hot", "opti", "cold"),
        species = mean$info["species"],
        l_rec, mean_rec) %>%
        tidyr::pivot_longer(cols = -c("clim", "species"),
                            names_to = "var", values_to = "value") %>%
        tidyr::separate(var, sep = "_", into = c("vital", "models", "BA")) %>%
        dplyr::mutate(BA = as.numeric(sub("^BA", "",BA)))

    if(minmax){
        res <- res %>%
            dplyr::mutate(new_mod = ifelse(models == "mean", "mean", "models")) %>%
            dplyr::group_by(new_mod, BA, clim, species, vital) %>%
            dplyr::summarise(rec_min = min(value), rec_max = max(value),
                             mean = mean(value), .groups ="drop")

    }

    return(res)
}


mix_clim <- function(
        path, files, species, climate_sp,
        BATOTcomp,BATOTSP, BATOTNonSP, minmax = FALSE){

    climate <- dplyr::filter(climate_sp, sp == species)


    dir.create(here("output/", species)) # create tmp dir with name
    file.copy(here(path, "output", files[species]), # rename old by copy
              here("output", species, "fit_sgr_all.Rds"))

    mean <- old_fit2fit(species, path = here(), mean = TRUE)
    mesh <- seq(90, get_maxdbh(mean)*1.1, length.out = 700)

    models <- readRDS(here("output", species, "fit_sgr_all.Rds"))

    res <- vector(mode = "list", length = 3)
    for(i in 1:3){
        clim <- c(unlist(climate[i,-c(7,8,10,11)]), BATOTcomp = BATOTcomp)
        res[[i]] <- compt_clim_resp(
            mean, models, clim, mesh, BATOTSP, BATOTNonSP, minmax = minmax
        )
    }

    unlink(here("output/", species), recursive = TRUE) # remove tmp dir

    res <- do.call(rbind, res)
    levels(res$clim) <- c("cold", "opti", "hot")
    return(res)
}

rec_mix_clim <- function(
        path, files, species, climate_sp,
        BATOTcomp,BATOTSP, BATOTNonSP, minmax = FALSE){

    climate <- dplyr::filter(climate_sp, sp == species)

    dir.create(here("output/", species)) # create tmp dir with name
    file.copy(here(path, "output", files[species]), # rename old by copy
              here("output", species, "fit_sgr_all.Rds"))

    mean <- old_fit2fit(species, path = here(), mean = TRUE)
    mesh <- seq(90, get_maxdbh(mean)*1.1, length.out = 700)

    models <- readRDS(here("output", species, "fit_sgr_all.Rds"))

    res <- vector(mode = "list", length = 3)
    for(i in 1:3){
        clim <- c(unlist(climate[i,-c(7,8,10,11)]))
        res[[i]] <- compt_clim_rec(
            mean, models, clim, mesh, BATOTcomp, BATOTNonSP, minmax = minmax
        )
    }

    unlink(here("output/", species), recursive = TRUE) # remove tmp dir

    res <- do.call(rbind, res)
    levels(res$clim) <- c("cold", "opti", "hot")
    return(res)
}


plot_ribbon <- function(data){
    data %>%
        tidyr::pivot_wider(names_from = models, values_from = value) %>%
        # mutate_if(is.numeric, ~replace(., . <= 0, NA)) %>%
        ggplot(aes(x = mesh, y = mean, ymin = min, ymax = max, color = clim, fill = clim)) +
        geom_ribbon(alpha = 0.2, linetype = 0) +
        geom_line( linewidth=0.5) +
        facet_wrap(species ~ vital, scales = "free") +
        scale_colour_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
        scale_fill_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
        scale_alpha_discrete(range = c(0.08, 1), name = "mean") +
        NULL
}

plot_lines <- function(data){
    data %>%
        # filter(value > 0) %>%
        mutate(alpha = models == "mean") %>%
        ggplot(aes(x = mesh, y = value, group = interaction(models,clim),
                   color = clim,  alpha = alpha)) +
        geom_line( linewidth=0.5) +
        facet_wrap(species ~ vital, scales = "free") +
        scale_colour_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
        scale_alpha_discrete(range = c(0.08, 1), name = "mean") +
        NULL
}


## Common dataset ####
path <- "~/IPM" # personnal place, dataset comes from ZENODO

files <- list.files(here(path, "output"))
files <- files[grepl("fit_sgr_all_.*.Rds", files)]
species <- sub("^(fit_sgr_all_)(.*)(.Rds)$", "\\2", files)
uspecies <- sub(" ", "_", species) # correct name
names(files) <- uspecies

climate_sp <- matreex::climate_species
BATOTNonSP <- 0

## Gr and Sv ####
start <- Sys.time()
### Compute all models relation to mesh ####

BATOTcomp <- BATOTSP <- 20

l_species <- l_species_lines <- vector("list", length(species)) %>%
    `names<-`(uspecies)

for(sp in uspecies){
    cat("doing :", sp, which(uspecies == sp), "/", length(uspecies), "\n")
    l_species[[sp]] <- mix_clim(path, files, sp, climate_sp,
                          BATOTcomp,BATOTSP, BATOTNonSP, minmax = TRUE)
}
data <- do.call(rbind, l_species)

for(sp in uspecies){
    cat("doing :", sp, which(uspecies == sp), "/", length(uspecies), "\n")
    l_species_lines[[sp]] <- mix_clim(path, files, sp, climate_sp,
                                BATOTcomp,BATOTSP, BATOTNonSP, minmax = FALSE)
}
data_lines <- do.call(rbind, l_species_lines)

### Graphic output ####
png(width = 1800, height = 800, filename = "dev/gr_ribbon.png")
dplyr::filter(data, vital =="gr") %>% plot_ribbon()
dev.off()
png(width = 1800, height = 800, filename = "dev/sv_ribbon.png")
dplyr::filter(data, vital =="sv") %>% plot_ribbon()
dev.off()

png(width = 1800, height = 800, filename = "dev/gr_lines.png")
dplyr::filter(data_lines, vital =="gr") %>% plot_lines()
dev.off()
png(width = 1800, height = 800, filename = "dev/sv_lines.png")
dplyr::filter(data_lines, vital =="sv", value != 0) %>% plot_lines()
dev.off()
Sys.time() - start
#Time difference of 5.596317 mins


## Rec ####

BATOTcomp <- BATOTSP <- 1:100

l_species <- l_species_lines <- vector("list", length(species)) %>%
    `names<-`(uspecies)

for(sp in uspecies){
    cat("doing :", sp, which(uspecies == sp), "/", length(uspecies), "\n")
    l_species[[sp]] <- rec_mix_clim(path, files, sp, climate_sp,
                                BATOTcomp,BATOTSP, BATOTNonSP, minmax = TRUE)
}
rec_data <- do.call(rbind, l_species)

for(sp in uspecies){
    cat("doing :", sp, which(uspecies == sp), "/", length(uspecies), "\n")
    l_species_lines[[sp]] <- rec_mix_clim(path, files, sp, climate_sp,
                                      BATOTcomp,BATOTSP, BATOTNonSP, minmax = FALSE)
}
rec_data_lines <- do.call(rbind, l_species_lines)

rec_data %>%
    ggplot(aes(x = BA, y = mean, ymin = rec_min, ymax = rec_max, color = clim, fill = clim)) +
    geom_ribbon(alpha = 0.2, linetype = 0) +
    geom_line( linewidth=0.5) +
    facet_wrap(species ~ vital, scales = "free") +
    scale_colour_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
    scale_fill_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
    scale_alpha_discrete(range = c(0.08, 1), name = "mean") +
    NULL

rec_data_lines %>%
    # filter(value > 0) %>%
    mutate(alpha = models == "mean") %>%
    ggplot(aes(x = mesh, y = value, group = interaction(models,clim),
               color = clim,  alpha = alpha)) +
    geom_line( linewidth=0.5) +
    facet_wrap(species ~ vital, scales = "free") +
    scale_colour_manual(values = c("cadetblue3", "chartreuse3", "coral2")) +
    scale_alpha_discrete(range = c(0.08, 1), name = "mean") +
    NULL
