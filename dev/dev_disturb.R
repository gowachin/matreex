
disturb_fun <- function(x, species, disturb = NULL, ...){

    dots <- list(...)
    qmd <- dots$qmd

    size <- species$IPM$mesh
    coef <- species$disturb_coef
    if(any(disturb$type %in% coef$disturbance)){
        coef <- subset(coef, disturbance == disturb$type)
    } else {
        stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                     sp_name(species), disturb$type))
    }

    # qmd <- matreex::QMD(size = size, n = x)
    logratio <-  log(size / qmd)
    dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
    logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
    Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                        coef$b * disturb$intensity^(coef$c * dbh.scaled))

    return(x* Pkill) # always return the mortality distribution
}

Picea_abies$disturb_fun <- disturb_fun

# plot(equil, type= "l", xlab = "size index", ylab = "n")
# lines(1:700, equil - disturb_fun(equil, Picea_abies, ex_disturb), col = "red")
# legend("topright", c("Distribution", "After Disturbance"),
#        lty = c(1, 1), col = c(1, 2))

time <- 12500
disturb <- data.frame(type = "storm", intensity = c(0.2, 0.4, 0.6, 0.8, 1),
                      IsSurv = FALSE, t = c(100, 2500, 5000, 7500, 10000))

# time <- 3000
# disturb <- data.frame(type = "storm", intensity = c(0, 0.4),
                      # IsSurv = FALSE, t = c(100, 2500))


load_all()
Picea_abies$init_pop <- def_init_k(equil * 0.03)
forest_ipm <- new_forest(species = list(Picea = Picea_abies))
set.seed(42)
memor <- sim_deter_forest.forest(forest_ipm, tlim = time,
                                 equil_dist = 250, equil_time = time,
                                 disturbance  = disturb,
                                 verbose = TRUE, correction = "cut") %>%
    tree_format()

memor %>%
    filter(var %in% c("BAsp"), ! equil, value != 0) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .4) + geom_point(size = .4) +
    NULL
