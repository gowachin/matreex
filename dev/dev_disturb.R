
disturb_fun <- function(x, species, disturb = NULL, ...){
    size <- species$IPM$mesh
    coef <- species$disturb_coef
    if(any(disturb$type %in% coef$disturbance)){
        coef <- subset(coef, disturbance == disturb$type)
    } else {
        stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                     sp_name(species), disturb$type))
    }


    qmd <- treeforce::QMD(size = size, n = x)
    logratio <-  log(size / qmd)
    dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
    logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
    Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled +
                        coef$b * disturb$intensity^(coef$c * dbh.scaled)) ^ disturb$duration

    return(x* Pkill) # always return the mortality distribution
}

Picea_abies$disturb_fun <- disturb_fun

# plot(equil, type= "l", xlab = "size index", ylab = "n")
# lines(1:700, equil - disturb_fun(equil, Picea_abies, ex_disturb), col = "red")
# legend("topright", c("Distribution", "After Disturbance"),
#        lty = c(1, 1), col = c(1, 2))

time <- 2500
disturb <- data.frame(dist =  FALSE, type = "none", intensity = 0,
                      duration = 1, IsSurv = TRUE, t = 1:time)

disturb[100, ] <- data.frame(dist = TRUE, type = "storm", intensity = 0.2,
                             duration = 1, IsSurv = FALSE, t = 100)
# disturb[100, ] <- data.frame(dist = TRUE, type = "storm", intensity = 0.5,
#                              duration = 1, IsSurv = FALSE, t = 100)
disturb[98:102, ]

load_all()
Picea_abies$init_pop <- def_init_k(equil * 0.03)
forest_ipm <- new_forest(species = list(Picea = Picea_abies))
set.seed(42)
memor <- sim_deter_forest.forest(forest_ipm, tlim = time,
                                 equil_dist = time, equil_time = time,
                                 disturbance  = disturb,
                                 verbose = TRUE, correction = "cut") %>%
    tree_format()

memor %>%
    filter(var %in% c("BAsp", "H", "N"), ! equil, value != 0) %>%
    ggplot(aes(x = time, y = value)) +
    facet_wrap(~ var, scales = "free_y") +
    geom_line(size = .4) + geom_point(size = .4) +
    NULL
