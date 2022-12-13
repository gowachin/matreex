new_disturb <- function(type = c("storm", "fire", "biotic", "snow", "other"),
                        intensity, duration = 1){

    disturb <- data.frame(type = type, intensity = instensity, duration = duration)

    return(disturb)
}
