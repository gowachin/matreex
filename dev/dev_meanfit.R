# checking meanfit

load_all()
library(purrr)

# l <- lapply(1:100, function(x) old_fit2fit("Yggdrasil", replicat = x)) # 15.41s
f_fit <- here(here::here(), "output", "Yggdrasil", "fit_sgr_all.Rds")
l <- readRDS(f_fit) # 0.16 sec

clim <- readRDS("output/Yggdrasil/plots_pred.Rds")
climate <- c(unlist(clim[2,]),  BATOTcomp = 20)

size <- seq(90,1500, by = 2)

mean <- old_fit2fit("Yggdrasil", mean = TRUE)

mean_gr <- exp_sizeFun(mean$gr$params_m, climate)
mean_sv <- exp_sizeFun(mean$sv$params_m, climate)

inv <- l[[1]]$sv$family$linkinv

l_gr <- map(l, ~ exp_sizeFun(.x$gr$params_m, climate)) %>%
    map(~ exec(.x, size), size = size)
l_sv <- map(l, ~ exp_sizeFun(.x$sv$params_m, climate)) %>%
    map(~ 1 - inv(exec(.x, size)), size = size)

col <- map_dbl(l, ~ sum(grep("sgddb", names(.x$sv$params_m))) > 0) + 1

col
# range(unlist(map(l_sv, range)))
par(mfrow = c(1, 2))
plot(1, type = "n", ylim = c(.9, 1), xlim = c(90, 1500),
     xlab = "size", ylab = "", main = "Survival")
for(i in seq_along(l_sv)) {
    lines(size, l_sv[[i]], col = "gray") #col[i])
    }
lines(size, 1 - inv(mean_sv(size)))


plot(1, type = "n", ylim = range(unlist(map(l_gr, range))), xlim = c(90, 1500),
     xlab = "size", ylab = "", main = "Growth")
for(i in seq_along(l_gr)) {
    lines(size, l_gr[[i]], col = "gray")
}
lines(size, mean_gr(size))
par(mfrow = c(1, 1))
