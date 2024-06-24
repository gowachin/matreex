# Compute mean and sd of species dist for young stand for Jasper

#FUNDIV <- read.csv("data/NewFundiv/FUNDIV_tree.csv")
FUNDIV <- read.csv("data/NewFundiv/FUNDIV_data.csv", sep = " ")

hist(FUNDIV$dbh2, breaks = 100)
abline(v = 100, col = "red")

FUNDIV <- FUNDIV[FUNDIV$dbh2 >100, , ]
FUNDIV <- FUNDIV[FUNDIV$treestatus %in% c(1,2), ]
 # 1 2 in growth alive
library(dplyr)
library(tidyr)
library(Hmisc)

# Compute dominant species in ba
fun_dom_sp <- function (sp,ba, Nha){
  res <- data.frame(sp = sp, ba_ha2 = ba, Nha = Nha) %>% group_by(sp) %>% summarise(BA = sum(ba_ha2), N = sum(Nha)) %>%
    arrange(desc(BA)) %>%ungroup() %>% mutate(BA_per = BA/sum(BA))
  return((res[1,]))
}

df <- FUNDIV %>% group_by(plotcode) %>%
  summarise(country = unique(country),
            MQD_log = sum(log(dbh2)*Nha2)/sum(Nha2),
            MQD_log_sd = sqrt(wtd.var(log(dbh2), Nha2)),
            domsp = list(fun_dom_sp(sp=species, ba = Nha2 *(dbh2/1000)^2/4*pi, Nha = Nha2)))%>%
  unnest_wider(domsp)

sp_sel <- names(sort(table(df$sp), decreasing = TRUE)[1:45])
c("Pinus sylvestris", "Picea abies", "Fagus sylvatica", "Quercus ilex",
  "Pinus pinaster", "Pinus halepensis", "Quercus robur", "Quercus petraea",
  "Pinus nigra", "Betula sp")
res <- df %>% filter(sp %in% sp_sel ,
                     BA_per >0.9,
                     MQD_log < log(120))%>%
  group_by(sp) %>% summarise(MQD_log = mean(MQD_log),
                             MQD_log_sd = mean(MQD_log_sd),
                             Nha = mean(N),
                             BA = mean(BA))

res <- recrut[, -1]
names(res) <- c("sp", "MQD_log", "MQD_log_sd", "Nha", "BAha")

# Plot of the simulated distribution
par(mfrow = c(2,5))
for (i in 1:10){
  hist(exp(rnorm(10000, mean = res$MQD_log[i], sd = res$MQD_log_sd[i])),
       main = res$sp[i], breaks = 40, probability= TRUE, ylim = c(0, 1))
  # size mesh for tree size
  mesh <- 100:150
  # compute relative probability density for each mesh based on lognormal proba
  prob_dens <- dnorm(log(mesh), mean = res$MQD_log[i], sd = res$MQD_log_sd[i])/
    (sum(dnorm(log(mesh), mean = res$MQD_log[i], sd = res$MQD_log_sd[i])))
  lines(mesh, prob_dens, col = "blue")
  abline(v = 99, col = "red")
  # compute basal area from size dist and mean density and compare to mean basal area in data
  print(res$sp[i])
  print(sum(mesh^2/(4*1000*1000) *pi *  prob_dens *res$Nha[i]))
  print(res$BA[i])
}
names(res) <- c("species", "mean_DBH_log", "sd_DBH_log", "Nha", "BAha")

write.csv(res, file = "ParamSizeDistYoungStand.csv")
#How to introduce this new cohort ? Start after the time lag of the species? As plantation shorter than the full time lag ? By how much ?10 years?
# A given percentage of teh time lag ? 50 % ?


#Not possible to get the structure in the lag Y years before by inverting the IPM because the IPM is non invertible (its determinant is zero)

# writing a function to use these parameters
matreex::distrib_planting
load_all()
data("fit_Picea_abies")
data("climate_species")
climate <- subset(climate_species, N == 2 & sp == "Picea_abies", select = -sp)

Picea_ipm <- make_IPM(
    species = "Picea_abies",
    climate = climate,
    fit = fit_Picea_abies,
    clim_lab = "optimum clim",
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Picea_abies) * 1.1),
    BA = 0:60, # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
)

fun <- def_init_planting("Picea_abies")
fun(Picea_ipm$mesh, 0.03)
