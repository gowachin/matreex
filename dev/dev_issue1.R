# Working on multiple species systems ####

document()
load_all()
spe <- "Yggdrasil"
# Yggdrasil <- old_ipm2species(spe, path = "tests/testthat/testdata/", replicat = 1)
Yggdrasil <- old_ipm2species(spe, path = here(), replicat = 1)
Ents <- Yggdrasil
Ents$info["species"] <- "Ents"

Forest <- forest(list(Yggdrasil))

# Forest$species$Yggdrasil$recruit_fun

set.seed(42)
res <- sim_deter_forest(Forest, tlim = 600, equil_time = 1e3,
                        correction = "cut",
                        verbose = TRUE)

times <- as.numeric(sub("t", "", colnames(res)))
plot(times, res[grepl("BA",rownames(res)),], ylab = "Total BA", xlab = "time",
     cex = 0.1, type = "b", ylim = c(0, 100))
text(700, res[grepl("BA",rownames(res)),ncol(res)],
     round(res[grepl("BA",rownames(res)),ncol(res)], 2), cex = .7, pos = 3)

set.seed(42)
Forest2 <- forest(list(Yggdrasil, Ents))
res2 <- sim_deter_forest(Forest2, tlim = 600, equil_time = 1e3,
                         correction = "cut",
                         verbose = TRUE)


value <- colSums(res2[grepl("BA",rownames(res2)),])
points(times, value, type = "b", cex = 0.1, col = "darkred")
text(700, value[ncol(res2)], round(value[ncol(res2)], 2), cex = .7, pos = 3)
points(times, res2[grepl("Ygg.*BA",rownames(res2)),],
       type = "b", cex = 0.1, col = "darkblue")
points(times, res2[grepl("Ent.*BA",rownames(res2)),],
       type = "b", cex = 0.1, col = "darkgreen")
legend("top", c("Yddgrasil solo", "Yggdrasil + Ents", "Yggdrasil", "Ents"),
       fill = c("black", "darkred", "darkblue", "darkgreen"))
text(700, res2[grepl("BA",rownames(res2)),ncol(res2)],
     round(res2[grepl("BA",rownames(res2)),ncol(res2)], 2), cex = .7, pos = 3)

#' La difference de BA s'expliquerais par une plus faible competition sur
#' les cellules de mesh petits car la competition extrasp s'applique
#' lors du recrutement.

eq_1 <- res[grepl("m",rownames(res)),ncol(res)]
eq_2_Ygg <- res2[grepl("Ygg.*m",rownames(res2)),ncol(res2)]
eq_2_Ent <- res2[grepl("Ent.*m",rownames(res2)),ncol(res2)]
eq_2 <- eq_2_Ent + eq_2_Ygg

plot(eq_2- - eq_1, type = "b", cex = .2,
     xlab = "mesh cell", ylab = "Net difference",
     main = paste("Difference in equilibrium distribution between",
                  "\nintra-specific and hetero-specific competition"))

