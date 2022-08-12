document()
load_all()
spe <- "Yggdrasil"
Yggdrasil <- old_ipm2species(spe, path = "tests/testthat/testdata/", replicat = 1)
o_Yggdrasil <- Yggdrasil
o_Yggdrasil$info["species"] <- "Pinea"

Forest <- forest(list(Yggdrasil, o_Yggdrasil))
Forest <- forest(list(Yggdrasil))

set.seed(42)
res <- sim_deter_forest(Forest, tlim = 10, equil_time = 1e3,
                 correction = "cut",
                 verbose = TRUE)

res[c(1:2, 700:702, 703,1402:1403, 1404:1405, 2103:2106, 2805:2806),
         c(1:3,30:31)]
res[c(1:2, 700:702, 703,1402:1403), c(1:3,30:31)]

plot(res[grepl("BA",rownames(res)),])
plot(res[grepl("N",rownames(res)),])

