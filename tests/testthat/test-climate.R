test_that("values2sqb works", {
  expect_equal(value2sqb(x = data.frame(sgdd = 2)), 
               data.frame(sgdd = 2, sgdd2 = 4, sgddb = 0.5))
  
  expect_equal(value2sqb(x = data.frame(sgdd = 3), inv_null = TRUE), 
               data.frame(sgdd = 3, sgdd2 = 9, sgddb = 0.25))
})

test_that("expand_clim works", {
  expect_equal(expand_clim(climate = data.frame(sgdd = 2, wai = 3), 
                           inv_null = c(sgdd = FALSE, wai = TRUE)), 
               data.frame(sgdd = 2, sgdd2 = 4, sgddb = 0.5, 
                          wai = 3, wai2 = 9, waib = 0.25))
})

test_that("clim_gradient works", {

  data("climate_species")
  climate <- subset(climate_species, sp == "Picea_abies", select = -sp)

  expect_equal(clim_gradient(
                climate,
                start_clim = "opti", end_clim = "opti", times = 2,
                noise = rnorm, sigma = c(1, 0.001), seed = 42,
                extrapolate = 0), 
              data.frame(sgdd = c(1446.0375866043, 1444.10192998576), 
                         sgdd2 = c(2091024.7018724, 2085430.3841886), 
                         sgddb = c(0.000691544956551425, 0.000692471894978951), 
                         wai = c(0.452301820781064, 0.452571554974688), 
                         wai2 = c(0.204576937081866, 0.204821012372207), 
                         waib = c(0.688562105817776, 0.688434243790094))
  )
  
  expect_equal(clim_gradient(
                climate,
                start_clim = "opti", end_clim = "opti", times = 2,
                noise = rnorm, sigma = c(1, 0.001), seed = 42,
                extrapolate = 50), 
              data.frame(sgdd = c(1446.0375866043, 1444.10192998576, 1445.02975656849), 
                         sgdd2 = c(2091024.7018724, 2085430.3841886, 2088110.9973684), 
                         sgddb = c(0.000691544956551425, 0.000692471894978951, 0.000692027271725321), 
                         wai = c(0.452571554974688, 0.452342960692868, 0.451832567853635), 
                         wai2 = c(0.204821012372207, 0.204614154088389, 0.20415266937321), 
                         waib = c(0.688434243790094, 0.688542601206901, 0.688784658880041))
  ) 
})
