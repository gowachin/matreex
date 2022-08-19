# Issue 1 : recrutment function comes with env ####

#' In recruitmenent function, the function is written and get values from
#' a binded environment coming from function definition.
load_all()
library(rlang)

species <- "Yggdrasil"
climatic <- 1
path <- here()
replicat <- 42

fIPM <- here(path, "output", species, paste0("IPM_Clim_", climatic, ".Rds"))
raw_IPM <- readRDS(assertFileExists(fIPM)) # NOTE 10" to load...
assertNumber(replicat, lower = 1, upper = length(raw_IPM))
raw_IPM <- raw_IPM[[replicat]]


func <- raw_IPM$RecFun
raw_IPM$RecFun
raw_IPM$rec$formula
raw_IPM$rec$params_m
raw_IPM$list_m

str(lobstr::ast(!!eval(parse(text = raw_IPM$rec$formula))))

e <- body(func)
ts <- e[[2]]
lobstr::ast(!!ts)
lobstr::sxp(ts)
rlang::parse_expr(!!ts)
parse(text  = parse(text = ts)[[2]])[[1]]

accumulate(ts, ~ parse(text = .x), .dir = "backward")

foo <- function(x) {
    browser()
    raw_IPM$RecFun(x)
}

foo(1)

environment(raw_IPM$RecFun)$K_i

# function (BATOTSP)
# {
#     as.numeric(K_i + K_BA * BATOTSP + K_logBA * log(BATOTSP))
# }
# <environment: 0x5616bc88bc18> # WHY IS THERE A BINDED env HERE ???

raw_IPM$rec$params_m
raw_IPM$list_m

# Browse[3]> K_i
# intercept
# -0.7360761
# Browse[3]> K_BA
# BATOTSP
# -0.0179456
# Browse[3]> K_logBA
# logBATOTSP
# -0.2630412

# Browse[3]> ls(envir = where("K_i"))
# [1] "IsFunc"               "K_BA"
# [3] "K_i"                  "K_i_inter"
# [5] "K_logBA"              "L"
# [7] "L1"                   "L2"
# [9] "list_covs"            "mu_rec"
# [11] "params_i"             "params_i_inter"
# [13] "params_i_no"          "params_logsize_inter"
# [15] "params_rec"           "params_size_inter"
# [17] "r_rec"

# Finding a way to build function from reg ####

#' Is it possible to write numeric values as is in a function ??
#' maybe with expr

x <- 2
library(rlang)
foo <- expr(x + 2)
foo
eval_bare(foo, x = 3)

x <- 1
y <- 2
z <- call2("<-", expr(x), y + 1)
z
eval(z)
x
y

call2("{",
      call2(
          "+", expr(x), {1}
      )
)

z <- 3
tot <- function(x, y){
    tmp <- z * y
    tmp <- x +  tmp
    tmp
    }
tot
tot(1, y = 0.5) # 2.5
formals(tot)
names(formals(tot))[2] <- "Y" # Change argument name !
formals(tot)
body(tot)[[4]] <- call2("<-", expr(tmp), call2("+", 1, expr(tmp)))
body(tot)[[5]] <- call2("return", expr(tmp))
body(tot)
body(tot) <- call2("{",
                   call2( "+", expr(x), call2("*", z, expr(Y))
                   ) )
tot
tot(1, Y = 0.5) # 2.5
# tot(1, y = 0.5)

test <- call2(
    "<-", expr(res), expr(function(x){1})
)
eval(test)
res
res()

function(x) {1}


lobstr::ast(foo <- function(x){x +1})

# it works !

# Final test ####

list_covs <- raw_IPM$list_m
params <- raw_IPM$rec$params_m

foo <- exp_recFun(params = raw_IPM$rec$params_m, list_covs = raw_IPM$list_m)
foo
foo(1, 2, 1:5, 0.03)

