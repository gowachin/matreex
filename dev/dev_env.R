
rm(list = ls())
base <- new.env()
base$X <- 1:5
base$nsp <- 0



foo <- function(){
    res <- new.env()

    res$X <- 1:5
    res$nsp <- 0

    return(res)
}

a <- foo()
b <- foo()

attach(base)
X
ls()

detach(base)
X


with(a,{
    X <- X +1
})
a$X <- a$X - 1

a$X

b$X

fuu <- function(x){
    attach(x)
}

l <- list(A = a
          B = b)
