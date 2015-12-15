setwd("C:/Users/estbo560/Documents/COURSES/Advanced Statistical Computing/Projet")

### libraries
library("devtools")
library(microbenchmark)
library("Rcpp")
library("mvtnorm")
library(gtools)
library(foreach)
library(doParallel)  
library(parallel)
library(xlsx)

# sourceCpp("operations.cpp")
Rcpp.package.skeleton("MyPackage", example_code = FALSE, attributes = TRUE,
                      cpp_files = c("operations.cpp"))
Rcpp::compileAttributes('./MyPackage')

install("MyPackage")
library("MyPackage")

numCores <- detectCores()  
cl <- makeCluster(numCores)  
registerDoParallel(cl)

### Comparisons of rmvnorm
microbenchmark(rmvnorm(5, mean=rep(0,500)),
               rmvnorm(5, mean=rep(0,1000)),
               rmvnorm(20, mean=rep(0,500)),
               rmvnorm(20, mean=rep(0,1000)),
               rmvnorm(50, mean=rep(0,500)),
               rmvnorm(50, mean=rep(0,1000)))

### Comparisons of permutations(n,k)
microbenchmark(permutations(20,4),
               permutations(30,4),
               permutations(50,4))


### Comparisons of operations on vectors/matrices
source("sim2.R")
x <- rnorm(10000)
y <- rnorm(10000)
u <- rnorm(10000)
v <- rnorm(10000)

n <- 10; p <- 500
X <- rmvnorm(n, mean = rep(0, p))
perm2 <- permutations(n,2)
perm4 <- permutations(n,4)
tperm2 <- t(perm2)
tperm4 <- t(perm4)

microbenchmark(App(x,y), A(x,y))
#etc.


### comparisons of functions E_i
n <- 10; p <- 500
X <- rmvnorm(n, mean = rep(0, p))
perm2 <- permutations(n,2)
perm4 <- permutations(n,4)
tperm2 <- t(perm2)
tperm4 <- t(perm4)
microbenchmark(E3(X,tperm2), E3pp(X,perm2))
microbenchmark(E4(X,tperm2), E4pp(X,perm2))
microbenchmark(E5(X,tperm4), E5pp(X,perm4))
microbenchmark(E6(X,tperm4), E6pp(X,perm4))


### Simulations

n <- 5
p <- 50
perm2 <- permutations(n,2)
perm4 <- permutations(n,4)
tperm2 <- t(perm2)
tperm4 <- t(perm4)

sim.pp <- function(n, p, Nsim, perm2, perm4){
  A.pp <- foreach(icount(Nsim), 
          .combine='rbind', 
          .export=c('rmvnorm',"E1", "E2"), 
          .packages=c("foreach", "MyPackage", "Rcpp")) %dopar% {
            X <- rmvnorm(n, mean = rep(0, p))
            c(E1(X),E2(X),E3pp(X,perm2),E4pp(X,perm2),E5pp(X,perm4),E6pp(X,perm4))
          }
  e <- colSums(A.pp)/Nsim/p
  v <- colSums((A.pp-e)^2)/Nsim/p^2
  df <- data.frame(rbind(e,v,A.pp), row.names=c("E","Var",1:Nsim))
  return(df)
}

# collect the results
res <- sim.pp(n,p,Nsim,perm2,perm4)
file_name <- paste('resultsI_',toString(n),toString(p),'.xls',sep='')
write.xlsx(res, file_name)

### Comparison simulations using sim.pp and sim
n <- 5; p <- 50
perm2 <- permutations(n,2)
perm4 <- permutations(n,4)
tperm2 <- t(perm2)
tperm4 <- t(perm4)
microbenchmark(sim.pp(n,p,100,perm2,perm4), sim(n,p,100,tperm2,tperm4))
microbenchmark(sim.pp(n,p,5000,perm2,perm4), sim(n,p,5000,tperm2,tperm4))

