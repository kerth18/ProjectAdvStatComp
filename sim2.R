library(gtools)
library(foreach)
library(doParallel)  
library(parallel)
library("mvtnorm")
library(xlsx)

numCores <- detectCores()  
cl <- makeCluster(numCores)  
registerDoParallel(cl)

#### OPERATIONS ON VECTORS ####

# Difference
D <- function(X,Y)
{
  if (identical(X,Y)==TRUE){d <- rep(0,length(X))
  }   else {d <- X-Y}  
  return(d)
}

# Akr
A <- function(X,Y)
{
  a <- crossprod(D(X,Y),D(X,Y))
  return(a)
}

# A2krls
Ad <- function(X,Y,U,V)
{
  a <- (crossprod(D(X,Y),D(U,V)))^2
  return(a)
}

#Bkr
B <- function(X,Y)
{
  b <- (crossprod(X,Y))^2
  return(b)
}

# D(k,r)
Dpar <- function(X,k,r)
{
  Y <- X[-c(k,r),]
  d <- X[k,] - colMeans(Y)
  return(d)
}

# C2kr
Cd <- function(X,k,r)
{
  c <- crossprod(X[k,],Dpar(X,k,r))*crossprod(X[r,],Dpar(X,r,k))
  return(c)
}

# B2krls
Bd <- function(X,Y,U,V)
{
  b <- Ad(X,Y,U,V) + Ad(X,U,Y,V) + Ad(X,V,U,Y)
  return(b)
}

Q <- function(ni)
{
  q <- ni*(ni-1)
  return(q)
}

P <- function(ni)
{
  p <- ni*(ni-1)*(ni-2)*(ni-3)
  return(p)
}

#### ESTIMATORS ####

E1 <- function(X)
{
  e1 <- sum(diag(cov(X)%*%cov(X)))
  return(e1)
}

E2 <- function(X)
{
  n <- dim(X)[1]
  e2 <- n^2/((n+2)*(n-1))*(sum(diag(cov(X)%*%cov(X))) - 1/n*(sum(diag(cov(X))))^2)
  return(e2)
}

E3 <- function(X,tperm2)
{
  n <- dim(X)[1]
  #perm <- t(permutations(n,2))
  e3 <- sum(foreach(i=tperm2, .combine='c', .export='B') %dopar% {B(X[i[1,],],X[i[2,],])})/Q(n)
  return(e3)
}

E4 <- function(X,tperm2)
{
  n <- dim(X)[1]
  #perm <- t(permutations(n,2))
  e4 <- sum(foreach(i=tperm2, .combine='c', .export=c('Cd','Dpar')) %dopar% {Cd(X,i[1,],i[2,])})/Q(n)
  return(e4)
}

E5 <- function(X,tperm4)
{
  n <- dim(X)[1]
  #perm <- t(permutations(n,4))
 e5 <- sum(foreach(i=tperm4, .combine='c', .export=c('Ad','D')) %dopar% {Ad(X[i[1,],],X[i[2,],],X[i[3,],],X[i[4,],])})/(4*P(n))
 return(e5)           
}

E6 <- function(X,tperm4)
{
  n <- dim(X)[1]
  #perm <- t(permutations(n,4))
  e6 <- sum(foreach(i=tperm4, .combine='c', .export=c('Ad','Bd','D')) %dopar% {Bd(X[i[1,],],X[i[2,],],X[i[3,],],X[i[4,],])})/(12*P(n))
  return(e6)
}

### Sigma = Identity Matrix ###
####### MVN distribution ######

Nsim <- 5000
n <- 10 #c(5,10,20,30,50)
p <- 300 #c(5,10,20,50,100,300,500,1000)
perm2 <- permutations(n,2)
perm4 <- permutations(n,4)
tperm2 <- t(perm2)
tperm4 <- t(perm4)

sim <- function(n, p, Nsim, tperm2, tperm4){
  A <- foreach(icount(Nsim), .combine='rbind', 
               .export=c('rmvnorm','permutations','D','B','A','Ad','Dpar','Cd','Bd','Q','P','E1','E2','E3','E4','E5','E6'), 
               .packages="foreach") %dopar% {
    X <- rmvnorm(n, mean = rep(0, p))
    c(E1(X),E2(X),E3(X,tperm2),E4(X,tperm2),E5(X,tperm4),E6(X,tperm4))
  }
  return(A)
  e <- colSums(A.pp)/Nsim/p
  v <- colSums((A.pp-e)^2)/Nsim/p^2
  df <- data.frame(rbind(e,v,A.pp), row.names=c("E","Var",1:Nsim))
  return(df)
}

# collect the results
res <- sim.pp(n,p,Nsim,perm2,perm4)
file_name <- paste('resultsI_',toString(n),toString(p),'.xls',sep='')
write.xlsx(res, file_name)


