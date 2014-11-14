#!/usr/bin/R
#
# A test area for DEA analysis, using various mock and prepared datasets
#
# Jeff Shrader
# First version: 2011-12-29
# Time-stamp: "2012-09-12 18:05:34 jgs"
#

rm(list = ls())
library(FEAR)
library(linprog)
library(lpSolveAPI)
library(rbenchmark)

## Basic Model from 'lpsolve' website section "More theoretical example"
## http://lpsolve.sourceforge.net/5.5/Python.htm

nsims <- 1000
nrep <- 1
## 1. linprog
t <- proc.time()
for(i in 1:nsims){
  f <- c(-4, -2, -1)
  A <- matrix(c(2, 1, 0,
                1, 0, 2,
                1, 1, 1,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1), nrow=3, ncol=9)
  e <- matrix(c('<', '<', '<', '<', '<', '<', '>', '>', '>'), nrow=1, ncol=9)
  b <- matrix(c(1, 2, 1,
                1, 1, 2,
                0, 0, 0), nrow=1, ncol=9)
}
proc.time() - t

### Pre-transposed constraints
t <- proc.time()
for(i in 1:nsims){
  opt <- lp('min', f, A, e, b, transpose.constraints=FALSE)
}
proc.time() - t
print(opt$objval)
rm(opt)

### Non-transposed constraints
tA <- t(A)
t <- proc.time()
for(i in 1:nsims){
  opt <- lp('min', f, tA, e, b)
}
proc.time() - t
print(opt$objval)
rm(opt)


### Using diagonal matrices for bounds
t <- proc.time()
for(i in 1:nsims){
  f <- c(4, 2, 1)
  a <- matrix(c(2, 1, 0,
                1, 0, 2,
                1, 1, 1), nrow=3, ncol=3)
  u <- diag(3)
  l <- diag(3)
  A <- cbind(a, u, l)
  e <- matrix(c('<', '<', '<', '<', '<', '<', '>', '>', '>'), nrow=1, ncol=9)
  b <- matrix(c(1, 2, 1,
                1, 1, 2,
                0, 0, 0), nrow=1, ncol=9)
}
proc.time() - t

t <- proc.time()
for(i in 1:nsims){
  opt <- lp('min', f, A, e, b, transpose.constraints=FALSE)
}
proc.time() - t
print(opt$objval)
rm(opt)


### Dropping non-negativity
t <- proc.time()
for(i in 1:nsims){
  f <- c(4, 2, 1)
  a <- matrix(c(2, 1, 0,
                1, 0, 2,
                1, 1, 1), nrow=3, ncol=3)
  u <- diag(3)
  A <- cbind(a, u)
  e <- matrix(c('<', '<', '<', '<', '<', '<'), nrow=1, ncol=6)
  b <- matrix(c(1, 2, 1,
                1, 1, 2), nrow=1, ncol=6)
}
proc.time() - t

for(j in 1:nrep){
  t <- proc.time()
  for(i in 1:nsims){
    opt <- lp('min', f, A, e, b, transpose.constraints=FALSE)
  }
  proc.time() - t
}
print(opt$objval)
rm(opt)


## 2. lp_solve API
t <- proc.time()

for(i in 1:nsims){
  lprec <- make.lp(0, 3)
  set.objfn(lprec, c(-4, -2, -1))
  add.constraint(lprec, c(2, 1, 0), '<', 1)
  add.constraint(lprec, c(1, 0, 2), '<', 2)
  add.constraint(lprec, c(1, 1, 1), '<', 1)
  set.bounds(lprec, lower = c(0, 0, 0), columns = c(1, 2, 3))
  set.bounds(lprec, upper = c(1, 1, 2), columns = c(1, 2, 3))
  solve(lprec)
}
proc.time() - t
get.objective(lprec)
get.variables(lprec)
rm(lprec)

## 3. Variations on set-up: Deleting lp
t <- proc.time()
for(i in 1:nsims){
  f <- c(4, 2, 1)
  A <- matrix(c(2, 1, 0,
                1, 0, 2,
                1, 1, 1,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1), nrow=3, ncol=9)
  e <- matrix(c('<', '<', '<', '<', '<', '<', '>', '>', '>'), nrow=1, ncol=9)
  b <- matrix(c(1, 2, 1,
                1, 1, 2,
                0, 0, 0), nrow=1, ncol=9)
  opt <- lp('max', f, A, e, b, transpose.constraints=FALSE)
  rm(opt)
}
proc.time() - t

## Large Model
nsims <- 3000
set.seed(40593)

sequence <- seq(100,nsims,5)
systime <- rep(0, length(sequence))
elapsetime <- systime

j <- 1
for(i in sequence){
  t <- proc.time()
  o <- i
  h <- 2
  f <- matrix(runif(o, 0, 50), ncol=o, nrow=1)
  A <- matrix(runif(o*h, 0, 50), ncol=h, nrow=o)
  b <- matrix(c(runif(h, 0, 50), rbinom(o, 3, .5), rep(0, o)), ncol=(h+(2*o)), nrow=1)
  u <- diag(o)
  l <- diag(o)
  A <- cbind(A, u, l)
  e <- c(rep('<',h+o), rep('>',o))
  opt <- lp('max', f, A, e, b, transpose.constraints=FALSE)
  tt <- proc.time() - t
  systime[j] <- tt[2]
  elapsetime[j] <- tt[3]
  j <- j+1
  rm(opt)
}
print(opt$objval)
time <- 1:length(sequence)
plot(time, elapsetime, type="l")
lines(time, systime, col="blue")
  

## Diagonalization
nsims <- 500

t <- proc.time()
for(i in 1:nsims){
  o <- 1000
  diag(o)
}
proc.time() - t

t <- proc.time()
for(i in 1:nsims){
  o <- 1000
  diag(x = 1, nrow=o, ncol=o)
}
proc.time() - t

t <- proc.time()
for(i in 1:nsims){
  o <- rep(1, 1000)
  diag(o)
}
proc.time() - t



## Transposition
nsims <- 2000
nrep <- 10
set.seed(40593)

sequence <- seq(100,nsims,50)

systimet <- rep(0, length(sequence)*nrep)
systimen <- systimet
elapsetimet <- systimet
elapsetimen <- systimet
yaxisn <- systimet

j <- 1
for(i in sequence){
  for(k in 1:nrep){
    o <- i
    h <- 10
    f <- matrix(runif(o, 0, 50), ncol=o, nrow=1)
    A <- matrix(runif(o*h, 0, 50), ncol=h, nrow=o)
    b <- matrix(c(runif(h, 0, 50), rbinom(o, 3, .5),
                  rep(0, o)), ncol=(h+(2*o)), nrow=1)
    u <- diag(o)
    l <- diag(o)
    A <- cbind(A, u, l)
    At <- t(A)
    e <- c(rep('<',h+o), rep('>',o))
    t <- proc.time()
    opt <- lp('max', f, A, e, b, transpose.constraints=FALSE)
    tt <- proc.time() - t
    rm(opt)  
    ttt <- proc.time()
    opt <- lp('max', f, At, e, b)
    tttt <- proc.time() - ttt
    rm(opt)
    systimet[j] <- tt[2]
    elapsetimet[j] <- tt[3]
    systimen[j] <- tttt[2]
    elapsetimen[j] <- tttt[3]
    yaxisn[j] <- i
    j <- j+1  
  }
}

plot(yaxisn, elapsetimen)
points(yaxisn, elapsetimet, col='blue')

plot(yaxisn, systimen)
points(yaxisn, systimet, col="blue")

data <- data.frame(systimet, elapsetimet, systimen, elapsetimen, yaxisn)
test <- aggregate(data, by=list(yaxisn), FUN="mean")

plot(test$yaxisn, test$systimen, type='l')
lines(test$yaxisn, test$systimet, col="blue")
plot(test$yaxisn, test$elapsetimen, type='l')
lines(test$yaxisn, test$elapsetimet, col="blue")


## Non-negativity
rm(list = ls())
library(linprog)
library(lpSolveAPI)
maxmat <- 3000
stepsize <- 50
nsims <- round(maxmat/stepsize)
nrep <- 10
set.seed(40593)

sequence <- seq(100,nsims*stepsize,stepsize)

systimet <- rep(0, length(sequence)*nrep)
systimen <- systimet
elapsetimet <- systimet
elapsetimen <- systimet
yaxisn <- systimet
optt <- systimet
optn <- systimen

j <- 1
for(i in sequence){
  for(k in 1:nrep){
    o <- i
    h <- 10
    f <- matrix(runif(o, 0, 50), ncol=o, nrow=1)
    A <- matrix(runif(o*h, -5, 50), ncol=h, nrow=o)
    bl <- matrix(c(runif(h, 0, 50), rbinom(o, 3, .5)), ncol=(h+o), nrow=1)
    b <- c(bl, rep(0, o))
    u <- diag(o)
    l <- diag(o)
    Al <- cbind(A, u)    
    A <- cbind(A, u, l)
    e <- c(rep('<',h+o), rep('>',o))
    el <- c(rep('<',h), rep('>',o))    
    t <- proc.time()
    opt <- lp('max', f, A, e, b, transpose.constraints=FALSE)
    tt <- proc.time() - t
    optt[j] <- opt$solution[1]
    rm(opt)  
    ttt <- proc.time()
    opt <- lp('max', f, Al, el, bl, transpose.constraints=FALSE)
    tttt <- proc.time() - ttt
    optn[j] <- opt$solution[1]    
    rm(opt)
    systimet[j] <- tt[2]
    elapsetimet[j] <- tt[3]
    systimen[j] <- tttt[2]
    elapsetimen[j] <- tttt[3]
    yaxisn[j] <- i
    j <- j+1  
  }
}

plot(yaxisn, elapsetimen)
points(yaxisn, elapsetimet, col='blue')

plot(yaxisn, systimen)
points(yaxisn, systimet, col="blue")

data <- data.frame(systimet, elapsetimet, systimen, elapsetimen, yaxisn)
test <- aggregate(data, by=list(yaxisn), FUN="mean")

plot(test$yaxisn, test$systimen, type='l')
lines(test$yaxisn, test$systimet, col="blue")
plot(test$yaxisn, test$elapsetimen, type='l')
lines(test$yaxisn, test$elapsetimet, col="blue")



## Binary vars
rm(list = ls())
library(linprog)
library(lpSolveAPI)
maxmat <- 1000
stepsize <- 50
nsims <- round(maxmat/stepsize)
nrep <- 10
set.seed(40593)

sequence <- seq(100,nsims*stepsize,stepsize)

systimet <- rep(0, length(sequence)*nrep)
systimen <- systimet
elapsetimet <- systimet
elapsetimen <- systimet
yaxisn <- systimet
optt <- systimet
optn <- systimen

j <- 1
for(i in sequence){
  for(k in 1:nrep){
    theta <- 1
    o <- i
    h <- 10
    f <- matrix(c(1, rep(0, o)), ncol=(o+1), nrow=1) #runif(o, -5, 5)), ncol=(o+1), nrow=1)
    c <- t(matrix(c(0, rep(1, o)), ncol=(o+1), nrow=1))
    A <- matrix(runif((o+1)*h, 0, 25), ncol=h, nrow=(o+1))
    A <- cbind(A, c)
    b <- matrix(c(runif(h, 0, 50), 1), ncol=(h+1), nrow=1)
    bb <- c(b, rep(1, o), rep(0, o))
    u <- rbind(rep(0, o), diag(o))
    l <- rbind(rep(0, o), diag(o))
    Ab <- cbind(A, u, l)
    eb <- c(rep('<',h), '=', rep('<',o), rep('>',o))
    e <- c(rep('<',h), '=')    
    t <- proc.time()
    opt <- lp('max', f, Ab, eb, bb, transpose.constraints=FALSE)
    tt <- proc.time() - t
    optt[j] <- opt$solution[1]
    rm(opt)  
    ttt <- proc.time()
    opt <- lp('max', f, A, e, b, binary.vec=2:(o+1), transpose.constraints=FALSE)
    tttt <- proc.time() - ttt
    optn[j] <- opt$solution[1]    
    rm(opt)
    systimet[j] <- tt[2]
    elapsetimet[j] <- tt[3]
    systimen[j] <- tttt[2]
    elapsetimen[j] <- tttt[3]
    yaxisn[j] <- i
    j <- j+1  
  }
}

plot(yaxisn, elapsetimen)
points(yaxisn, elapsetimet, col='blue')

plot(yaxisn, systimen)
points(yaxisn, systimet, col="blue")

data <- data.frame(systimet, elapsetimet, systimen, elapsetimen, yaxisn)
test <- aggregate(data, by=list(yaxisn), FUN="mean")

plot(test$yaxisn, test$systimen, type='l')
lines(test$yaxisn, test$systimet, col="blue")
plot(test$yaxisn, test$elapsetimen, type='l')
lines(test$yaxisn, test$elapsetimet, col="blue")


## Parallelization
rm(list = ls())
library(linprog)

SimLP <- function(ns, matsize, par='Y'){
  set.seed(40593)  
  ldir <- list()
  lf <- list()
  lA <- list()
  le <- list()
  lb <- list()
  dir <- 'max'
  for(i in 1:ns){
    t <- proc.time()
    o <- matsize
    h <- 20
    f <- matrix(runif(o, 0, 50), ncol=o, nrow=1)
    A <- matrix(runif(o*h, 0, 50), ncol=h, nrow=o)
    b <- matrix(c(runif(h, 0, 50), rbinom(o, 3, .5), rep(0, o)), ncol=(h+(2*o)), nrow=1)
    u <- diag(o)
    l <- diag(o)
    A <- cbind(A, u, l)
    e <- c(rep('<',h+o), rep('>',o))
    ldir[i] <- list(dir=dir)
    lf[i] <- list(f=f)
    lA[i] <- list(A=A)
    le[i] <- list(e=e)
    lb[i] <- list(b=b)
  }
  
  if(par=='M'){
    opt <- foreach(i=1:ns, .combine=c, .packages='linprog') %dopar%
    {
      lp(ldir[[i]], lf[[i]], lA[[i]], le[[i]], lb[[i]], transpose.constraints=FALSE)$solution
      ##temp.s
    }
  }else if(par=='V'){
    opt <- mapply(lp, ldir, lf, lA, le, lb, transpose.constraints=FALSE)['solution',]
  }else if(par=='L'){
    opt <- list()
    for(i in 1:ns){
      opt[[i]] <- lp(ldir[[i]], lf[[i]], lA[[i]], le[[i]], lb[[i]], transpose.constraints=FALSE)$solution
    }
  }
  return(opt)
}

## Test the differnt methods
nsims <- 1000
matsize <- 10

require(foreach)
require(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

system.time(out1 <- SimLP(nsims, matsize, par='M'))
stopCluster(cl)

system.time(out2 <- SimLP(nsims, matsize, par='V'))
system.time(out3 <- SimLP(nsims, matsize, par='L'))


## Other parallelization stuff
rm(list = ls())
library(parallel)
runs <- 1e6
manyruns <- function(n) mean(unlist(lapply(X=1:(runs/4), FUN=onerun)))
 

cores <- 4
cl <- makeCluster(cores)
 
# Send function to workers
tobeignored <- clusterEvalQ(cl, {
    onerun <- function(.){ # Function of no arguments
        doors <- 1:3
        prize.door <- sample(doors, size=1)
        choice <- sample(doors, size=1)
        if (choice==prize.door) return(0) else return(1) # Always switch
    }
; NULL
})
 
# Send runs to the workers
tobeignored <- clusterEvalQ(cl, {runs <- 1e6; NULL})
runtime <- system.time({
    avg <- mean(unlist(clusterApply(cl=cl, x=rep(runs, 4), fun=manyruns)))
})[3]
stopCluster(cl)
 
cbind(avg, runtime)


## Another example
require(foreach)
require(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

x=iris[which(iris[,5] != "setosa"),c(1,5)]

trials = 1000 
system.time( 
  out <- foreach(icount(trials), .combine=cbind) %do% 
  {  
    ind=sample(100,100,replace=TRUE) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    coefficients(results1) 
  })[3] 

system.time( 
  out <- foreach(icount(trials), .combine=cbind) %dopar% 
  {  
    ind=sample(100,100,replace=TRUE) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    coefficients(results1) 
  })[3] 

out <- matrix(0, nrow=2, ncol=1000)
system.time( 
  for(i in 1:trials){
    ind=sample(100,100,replace=TRUE) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    out[,i] <- coefficients(results1)
  })[3] 

out <- matrix(0, nrow=2, ncol=1000)
system.time( 
  test<-foreach(i=1:trials) %dopar% 
  {  
    ind=sample(100,100,replace=TRUE) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    results1 = glm(x[ind,2]~x[ind,1],family=binomial(logit)) 
    ##out[,i] <- coefficients(results1)
    coefficients(results1)
  })[3] 

stopCluster(cl)

