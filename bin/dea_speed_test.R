#!/usr/bin/R
#
# A test area for DEA analysis, using various mock and prepared datasets
#
# Jeff Shrader
# First version: 2011-12-29
# Time-stamp: "2014-11-14 15:34:31 jeff.shrader"
#

## Testing each of the basic functions
rm(list = ls())
library(FEAR)
source("dea_package.R")
nsims <- 20
set.seed(2982)
J <- 100
nv <- 1
nf <- 2
nu <- 2
x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
l <- length(y.m)

timer <- matrix(0, nrow=nsims, ncol=30)
for(i in 1:nsims){
  timer[i, 1:3] <- proc.time()[1:3]

  timer[i, 4:6] <- proc.time()[1:3]
  testv <- f.dea(f.inputs = x.m.f, outputs=y.m, tech="V", orientation="OUT", report.z="YES", slack="YES")
  timer[i, 7:9] <- proc.time()[1:3]
  testv <- f.dea(f.inputs = x.m.f, outputs=y.m, tech="V", orientation="OUT", report.z="YES", slack="NO")
  timer[i, 10:12] <- proc.time()[1:3]
  testf <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  timer[i, 13:15] <- proc.time()[1:3]
  testnc <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                  report.z="YES", slack="NO", convex="NO")

  timer[i, 16:18] <- proc.time()[1:3]
  act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)
  outf <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testf$theta, z=testf$z)
  outv <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testv$theta, z=testv$z)
  outnc <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testnc$theta, z=testnc$z)

  timer[i, 19:21] <- proc.time()[1:3]
  jim.c <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                 f.in.opt=outf$x.f.s, v.in.opt=outf$x.v.s, out.opt=outf$u.s,
                 convex='YES')
  timer[i, 22:24] <- proc.time()[1:3]
  jim <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
               f.in.opt=outf$x.f.s, v.in.opt=outf$x.v.s, out.opt=outf$u.s,
               convex='NO')
  timer[i, 25:27] <- proc.time()[1:3]
  jim.nc <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                  f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
                  convex='NO')

  timer[i, 28:30] <- proc.time()[1:3]
}

timer <- data.frame(timer)
names(timer) <- c('1','2','3','4','5','6',
                  'DEA_slack','x1','x2','DEA_no_slack','x3','x4','DEA_convex','x5','x6',
                  'DEA_non_convex','x7','x8','Aggregation','x9','x10','JIM_convex','x11','x12',
                  'JIM_non_convex_1','x13','x14','JIM_non_convex_2','x15','x16')
timer.out <- timer
for(j in seq(7,30,3)){
  timer.out[,j:(j+2)] <- timer[,j:(j+2)] - timer[,(j-3):(j-1)]
}
colMeans(timer.out[,7:30])



## Setup for testing differences in non-convex JIM
rm(list = ls())
source("dea_package.R")
nsims <- 20
J <- 50
nv <- 1
nf <- 2
nu <- 2

set.seed(2982)
jim.1 <- rep(0, nsims)
veri.1 <- jim.1
system.time(for(i in 1:nsims){
  x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
  x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
  x.m.f[,1] <- x.m.f[,1]*x.m.v[,1]
  y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
  veri.1[i] <- mean(y.m)
  l <- length(y.m)

  testnc <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                  report.z="YES", slack="NO", convex="NO")
  act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)
  outnc <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testnc$theta, z=testnc$z)
  jim.nc <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                 f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
                 convex='NO')
  jim.1[i] <- jim.nc$theta
})

set.seed(2982)
jim.2 <- rep(0, nsims)
veri.2 <- jim.2
system.time(for(i in 1:nsims){
  x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
  x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
  x.m.f <- x.m.f*x.m.v[,1]
  y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
  veri.2[i] <- mean(y.m)
  l <- length(y.m)

  testnc <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                  report.z="YES", slack="NO", convex="NO")
  act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)
  outnc <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testnc$theta, z=testnc$z)
  jim.nc <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                 f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
                 convex='NO')
  jim.2[i] <- jim.nc$theta
})

set.seed(2982)
jim.3 <- rep(0, nsims)
veri.3 <- jim.3
system.time(for(i in 1:nsims){
  x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
  x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
  x.m.f <- x.m.f
  y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
  veri.3[i] <- mean(y.m)
  l <- length(y.m)

  testnc <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                  report.z="YES", slack="NO", convex="NO")
  act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)
  outnc <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testnc$theta, z=testnc$z)
  jim.nc <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                 f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
                 convex='NO')
  jim.3[i] <- jim.nc$theta
})



system.time(for(i in 1:nsims){
  jim.nc <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                  f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
                  convex='NO')
})/nsims


require(foreach)
require(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)

system.time( 
  jim.nc<-foreach(i=1:nsims, .combine=c, .packages='linprog') %dopar% 
  {  
  f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
         f.in.opt=outnc$x.f.s, v.in.opt=outnc$x.v.s, out.opt=outnc$u.s,
         convex='NO')
})/nsims




## Testing option effects
rm(list = ls())
library(FEAR)
source("dea_package.R")

J <- 100
nv <- 10
nf <- 2
nu <- 2
x.m.f <- t(matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J))
x.m.v <- t(matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J))
y.m <- t(matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J))
l <- length(y.m)


nsims <- 100
timer <- matrix(0, nrow=nsims, ncol=4)
for(i in 1:nsims){
  ts <- proc.time()
  testv <- f.dea(f.inputs = x.m.f, outputs=y.m,
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  temp <- proc.time() - ts
  timer[i, 1:2] <- temp[2:3]
  ts <- proc.time()
  testf <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m,
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  temp <- proc.time() - ts
  timer[i, 3:4] <- temp[2:3]
}
summary(timer)


nsims <- 100
timer <- matrix(0, nrow=nsims, ncol=4)
for(i in 1:nsims){
  ts <- proc.time()
  testv <- f.dea.c(f.inputs = x.m.f, outputs=y.m,
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  temp <- proc.time() - ts
  timer[i, 1:2] <- temp[2:3]
  ts <- proc.time()
  testf <- f.dea.c(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m,
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  temp <- proc.time() - ts
  timer[i, 3:4] <- temp[2:3]
}
summary(timer)


## Testing if statements versus matrix creation
nsims <- 10000
ts <- proc.time()
for(i in 1:nsims){
  test <- matrix(0, ncol=nv, nrow=J)
}
proc.time() - ts

ts <- proc.time()
for(i in 1:nsims){
  if(nv == 11){
    test <- matrix(0, ncol=nv, nrow=J)
  }
}
proc.time() - ts


## Testing compiling
rm(list = ls())
library(FEAR)
source("dea_package.R")
library(compiler)
library(rbenchmark)

f.dea.c <- cmpfun(f.dea)

J <- 100
nv <- 2
nf <- 2
nu <- 2
x.m.f <- t(matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J))
x.m.v <- t(matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J))
y.m <- t(matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J))
l <- length(y.m)

nsims <- 100

benchmark(f.dea(f.inputs = x.m.f, outputs=y.m,
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO"),
          f.dea.c(f.inputs = x.m.f, outputs=y.m,
                  tech="V", orientation="OUT",
                  report.z="YES", slack="NO"),
          columns=c('test','replications','elapsed','relative'),
          order='relative',
          replications=nsims)


## Testing differences in non-convex JIM
rm(list = ls())
library(FEAR)
source("dea_package.R")
nsims <- 4

set.seed(4920)

J <- 100
nv <- 1
nf <- 2
nu <- 2
x.m.f <- t(matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J))
x.m.v <- t(matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J))
y.m <- t(matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J))
l <- length(y.m)

timer <- matrix(0, nrow=nsims, ncol=33)
for(i in 1:nsims){
  timer[i, 1:3] <- proc.time()[1:3]
  timer[i, 4:6] <- proc.time()[1:3]
  timer[i, 7:9] <- proc.time()[1:3]
  timer[i, 10:12] <- proc.time()[1:3]
  testf <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                 report.z="YES", slack="NO")
  timer[i, 13:15] <- proc.time()[1:3]
  testnc <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                  report.z="YES", slack="NO", convex="NO")

  timer[i, 16:18] <- proc.time()[1:3]
  act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)
  outf <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testf$theta, z=testf$z)
  outnc <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testnc$theta, z=testnc$z)

  timer[i, 19:21] <- proc.time()[1:3]
  jim.c <- f.jim(f.in.total=t(act$x.f.a), v.in.total=t(act$x.v.a), out.total=t(act$u.a),
                 f.in.opt=t(outf$x.f.s), v.in.opt=t(outf$x.v.s), out.opt=t(outf$u.s),
                 convex='YES')
  timer[i, 22:24] <- proc.time()[1:3]
  jim.nc1 <- f.jim(f.in.total=t(act$x.f.a), v.in.total=t(act$x.v.a), out.total=t(act$u.a),
                   f.in.opt=t(outnc$x.f.s), v.in.opt=t(outnc$x.v.s), out.opt=t(outnc$u.s),
                   convex='NO')
  timer[i, 25:27] <- proc.time()[1:3]
  jim.nc2 <- f.jim(f.in.total=t(act$x.f.a), v.in.total=t(act$x.v.a), out.total=t(0.25*act$u.a),
                   f.in.opt=t(outnc$x.f.s), v.in.opt=t(outnc$x.v.s), out.opt=t(0.25*outnc$u.s),
                   convex='NO')
  timer[i, 28:30] <- proc.time()[1:3]
  jim.nc3 <- f.jim(f.in.total=t(1.5*act$x.f.a), v.in.total=t(1.5*act$x.v.a), out.total=t(0.05*act$u.a),
                   f.in.opt=t(outnc$x.f.s), v.in.opt=t(outnc$x.v.s), out.opt=t(0.15*outnc$u.s),
                   convex='NO')
  timer[i, 31:33] <- proc.time()[1:3]  
}

timer <- data.frame(timer)
names(timer) <- c('1','2','3','4','5','6',
                  'DEA_slack','x1','x2','DEA_no_slack','x3','x4','DEA_convex','x5','x6',
                  'DEA_non_convex','x7','x8','Aggregation','x9','x10','JIM_convex','x11','x12',
                  'JIM_non_convex_1','x13','x14','JIM_non_convex_2','x15','x16','JIM_3','x17','x18')
timer.out <- timer
for(j in seq(7,33,3)){
  timer.out[,j:(j+2)] <- timer[,j:(j+2)] - timer[,(j-3):(j-1)]
}
colMeans(timer.out[,7:33])



## Russell
rm(list = ls())
library(FEAR)
source("dea_package.R")
nsims <- 4

set.seed(4620)

J <- 10
nv <- 1
nf <- 2
nu <- 2
x.m.f <- t(matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J))
x.m.v <- t(matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J))
y.m <- t(matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J))
l <- length(y.m)

testf <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
               report.z="YES", slack="NO")
testr <- f.russell(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", convex='YES')
print(testf$theta)
print(testr$rmo)


## Transposition
rm(list = ls())
library(FEAR)
source("dea_package.R")
library(compiler)
library(rbenchmark)

J <- 100
nv <- 10
nf <- 10
nu <- 10
x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
l <- length(y.m)

nsims <- 10

benchmark(f.dea(f.inputs = t(x.m.f), outputs=t(y.m),
                 tech="V", orientation="OUT",
                 report.z="YES", slack="NO"),
          f.dea.t(f.inputs = x.m.f, outputs=y.m,
                  tech="V", orientation="OUT",
                  report.z="YES", slack="NO"),
          columns=c('test','replications','elapsed','relative'),
          order='relative',
          replications=nsims)


## Multi-processor
rm(list = ls())
require(foreach)
require(doParallel)
library('linprog')
load("cu_interim_full.RData")
source("dea_package.R")

trials <- 1

cl <- makeCluster(6)
registerDoParallel(cl)


system.time( 
  test<-foreach(i=1:trials, .combine=c, .packages='linprog') %dopar% 
  {  
    theta.hat.nc.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                             orientation="OUT",
                             report.z="YES", convex='NO')
    opt <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                          theta=theta.hat$theta,
                          z=theta.hat$z)
    opt.ex <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                             theta=theta.hat.ex,
                             z=theta.hat.te$z)
    jim2 <- f.jim2(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
                 out.total=act$u.a,
                 f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
                 out.opt=opt$u.s, convex='YES')
    jim <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
                 out.total=act$u.a,
                 f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
                 out.opt=opt$u.s, convex='YES')
      
    theta.hat.c.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                             orientation="OUT",
                             report.z="YES", convex='YES')
    c(theta.hat.nc.te,theta.hat.c.te)
  })

stopCluster(cl)

system.time( 
  test3<-foreach(i=1:trials, .combine=c, .packages='linprog') %do% 
  {  
    theta.hat.nc.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                             orientation="OUT",
                             report.z="YES", convex='NO')
    opt <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                          theta=theta.hat$theta,
                          z=theta.hat$z)
    opt.ex <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                             theta=theta.hat.ex,
                             z=theta.hat.te$z)
    jim <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
                 out.total=act$u.a,
                 f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
                 out.opt=opt$u.s, convex='NO')
    
    theta.hat.c.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                             orientation="OUT",
                             report.z="YES", convex='YES')
    c(theta.hat.nc.te,theta.hat.c.te)
  })

test2 <- list()
system.time( 
  for(i in 1:trials){
      theta.hat.nc.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                               orientation="OUT",
                               report.z="YES", convex='NO')
    opt <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                          theta=theta.hat$theta,
                          z=theta.hat$z)
    opt.ex <- optimal.in.out(f.inputs=x.f, v.inputs=x.v, outputs=u,
                             theta=theta.hat.ex,
                             z=theta.hat.te$z)
    jim <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
                 out.total=act$u.a,
                 f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
                 out.opt=opt$u.s, convex='NO')
      
      theta.hat.c.te <- f.dea(f.inputs=cbind(x.f,x.v), outputs=u, tech="V",
                              orientation="OUT",
                              report.z="YES", convex='YES')
      test2[[i]] <- c(theta.hat.nc.te,theta.hat.c.te)
    })



ns <- 1
system.time(for(i in 1:ns){
  jim2 <- f.jim2(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
                 out.total=act$u.a,
                 f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
                 out.opt=opt$u.s, convex='NO')
})
system.time(for(i in 1:ns){
jim <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a,
             out.total=act$u.a,
             f.in.opt=opt$x.f.s, v.in.opt=opt$x.v.s,
             out.opt=opt$u.s, convex='NO')
})
