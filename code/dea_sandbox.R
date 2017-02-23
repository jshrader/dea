#!/usr/bin/R
#
# A test area for DEA analysis, using various mock and prepared datasets
#
# Jeff Shrader
# First version: 2011-12-29
# Time-stamp: "2014-11-14 15:32:04 jeff.shrader"
#

# Preliminaries:
library(FEAR)
source("dea_package.R")

## Walking through some examples, comparing my code to FEAR
# Make a simple dataset with fixed inputs
x.m.f <- matrix(c(10, 5, 3, 4, 5), ncol=1, nrow=5)
# Variable inputs
x.m.v.1 <- matrix(c(2, 50, 3, 2, 2), ncol=1, nrow=5)
# Outputs
y.m <- matrix(c(5, 0, 6, 1, 1,
                1, 4, 2, 0, 0), ncol=2, nrow=5)
# My code
cu <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v.1,
            outputs=y.m,
            tech="V", orientation="OUT", report.z="YES", convex="YES")
cu
# dea.graph(x.m.f, y.m, cu$theta)
dhat.m <- dea(XOBS=t(x.m.f), YOBS=t(y.m), ORIENTATION=2, RTS=1)
dhat.m
# dea.graph(x.m.f, y.m, dhat.m)

# You can see that the thetas are the same! Note that variable inputs are
# treated the same in FEAR as if you left them out (indeed, you can do the
# same in my code unless you want to run f.jim afterwards).


## Running more complicated models
data(ccr)
# Outputs
x <- t(matrix(ccr$x1, ncol=1, nrow=70))
# Inputs
y <- t(matrix(ccr$y1, ncol=1, nrow=70))

# Estimating DEA with output orientation, one input, one output
dhat.vrs <- dea(XOBS=x, YOBS=y, ORIENTATION=2)

# plotting results a la Fare et al
plot(x, y)
eff.res <- matrix(c(x, y*(1/dhat.vrs)), ncol=2, nrow=70)[order(x),]
points(eff.res[,1], eff.res[,2], col='blue', type='l')

# Mock data from Table 2.1, pg 26 in Cooper, Seiford, and Tone
x.m <- t(matrix(c(2, 3, 3, 4, 5, 5, 6, 8), ncol=1, nrow=8))
y.m <- t(matrix(c(1, 3, 2, 3, 4, 2, 3, 5), ncol=1, nrow=8))

# Calculating efficiency of mock data with FEAR
dhat.m <- dea(XOBS=x.m, YOBS=y.m, ORIENTATION=2, RTS=1)
plot(x.m, y.m)
eff.res <- matrix(c(x.m, y.m*(1/dhat.m)), ncol=2, nrow=length(x.m))[order(x.m),]
points(eff.res[,1], eff.res[,2], col='blue', type='l')

## Mock data for fixed and variable inputs
x.m.f <- t(matrix(c(2, 3, 3, 4, 5, 5, 6, 8), ncol=1, nrow=8))
x.m.v <- t(matrix(c(1, 3, 5, 2, 3, 5, 1, 4), ncol=1, nrow=8))
y.m <- t(matrix(c(1, 3, 2, 3, 4, 2, 3, 5), ncol=1, nrow=8))

f.dea(inputs = rbind(x.m.f, x.m.v), outputs=y.m, tech="V")
f.dea(inputs=x.m.f, outputs=y.m, tech="V")
f.dea(inputs=x.m.v, outputs=y.m, tech="V")
viur.1vin(var.inputs=x.m.v, fix.inputs=x.m.f, outputs=y.m, tech="V")

x.m.f2 <- t(matrix(c(2, 0, 3, 4, 9, 5, 10, 8,
                     0, 3, 15, 2, 6, 7, 22, 1), ncol=2, nrow=8))
x.m.v2 <- t(matrix(c(1, 3, 5, 2, 3, 5, 1, 4,
                     10, 100, 13, 16, 17, 19, 16, 200), ncol=2, nrow=8))
y.m2 <- t(matrix(c(1, 3, 2, 3, 4, 50, 3, 59,
                   10, 14, 8, 4, 7, 34, 1, 6), ncol=2, nrow=8))
f.dea(inputs = rbind(x.m.f2, x.m.v2), outputs=y.m2, tech="V")
f.dea(inputs=x.m.f2, outputs=y.m2, tech="V")
f.dea(inputs=x.m.v2, outputs=y.m2, tech="V")
theta.viur <- dea.viur(var.inputs=x.m.v2, inputs=x.m.f2, outputs=y.m2, tech="V")
theta.dea <- dea.viur(inputs=rbind(x.m.f2,x.m.v2), outputs=y.m2, tech="V")
pc.utilization(theta.all=theta.dea, theta.var=theta.viur[,1])




# My own DEA solver
library(linprog)

dea.f <- function(inputs, outputs, tech="V"){
  # Function to calculate efficiency of firms based on Fare 1994 ch.4
  #
  # input/outputs in the form of rows = in/out, colums = firms
  #
  # To Do:
  # . Allow for constant returns to scale
  
  J <- length(inputs[1,])
  N <- length(inputs[,1])
  M <- length(outputs[,1])
  f.obj <- c(1, rep(0, J))
  f.in  <- cbind(0, inputs)
  f.neg <- cbind(0, diag(J))
  
  # Technology type determines whether additional constraints are put on
  # z's, the intensity variables
  if(tech == "N"){
    tec.dir <- "<="
  }else if(tech == "V"){
    tec.dir <- "="
  }else{
    tec.dir <- ""
  }
  f.int <- t(matrix(c(0, rep(1, J))))
                 
  f.dir <- c(rep("<=",M), rep("<=",N), rep(">=",J), tec.dir)
  f.rhs <- rep(0, M+N+J)

  theta <- rep(0, J)
  for(j in 1:J){
    f.out <- cbind(outputs[,j], -outputs)
    f.con <- rbind(f.out, f.in, f.neg, f.int)
    f.rhs <- c(rep(0,M) , inputs[,j], rep(0,J), 1)
    theta[j] <- 1/lp("max", f.obj, f.con, f.dir, f.rhs)$solution[1]
  }
  return(theta)
}

dhat.t <- dea.f(inputs = x.m, outputs = y.m, tech="V")

# Further tests using mock data from Cooper
x.3.2 <- t(matrix(c(4, 7, 8, 4, 2, 10, 3,
                    3, 3, 1, 2, 4, 1, 7
                    ), ncol=2, nrow=7))
y.3.2 <- t(matrix(c(1, 3, 4, 8, 4, 5, 2,
                    0, 1, 3, 0, 2, 1, 3), ncol=2, nrow=7)) # slightly changed
dhat.t <- dea.f(inputs=x.3.2, outputs=y.3.2, tech="V")
dhat.t2 <- dea(XOBS=x.3.2, YOBS=y.3.2, ORIENTATION=2)

# Test against larger dataset
dhat.t <- dea.f(inputs=x, outputs=y, tech="V")
dhat.t2 <- dea(XOBS=x, YOBS=y, ORIENTATION=2)








## Manual tests of alternative methods

# DMU A Cooper method
f.obj <- c(1, 0)
f.con <- matrix(c(0, 2,
                  1, -2,
                  3, -3,
                  2, -3,
                  3, -4,
                  4, -5,
                  2, -5,
                  3, -6,
                  5, -8
                  ),
                nrow=9, byrow=TRUE)
f.dir <- c("=","<=","<=","<=","<=","<=","<=","<=","<=")
f.rhs <- c(1,0,0,0,0,0,0,0,0)
lp("max", f.obj, f.con, f.dir, f.rhs)$solution*1


# DMU A Fare method
f.obj <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
f.con <- matrix(c(5, -1, -3, -2, -3, -4, -2, -3, -5,
                  0,  2,  3,  3,  4,  5,  5,  6,  8,
                  0,  1,  0,  0,  0,  0,  0,  0,  0,
                  0,  0,  1,  0,  0,  0,  0,  0,  0,
                  0,  0,  0,  1,  0,  0,  0,  0,  0,
                  0,  0,  0,  0,  1,  0,  0,  0,  0,
                  0,  0,  0,  0,  0,  1,  0,  0,  0,
                  0,  0,  0,  0,  0,  0,  1,  0,  0,
                  0,  0,  0,  0,  0,  0,  0,  1,  0,
                  0,  0,  0,  0,  0,  0,  0,  0,  1
                  ),
                nrow=10, byrow=TRUE)
f.dir <- c("<=","<=",">=",">=",">=",">=",">=",">=",">=",">=")
f.rhs <- c(0,8,0,0,0,0,0,0,0,0)
1/lp("max", f.obj, f.con, f.dir, f.rhs)$solution[1]

             
