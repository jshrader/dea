#!/usr/bin/R
## A collection of DEA and efficiency analysis tools. Created originally for
## analyzing purse seine tuna fishers in the Pacific.
##
## Jeff Shrader
## Created On: 2011-01-04
## Current Version:
##   Time-stamp: "2012-09-21 15:27:48 jgs"
##
## To Do:
## . allow for more than 1 input and output in graphs -- either through
##   heat map/multidimensional graphs or by just graphing inputs or outputs
##    -- see, for example, Fare 1994 pg 29.
## . Speed is a big issue as is size of the vectors that I need to assign
## .. Currently, I can only handle about 1000 observations at a time with my
##    own functions. FEAR does better, being written in FORTRAN, but it does
##    not calculate VIUR, Russel non-radial measures, or Johansen industry models
##    directly. Some of these can be kludged out by repeatedly running dea()
##    with the appropriate inputs and outputs, but some things like non-
##    convexity cannot be done. 
## . Constant returns to scale doesn't work because the variables are still
##   added to the constraint but the inequality is removed.
##

#### Preliminaries ####
library(linprog)
library(foreach)
library(doParallel)

#### DEA calculations ####
f.dea <- function(f.inputs, v.inputs=matrix(NA, ncol=1, nrow=1),
                  outputs, tech='V', orientation='OUT',
                  report.z='NO', slack='NO',
                  convex='YES'){
  ## Function to calculate efficiency of firms based on Fare, Grosskopf,
  ## and Lovell, 1994. Input orientation is described in ch.3 eq
  ## 3.1.12. Output orientation is described in ch.4 eq 4.1.11.
  ## Variable input utilization is described in ch. 10, eq 10.3.3
  ## and 10.3.4/5.
  ##
  ## For slack variables, I consulted the margin notes of Dale's Fare 1994
  ## book. I am sure there is another source one could find.
  ##
  ## input/outputs in the form of columns = in/out, rows = firms

  ## Testing
  if(is.na(v.inputs[1,1])==FALSE){ # Testing for presence of variable inputs
    has.v <- 1
    if(orientation == "IN"){
      stop(cat("Variable input model can only be estimated with output orientation.","\n",
               "Please choose 'OUT' for your orientation","\n"))
    }
  }else{
    has.v <- 0
  }

  ## Index initialization 
  J <- length(f.inputs[,1]) # Number of firms
  M <- length(outputs[1,]) # Number of outputs  
  N.f <- length(f.inputs[1,]) # Number of fixed inputs
  if(has.v == 1){
    N.v <- length(v.inputs[1,]) # Number of variable inputs
    f.in.v <- rbind(0, v.inputs, matrix(0, nrow=N.v, ncol=N.v)) # Variable inputs,
                          # needs replacement of jth term of variable input matrix
    if(slack=='YES'){
      f.in.v <- rbind(f.in.v, matrix(0, nrow=(M+N.f), ncol=N.v))
    }
  }else{
    N.v <- 0
  }

  ## Objective function initialization
  if(slack=='NO'){
    f.obj <- c(1, rep(0, J + N.v)) # Objective function (1*theta + 0*z_j + 0*lambda_jn)
  }else{
    ## Objective function (1*theta + 0*z_j + 0*lambda_jn + 0*s + 0*e)    
    f.obj <- c(1, rep(0, J + N.v + M + N.f))
  }

  ## Constraints
  if(orientation == "IN"){
    f.out <- rbind(0, -outputs)
  }else if(orientation == "OUT"){
    if(has.v == 0){
      f.in.f  <- rbind(0, f.inputs) # Just Fixed inputs
    }else{
      f.in.f  <- rbind(0, f.inputs, matrix(0, nrow=N.v, ncol=N.f)) # Fixed inputs with variable
    }
    if(slack=='YES'){
      f.in.f <- rbind(f.in.f, matrix(0, nrow=M, ncol=N.f), diag(N.f))
    }
  }
  if(convex=='YES'){
    f.neg <- rbind(0, diag(J + N.v))
  }else{
    if(has.v==1){
      f.neg <- rbind(0, matrix(0, nrow=J, ncol=N.v), diag(N.v))      
    }
  }
  f.int <- matrix(c(0, rep(1, J), rep(0, N.v)))
  if(slack=='YES'){
    f.neg <- rbind(f.neg, matrix(0, nrow=(M+N.f), ncol=(J+N.v)))
    f.int <- matrix(c(f.int, rep(0, (M+N.f))))
  }

  ## Technology type determines whether additional constraints are put on
  ## z's, the intensity variables
  if(slack=='NO'){
    f.dir <- c(rep("<=",M), rep("<=",N.f), rep("=",N.v))#, rep(">=",J+N.v))
  }else{
    f.dir <- c(rep("=",M), rep("=",N.f), rep("=",N.v))#, rep(">=",J+N.v))    
  }
  if(convex=='YES'){
    f.dir <- c(f.dir, rep(">=", J+N.v))
  }else{
    f.dir <- c(f.dir, rep(">=", N.v))
  }
  if(tech == "N"){         # Non-increasing
    f.dir <- c(f.dir, "<=")
  }else if(tech == "V"){   # Variable
    f.dir <- c(f.dir, "=")
  }else if(tech == "C"){   # Constant
    f.dir <- f.dir
  }else{
    stop(cat("Tech type is not recognized. Use one of the following:","\n",
             " N (non-increasing), V (variable), or C (constant)","\n"))
  }

  theta <- rep(0, J)
  lambda <- matrix(0, ncol=N.v, nrow=J)
  z <- matrix(0, ncol=J, nrow=J)
  s <- matrix(0, ncol=M, nrow=J)
  e <- matrix(0, ncol=N.f, nrow=J)
  for(j in 1:J){
    if(orientation == "IN"){
      f.in.f <- rbind(-f.inputs[j,], f.inputs)
      if(convex == 'YES'){
        f.rhs <- c(-outputs[j,] , rep(0,N.f), rep(0,J), 1)
        f.con <- cbind(f.out, f.in.f, f.neg, f.int)
        optimum <- lp("min", f.obj, f.con, f.dir, f.rhs,
                      transpose.constraints=FALSE)
      }else{
        f.rhs <- c(-outputs[j,] , rep(0,N.f), 1)
        f.con <- cbind(f.out, f.in.f, f.int)
        optimum <- lp("min", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)),
                      transpose.constraints=FALSE)
      }
      ## print(optimum)
      theta[j] <- 1/optimum$solution[1]
    }else if(orientation == "OUT"){
      if(has.v == 0){
        if(slack=='NO'){
          f.out <- rbind(outputs[j,], -outputs)
        }else{
          f.out <- rbind(outputs[j,], -outputs, diag(M),
                         matrix(0, nrow=N.f, ncol=M))
        }
        if(convex=='YES'){
          f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0,J), 1)
          f.con <- cbind(f.out, f.in.f, f.neg, f.int)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs,
                        transpose.constraints=FALSE)
        }else{
          f.rhs <- c(rep(0,M) , f.inputs[j,], 1)
          f.con <- cbind(f.out, f.in.f, f.int)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)),
                        transpose.constraints=FALSE)
        }
        ## print(optimum)
        theta[j] <- 1/optimum$solution[1]
        z[j,] <- optimum$solution[2:(1+J)]
        if(slack=="YES"){
          s[j,] <- optimum$solution[(2+J):(1+J+M)]
          e[j,] <- optimum$solution[(2+J+M):(1+J+M+N.f)]          
        }
      }else{
        if(slack=='NO'){
          f.out <- rbind(outputs[j,], -outputs, matrix(0, nrow=N.v, ncol=M))
        }else{
          f.out <- rbind(outputs[j,], -outputs,
                         matrix(0, nrow=N.v, ncol=M),
                         diag(M), matrix(0, nrow=N.f, ncol=M))
        }
        f.in.v[(1+J+1):(1+J+N.v),] <- -diag(N.v)*v.inputs[j,]
        if(convex=='YES'){
          f.con <- cbind(f.out, f.in.f, f.in.v, f.neg, f.int)
          f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0, N.v+N.v+J), 1)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs,
                        transpose.constraints=FALSE)
        }else{
          f.con <- cbind(f.out, f.in.f, f.in.v, f.neg, f.int)
          f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0, N.v+N.v), 1)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)),
                        transpose.constraints=FALSE)
        }
        ## print(optimum)        
        theta[j] <- 1/optimum$solution[1]
        lambda[j,] <- optimum$solution[(1+J+1):(1+J+N.v)]
        z[j,] <- optimum$solution[2:(1+J)]
        if(slack=="YES"){
          s[j,] <- optimum$solution[(2+J+N.v):(1+J+N.v+M)]
          e[j,] <- optimum$solution[(2+J+N.v+M):(1+J+N.v+M+N.f)]          
        }
      }
    }
  }
  ##return(optimum)  
  rm(optimum)
  ## print(f.obj)
  ## print(f.con)
  ## print(f.dir)
  ## print(f.rhs)
  if(orientation == 'OUT'){
    if(report.z == 'NO' & has.v == 0 & slack=='NO'){
      list(theta=theta)
    }else if(report.z == 'YES' & has.v == 0 & slack=='NO'){
      list(theta=theta, z=z)
    }else if(report.z == 'NO' & has.v == 1 & slack=='NO'){
      list(theta=theta, lambda=lambda)
    }else if(report.z == 'YES' & has.v == 1 & slack=='NO'){
      list(theta=theta, z=z, lambda=lambda)
    }else if(report.z == 'YES' & has.v == 0 & slack=='YES'){
      list(theta=theta, z=z, s=s, e=e)
    }else if(report.z == 'YES' & has.v == 1 & slack=='YES'){
      list(theta=theta, z=z, lambda=lambda, s=s, e=e)
    }else if(report.z == 'NO' & has.v == 1 & slack=='YES'){
      list(theta=theta, lambda=lambda, s=s, e=e)
    }else if(report.z == 'NO' & has.v == 0 & slack=='YES'){
      list(theta=theta, s=s, e=s)
    }
  }else if(orientation == 'IN'){
    list(theta=theta)
  }
}


f.russell <- function(f.inputs, v.inputs=matrix(NA, ncol=1, nrow=1),
                     outputs, tech='V', convex='YES'){
  ## Function to calculate efficiency using the non-radial Russel measure
  ## as presented in Fare et al 1994 pg 116.
  ##
  ## input/outputs in the form of columns = in/out, rows = firms

  ## Testing
  if(is.na(v.inputs[1,1])==FALSE){ # Testing for presence of variable inputs
    has.v <- 1
  }else{
    has.v <- 0
  }

  ## Index initialization 
  J <- length(outputs[,1]) # Number of firms
  M <- length(outputs[1,]) # Number of outputs  
  N.f <- length(f.inputs[1,]) # Number of fixed inputs
  if(has.v == 1){
    N.v <- length(v.inputs[1,]) # Number of variable inputs
    f.in.v <- rbind(matrix(0, nrow=M, ncol=N.v), v.inputs,
                    matrix(0, nrow=N.v, ncol=N.v)) # Variable inputs,
                          # needs replacement of jth term of variable input matrix
  }else{
    N.v <- 0
  }

  ## Objective function initialization
  f.obj <- c(rep(1, M), rep(0, J + N.v)) # Objective function ((1/M)*theta_m + 0*z_j + 0*lambda_jn)

  ## Constraints
  if(has.v == 0){
    f.in.f  <- rbind(matrix(0, nrow=M, ncol=N.f), f.inputs) # Just Fixed inputs
  }else{
    f.in.f  <- rbind(matrix(0, nrow=M, ncol=N.f), f.inputs,
                     matrix(0, nrow=N.v, ncol=N.f)) # Fixed inputs with variable
  }
  if(convex=='YES'){
    f.neg <- rbind(matrix(0, nrow=M, ncol=(J+N.v)), diag(J + N.v))
  }else{
    if(has.v==1){
      f.neg <- rbind(matrix(0, nrow=M, ncol=N.v)
                     , matrix(0, nrow=J, ncol=N.v), diag(N.v))      
    }
  }
  f.int <- matrix(c(rep(0, M), rep(1, J), rep(0, N.v)))

  ## Technology type determines whether additional constraints are put on
  ## z's, the intensity variables
  f.dir <- c(rep("<=",M), rep("<=",N.f), rep("=",N.v))
  if(convex=='YES'){
    f.dir <- c(f.dir, rep(">=", J+N.v))
  }else{
    f.dir <- c(f.dir, rep(">=", N.v))
  }
  if(tech == "N"){         # Non-increasing
    f.dir <- c(f.dir, "<=")
  }else if(tech == "V"){   # Variable
    f.dir <- c(f.dir, "=")
  }else if(tech == "C"){   # Constant
    f.dir <- f.dir
  }else{
    stop(cat("Tech type is not recognized. Use one of the following:","\n",
             " N (non-increasing), V (variable), or C (constant)","\n"))
  }


  rmo <- matrix(0, ncol=M, nrow=J)
  lambda <- matrix(0, ncol=N.v, nrow=J)
  z <- matrix(0, ncol=J, nrow=J)
  for(j in 1:J){
    if(has.v == 0){
      f.out <- rbind(diag(M)*outputs[j,], -outputs)
      if(convex=='YES'){
        f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0,J), 1)
        f.con <- cbind(f.out, f.in.f, f.neg, f.int)
        optimum <- lp("max", f.obj, f.con, f.dir, f.rhs,
                      transpose.constraints=FALSE)
      }else{
        f.rhs <- c(rep(0,M) , f.inputs[j,], 1)
        f.con <- cbind(f.out, f.in.f, f.int)
        optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=((M+1):(1+J)),
                      transpose.constraints=FALSE)
      }
    }else{
      f.out <- rbind(diag(M)*outputs[j,], -outputs, matrix(0, nrow=N.v, ncol=M))
      f.in.v[(M+J+1):(M+J+N.v),] <- -diag(N.v)*v.inputs[j,]
      if(convex=='YES'){
        f.con <- cbind(f.out, f.in.f, f.in.v, f.neg, f.int)
        f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0, N.v+N.v+J), 1)
        optimum <- lp("max", f.obj, f.con, f.dir, f.rhs,
                      transpose.constraints=FALSE)
      }else{
        f.con <- cbind(f.out, f.in.f, f.in.v, f.neg, f.int)
        f.rhs <- c(rep(0,M) , f.inputs[j,], rep(0, N.v+N.v), 1)
        optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=((M+1):(1+J)),
                      transpose.constraints=FALSE)
      }
      ## print(optimum$objval)
      lambda[j,] <- optimum$solution[(1+J+M):(M+J+N.v)]
    }
    rmo[j,] <- optimum$solution[1:M]
    z[j,] <- optimum$solution[(1+M):(J+M)]    
  }
  theta <- 1/rowMeans(rmo, na.rm=TRUE)
  ## print(f.obj)
  ## print(f.con)
  ## print(f.dir)
  ## print(f.rhs)
  if(has.v == 0){
    list(theta=theta, rmo=rmo, z=z)
  }else if(has.v == 1){
    list(theta=theta, rmo=rmo, z=z, lambda=lambda)
  }
  ##return(optimum)
}


pc.utilization <- function(theta.all, theta.var){
  pcum <- theta.var/theta.all # We define it as the recipricol of Fare's measure because I think
                                        # it makes more sense to have output-oriented measures be less
                                        # than 1 since it is the current output relative to frontier
                                        # output, which is always at least as large as current
  return(data.frame(pcum))
}


#### Johansen Industry Model Estimation ####
optimal.in.out <- function(f.inputs, v.inputs=matrix(NA, ncol=1, nrow=1),
                           outputs, theta, z, flow=c(0), flow.var=1){
  ## You must supply output oriented intensity variables and slack
  ## variables to this function by running
  ## f.dea(., ., ., tech="V", orientation="OUT", report.z="YES", slack="YES")
  ##
  ## f.inputs: fixed inputs, same as provided to f.dea
  ## v.inputs: variable inputs, same as provided to f.dea
  ## outputs:
  ## theta: firm efficiency measure, output from f.dea
  ## z: firm intensity variables, output from f.dea
  ## flow: an optional list of index values for flow variables in f.inputs
  ## flow.var: optional index for which variable is the flow variable. If
  ##           unspecified but flow != 0, then the first variable is
  ##           assumed to be the flow variable.
  ##
  theta <- as.matrix(theta)
  z <- as.matrix(z)
  v.inputs <- as.matrix(v.inputs)
  f.inputs <- as.matrix(f.inputs)
  outputs <- as.matrix(outputs)

  ## If column dim == 1, then it is radial DEA output
  ## If dim > 1, then it is Russell non-radial measure  
  if(dim(theta)[2] > 1){
    theta <- dim(theta)[1]/theta
  }else{
    theta <- c(theta)
  }
  u.s <- outputs/theta
  ##u.s <- ifelse(is.na(u.s), 0, u.s)
  x.f.s <- f.inputs

  if(is.na(v.inputs[1,1])==FALSE){
    ## Has variable inputs
    x.v.s <- z%*%v.inputs
    if(flow[1]!=0){
      ## TO DO!
      ## x.f.s[,flow] <- f.inputs[,flow]/(theta*x.v.s[,flow.var])
      ## x.f.s[,flow] <- f.inputs[is.na(theta) == FALSE,flow]/(theta[is.na(theta) == FALSE]*v.inputs[is.na(theta) == FALSE,flow.var])      
    }
    list(u.s=u.s, x.f.s=x.f.s, x.v.s=x.v.s)
  }else{
    list(u.s=u.s, x.f.s=x.f.s)
  }
}

industry.observed <- function(f.inputs, v.inputs=matrix(NA, ncol=1, nrow=1),
                              outputs, flow=c(0), flow.var=1){
  u.a <- colSums(outputs)
  x.f.a <- colSums(f.inputs)
  if(is.na(v.inputs[1,1])==FALSE){
    x.v.a <- colSums(v.inputs)
    if(flow[1]!=0){
      x.f.a[flow] <- colSums(f.inputs[,flow]/v.inputs[,flow.var])
    }
    list(u.a=u.a, x.f.a=x.f.a, x.v.a=x.v.a)
  }else{
    list(u.a=u.a, x.f.a=x.f.a)
  }
}

f.jim <- function(f.in.total, v.in.total=matrix(NA, ncol=1, nrow=1), out.total,
                  f.in.opt, v.in.opt=matrix(NA, ncol=1, nrow=1), out.opt,
                  tac=0, fd.max=0, ivq=NULL, ivq.var=NULL, convex='YES'){
  ## Function to calculate Johansen industry model and minimum fixed
  ## inputs following Kerstens et al 2006.
  ##
  ## input/outputs come from the industry.observed and optimal.in.out
  ## functions. 
  ##
  ## tac should be specified for each output
  ## fd.max fixes a maximum allotment of the first variable input for all
  ##  vessels
  ## ivq specifies a vessel specific restriction on outputs
  ## ivq.var gives the variable index that the ivq should be applied to
  ## convex determines whether the frontier is convex (YES) or no (NO)

  ## tin <- proc.time()
  
  ## Testing
  if(!is.na(v.in.opt[1,1])){ # Testing for presence of variable inputs
    has.v <- 1
  }else{
    has.v <- 0
  }

  ## Index specification
  J <- length(out.opt[,1]) # Number of firms in industry
  N.f <- length(f.in.total) # Number of fixed inputs
  N.v <- length(v.in.total) # Number of variable inputs
  M <- length(out.total) # Number of outputs

  ## Setting TAC if appropriate
  if(tac[1] > 0){
    for(i in 1:M){
      out.total[i] <- min(out.total[i], tac[i])
    }
  }

  ## Objective function
  f.obj <- c(1, rep(0, J + N.v)) # Objective function (1*theta + 0*w_j + 0*X_v)

  ## Constraints
  f.out <- rbind(0, out.opt, matrix(0, nrow=N.v, ncol=M)) # Outputs
  f.in.f  <- rbind(matrix(-f.in.total, nrow=1, ncol=N.f), f.in.opt,
                   matrix(0, nrow=N.v, ncol=N.f)) # Fixed inputs
  if(length(v.in.total)==1){
    f.in.v <- rbind(0, v.in.opt, -v.in.total) # Variable inputs
  }else if(length(v.in.total)>1){
    ## TO DO: Check that this has the correct dimension
    f.in.v <- rbind(0, v.in.opt, -diag(as.vector(t(v.in.total))))
                                        # Variable inputs
  }
  f.neg.w <- rbind(0, diag(J), matrix(0, nrow=N.v, ncol=J))
  f.one.w <- f.neg.w
  f.neg.t <- matrix(c(1, rep(0, J), rep(0, N.v)), nrow=(1+J+N.v), ncol=1)
  if(fd.max > 0){
    f.day <- rbind(0, v.in.opt[,1]*diag(J), matrix(0, nrow=N.v, ncol=J))
  }
  if(!is.null(ivq)){
    f.ivq <- rbind(0, out.opt[,ivq.var]*diag(J), matrix(0, nrow=N.v, ncol=J))
  }

  ## Direction of constraint (inequality/equality and in what direction)
  if(convex=='YES'){
    f.dir <- c(rep(">=",M), rep("<=",N.f), rep("<=",N.v), rep(">=",J),
               rep("<=",J), rep(">=",1))
  }else{
    f.dir <- c(rep(">=",M), rep("<=",N.f), rep("<=",N.v), rep(">=",1))
  }
  if(fd.max > 0){
    f.dir <- c(f.dir, rep("<=",J))
  }
  if(!is.null(ivq)){
    f.dir <- c(f.dir, rep("<=",J))
  }
  

  ## Right hand side values for each constraint
  if(convex=='YES'){
    f.rhs <- c(out.total, rep(0,N.f), rep(0,N.v), rep(0,J), rep(1,J), 0)
  }else{
    f.rhs <- c(out.total, rep(0,N.f), rep(0,N.v), 0)
  }
  if(fd.max > 0){
    f.rhs <- c(f.rhs, rep(fd.max, J))
  }
  if(!is.null(ivq)){
    f.rhs <- c(f.rhs, ivq)
  }


  theta <- 0
  w <- rep(0,J)
  Xv <- rep(0, N.v)

  ## Fitting
  if(convex=='YES'){
    f.con <- cbind(f.out, f.in.f, f.in.v, f.neg.w, f.one.w, f.neg.t)
  }else{
    f.con <- cbind(f.out, f.in.f, f.in.v, f.neg.t)
  }
  if(fd.max > 0){
    f.con <- cbind(f.con, f.day)
  }
  if(!is.null(ivq)){
    f.con <- cbind(f.con, f.ivq)
  }
  ## print('obj')
  ## print(f.obj)
  ## print('constraint')
  ## print(f.con)
  ## print('f.dir')
  ## print(f.dir)
  ## print('rhs')
  ## print(f.rhs)
  ## print(proc.time() - tin)
  ## tin <- proc.time()
  if(convex=='YES'){
    optimum <- lp("min", f.obj, f.con, f.dir, f.rhs,
                  transpose.constraints=FALSE)
  }else{
    optimum <- lp("min", f.obj, f.con, f.dir, f.rhs, binary.vec=2:(J+1),
                  transpose.constraints=FALSE)
  }
  ## print(proc.time() - tin)
  ## print(optimum)
  theta <- optimum$solution[1]
  w <- optimum$solution[2:(1+J)]
  Xv <- optimum$solution[(2+J):(1+J+N.v)]
  rm(optimum)
  return(list(theta=theta, w=w, Xv=Xv))
}


#### Graphing ####
dea.graph <- function(inputs, outputs, theta){
  ## Simple graphing along the lines of Fare 1994
  ## To Do:
  ## . add support for more dimensions
  ## .. perhaps by doing marginal plots? Grid of plots?
  plot(inputs[,1], outputs[,1])
  print(length(inputs[,1]))
  print(length(outputs[,1]))
  J <- length(inputs[,1])
  eff.res <- matrix(c(inputs[,1], outputs[,1]*(1/theta)), ncol=2, nrow=J)[order(inputs[,1]),]
  points(eff.res[,1], eff.res[,2], col='blue', type='l')
}

#### Tables ####
dea.tex.table <- function(est, header=FALSE){
  ## Making TeX table out of a generic matrix
  ## To Do:
  ## . allow arbitrary input width
  ## . break long tables appropriately
  J <- length(matrix(est)[,1]) #number of DMUs
  n <- length(matrix(est)[1,]) #number of columns in table
  table.in=matrix(nrow=J,ncol=n)
  table.in[,1]=c(1:n)
  for(i in 2:n+1){
    table.in[,i]=est[,1]
  }

  table.in[1:9,1]=paste(" ",table.in[1:9,1],sep="")
  table.in[,2]=ifelse(nchar(table.in[,2])==1,
            paste(table.in[,2],".",sep=""), table.in[,2])

  table.in[,2:7]=paste(table.in[,2:7],"000000",sep="")
  table.in[,c(2:3,5:7)]=substr(table.in[,c(2:3,5:7)],1,6)
  table.in[,4]=substr(table.in[,4],1,7)
  table.in=paste(table.in[,1]," & ",
    table.in[,2]," & ",
    table.in[,3]," & ",
    table.in[,4]," & ",
    table.in[,5]," & ",
    table.in[,6]," & ",
    table.in[,7]," \\",sep="")
}





#### Deprecated Functions ####

f.dea.t <- function(f.inputs, v.inputs=matrix(NA, ncol=1, nrow=1),
                  outputs, tech='V', orientation='OUT',
                  report.z='NO', slack='NO',
                  convex='YES'){
  print("Deprecated: use f.dea instead and reverse the orientation of your input matrices")
  ## Function to calculate efficiency of firms based on Fare, Grosskopf,
  ## and Lovell, 1994. Input orientation is described in ch.3 eq
  ## 3.1.12. Output orientation is described in ch.4 eq 4.1.11.
  ## Variable input utilization is described in ch. 10, eq 10.3.3
  ## and 10.3.4/5.
  ##
  ## For slack variables, I consulted the margin notes of Dale's Fare 1994
  ## book. I am sure there is another source one could find.
  ##
  ## input/outputs in the form of rows = in/out, colums = firms

  ## Testing
  if(is.na(v.inputs[1,1])==FALSE){ # Testing for presence of variable inputs
    has.v <- 1
    if(orientation == "IN"){
      stop(cat("Variable input model can only be estimated with output orientation.","\n",
               "Please choose 'OUT' for your orientation","\n"))
    }
  }else{
    has.v <- 0
  }

  ## Index initialization 
  J <- length(f.inputs[1,]) # Number of firms
  M <- length(outputs[,1]) # Number of outputs  
  N.f <- length(f.inputs[,1]) # Number of fixed inputs
  if(has.v == 1){
    N.v <- length(v.inputs[,1]) # Number of variable inputs
    f.in.v <- cbind(0, v.inputs, matrix(0, nrow=N.v, ncol=N.v)) # Variable inputs,
                          # needs replacement of jth term of variable input matrix
    if(slack=='YES'){
      f.in.v <- cbind(f.in.v, matrix(0, nrow=N.v, ncol=(M+N.f)))
    }
  }else{
    N.v <- 0
  }

  ## Objective function initialization
  if(slack=='NO'){
    f.obj <- c(1, rep(0, J + N.v)) # Objective function (1*theta + 0*z_j + 0*lambda_jn)
  }else{
    ## Objective function (1*theta + 0*z_j + 0*lambda_jn + 0*s + 0*e)    
    f.obj <- c(1, rep(0, J + N.v + M + N.f))
  }

  ## Constraints
  if(orientation == "IN"){
    f.out <- cbind(0, -outputs)
  }else if(orientation == "OUT"){
    if(has.v == 0){
      f.in.f  <- cbind(0, f.inputs) # Just Fixed inputs
    }else{
      f.in.f  <- cbind(0, f.inputs, matrix(0, nrow=N.f, ncol=N.v)) # Fixed inputs with variable
    }
    if(slack=='YES'){
      f.in.f <- cbind(f.in.f, matrix(0, nrow=N.f, ncol=M), diag(N.f))
    }
  }
  if(convex=='YES'){
    f.neg <- cbind(0, diag(J + N.v))
  }else{
    if(has.v==1){
      f.neg <- cbind(0, matrix(0, nrow=N.v, ncol=J), diag(N.v))      
    }
  }
  f.int <- t(matrix(c(0, rep(1, J), rep(0, N.v))))
  if(slack=='YES'){
    f.neg <- cbind(f.neg, matrix(0, nrow=(J+N.v), ncol=(M+N.f)))
    f.int <- c(f.int, rep(0, (M+N.f)))    
  }

  ## Technology type determines whether additional constraints are put on
  ## z's, the intensity variables
  if(slack=='NO'){
    f.dir <- c(rep("<=",M), rep("<=",N.f), rep("=",N.v))#, rep(">=",J+N.v))
  }else{
    f.dir <- c(rep("=",M), rep("=",N.f), rep("=",N.v))#, rep(">=",J+N.v))    
  }
  if(convex=='YES'){
    f.dir <- c(f.dir, rep(">=", J+N.v))
  }else{
    f.dir <- c(f.dir, rep(">=", N.v))
  }
  if(tech == "N"){         # Non-increasing
    f.dir <- c(f.dir, "<=")
  }else if(tech == "V"){   # Variable
    f.dir <- c(f.dir, "=")
  }else if(tech == "C"){   # Constant
    f.dir <- f.dir
  }else{
    stop(cat("Tech type is not recognized. Use one of the following:","\n",
             " N (non-increasing), V (variable), or C (constant)","\n"))
  }


  theta <- rep(0, J)
  lambda <- matrix(0, ncol=N.v, nrow=J)
  z <- matrix(0, ncol=J, nrow=J)
  s <- matrix(0, ncol=M, nrow=J)
  e <- matrix(0, ncol=N.f, nrow=J)
  for(j in 1:J){
    if(orientation == "IN"){
      f.in.f <- cbind(-f.inputs[,j], f.inputs)
      if(convex == 'YES'){
        f.rhs <- c(-outputs[,j] , rep(0,N.f), rep(0,J), 1)
        f.con <- rbind(f.out, f.in.f, f.neg, f.int)
        optimum <- lp("min", f.obj, f.con, f.dir, f.rhs)
      }else{
        f.rhs <- c(-outputs[,j] , rep(0,N.f), 1)
        f.con <- rbind(f.out, f.in.f, f.int)
        optimum <- lp("min", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)))
      }
      ## print(optimum)
      theta[j] <- 1/optimum$solution[1]
    }else if(orientation == "OUT"){
      if(has.v == 0){
        if(slack=='NO'){
          f.out <- cbind(outputs[,j], -outputs)
        }else{
          f.out <- cbind(outputs[,j], -outputs, diag(M), matrix(0, nrow=M, ncol=N.f))
        }
        if(convex=='YES'){
          f.rhs <- c(rep(0,M) , f.inputs[,j], rep(0,J), 1)
          f.con <- rbind(f.out, f.in.f, f.neg, f.int)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs)
        }else{
          f.rhs <- c(rep(0,M) , f.inputs[,j], 1)
          f.con <- rbind(f.out, f.in.f, f.int)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)))
        }
        ## print(optimum)
        theta[j] <- 1/optimum$solution[1]
        z[j,] <- optimum$solution[2:(1+J)]
        if(slack=="YES"){
          s[j,] <- optimum$solution[(2+J):(1+J+M)]
          e[j,] <- optimum$solution[(2+J+M):(1+J+M+N.f)]          
        }
      }else{
        if(slack=='NO'){
          f.out <- cbind(outputs[, j], -outputs, matrix(0, nrow=M, ncol=N.v))
        }else{
          f.out <- cbind(outputs[, j], -outputs,
                         matrix(0, nrow=M, ncol=N.v),
                         diag(M), matrix(0, nrow=M, ncol=N.f))
        }
        f.in.v[, (1+J+1):(1+J+N.v)] <- -diag(N.v)*v.inputs[, j]
        if(convex=='YES'){
          f.con <- rbind(f.out, f.in.f, f.in.v, f.neg, f.int)
          f.rhs <- c(rep(0,M) , f.inputs[,j], rep(0, N.v+N.v+J), 1)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs)
        }else{
          f.con <- rbind(f.out, f.in.f, f.in.v, f.neg, f.int)
          f.rhs <- c(rep(0,M) , f.inputs[,j], rep(0, N.v+N.v), 1)
          optimum <- lp("max", f.obj, f.con, f.dir, f.rhs, binary.vec=(2:(1+J)))
        }
        ## print(optimum)        
        theta[j] <- 1/optimum$solution[1]
        lambda[j,] <- optimum$solution[(1+J+1):(1+J+N.v)]
        z[j,] <- optimum$solution[2:(1+J)]
        if(slack=="YES"){
          s[j,] <- optimum$solution[(2+J+N.v):(1+J+N.v+M)]
          e[j,] <- optimum$solution[(2+J+N.v+M):(1+J+N.v+M+N.f)]          
        }
      }
    }
  }
  ##return(optimum)  
  rm(optimum)
  ## print(f.obj)
  ## print(f.con)
  ## print(f.dir)
  ## print(f.rhs)
  if(orientation == 'OUT'){
    if(report.z == 'NO' & has.v == 0 & slack=='NO'){
      list(theta=theta)
    }else if(report.z == 'YES' & has.v == 0 & slack=='NO'){
      list(theta=theta, z=z)
    }else if(report.z == 'NO' & has.v == 1 & slack=='NO'){
      list(theta=theta, lambda=lambda)
    }else if(report.z == 'YES' & has.v == 1 & slack=='NO'){
      list(theta=theta, z=z, lambda=lambda)
    }else if(report.z == 'YES' & has.v == 0 & slack=='YES'){
      list(theta=theta, z=z, s=s, e=e)
    }else if(report.z == 'YES' & has.v == 1 & slack=='YES'){
      list(theta=theta, z=z, lambda=lambda, s=s, e=e)
    }else if(report.z == 'NO' & has.v == 1 & slack=='YES'){
      list(theta=theta, lambda=lambda, s=s, e=e)
    }else if(report.z == 'NO' & has.v == 0 & slack=='YES'){
      list(theta=theta, s=s, e=s)
    }
  }else if(orientation == 'IN'){
    list(theta=theta)
  }
}
