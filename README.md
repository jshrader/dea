dea
===

The main code is in two files. The first, `dea_package.R`, contains everything you need to run DEA and industry model estimates. The second file contains examples for code testing. 

I based my DEA syntax on the package FEAR (http://www.clemson.edu/economics/faculty/wilson/Software/FEAR/fear.html), so if you use that package, things should be familiar. 

The papers in the documentation folder are: the one on which I based my coding of the industry model (Kerstens et al 2005), the paper on which I based the DEA code (Reid et al 2003), and the paper Dale and I wrote on the EPO (Shrader and Squires 2013).

I don't have an official help file, so here is a quick run-down of the four commands and syntax you need:

### Dependencies
You will need the following R packages:
* linprog
* foreach
* doParallel (My code currently requires this one, but that might be a mistake)

In addition, to run the `dea_sandbox.R` you will need FEAR so that results can be compared.

### `f.dea`: 
This command calculates efficient frontiers from an input or output oriented perspective for fixed or variable (or both) inputs. To calculate technical efficiency, include your variable inputs as if they were fixed. To calculate capacity utilization, include your variable inputs as variable inputs. See my and Dale's paper for more details.

Syntax: 

`output <- f.dea(f.inputs, <v.inputs>, outputs, tech='V', orientation='OUT', report.z='NO', slack='NO', convex='YES')`

Variables:

`f.inputs` = an (n by k_f) matrix of fixed inputs (e.g. tonnage) where n is the number of firms and k is the number of fixed inputs.

`v.inputs` (optional) = an (n by k_v) matrix of variable inputs (e.g. days).

`outputs` = an (n by j) matrix of outputs where j is the number of outputs.

`tech` = 'V', 'N', or 'C' for variable, non-increasing, or constant technology. To my mind, there is no reason in practice to run anything other than 'V'. If you want to run 'C', use FEAR because my code will probably break.

`orientation` = 'OUT' or 'IN' for output or input orientation. Output orientation will maximize outputs while input orientation minimizes inputs. To run the industry model, this must be 'OUT'.

`report.z` = Set this to 'YES' for all runs! I don't know why I was silly and set the default to 'NO'

`slack`, and `convex` should be left at their defaults for standard analyses.

`output` = a list containing *theta* and *lambda*. *Theta* gives the relative efficiency of the vessel (on a scale between 0 and 1) and *lambda* gives the weights on each efficient vessel that the vessel in question should place in order to reach the frontier. Another way of saying this is that an inefficient vessel should strive to be a convex combination of some efficient vessels, and *lambda* tells you the weights in that convex combination. An efficient vessel should always have a lambda equal to 1 for themselves. 

Example:
```R
# Specify fixed inputs
x.m.f <- matrix(c(10, 5, 3, 4, 5), ncol=1, nrow=5)
# Specify variable inputs
x.m.v <- matrix(c(2, 50, 3, 2, 2, 1, 3, 5, 2, 1), ncol=2, nrow=5)
# Specify outputs
y.m <- matrix(c(5, 0, 6, 1, 1,
                1, 4, 2, 0, 0), ncol=2, nrow=5)
# Run 
cu <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v,
            outputs=y.m,
            tech="V", orientation="OUT", report.z="NO", slack="NO", convex="YES")
```


### `f.jim`:
This command calculates the Johansen industry model based optimal input values for the whole fleet. You need to run `f.dea`, `industry.observed`, and `optimal.in.out` before running `f.jim`. Almost all of the inputs are provided by the helper functions.

Syntax: 

`output <- f.jim(f.in.total, <v.in.total>, out.total, f.in.opt, <v.in.opt>, out.opt, tac=0, fd.max=0, ivq=NULL, ivq.var=NULL, convex='YES')`

Variables:

`f.in.total` = provided by industry.observed

`v.in.total` = provided by industry.observed

`out.total` = provided by industry.observed

`f.in.opt` = provided by optimal.in.out

`v.in.opt` = provided by optimal.in.out

`out.opt` = provided by optimal.in.out

`tac` = A total allowable catch for the fishery for each output (if set to 0 or the observed value, then no constraint)

`fd.max` = Fishing day maximum. To use, you must specify days as your first variable output.

`ivq` = Individual vessel quota

`ivq.var` = Column number of the species on which the IVQ applies (could be all).

`convex` = Leave this alone unless you are really interested in non-convex optimization. If you use this option, it will make your code take much, much longer.

`output` = a list containing *theta*, *W*, and *Xv*. *Theta* is the industry-wide input efficiency. *W* is the scaling that should be applied to each vessel in order to reach industry efficiency. Thus, the sum of *W* is the number of vessels you would want in the optimal fishery. *Xv* is the variable input scaling. I think it is the average of *W*.

Example:
```R
# Call the package
source("path/to/dea_package.R")
# Generate random data
set.seed(2982)
J <- 100
nv <- 1
nf <- 2
nu <- 2
x.m.f <- matrix(runif(nf*J, 0, 50), ncol=nf, nrow=J)
x.m.v <- matrix(runif(nv*J, 0, 50), ncol=nv, nrow=J)
y.m <- matrix(runif(nu*J, 0, 50), ncol=nu, nrow=J)
l <- length(y.m)

# Get the total inputs and outputs for the fishery
act <- industry.observed(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m)

# Calculate the DEA estimates
testf <- f.dea(f.inputs = x.m.f, v.inputs=x.m.v, outputs=y.m, tech="V", orientation="OUT",
                 report.z="YES", slack="NO")

# Calculate optimal inputs and outputs for the fishery based on the DEA step
outf <- optimal.in.out(f.inputs=x.m.f, v.inputs=x.m.v, outputs=y.m, theta=testf$theta, z=testf$z)

# Run the industry model
jim.c <- f.jim(f.in.total=act$x.f.a, v.in.total=act$x.v.a, out.total=act$u.a,
                 f.in.opt=outf$x.f.s, v.in.opt=outf$x.v.s, out.opt=outf$u.s,
                 convex='YES')
```

### `optimal.in.out`: 
A helper function for `f.jim`.

Syntax: 

`output <- optimal.in.out(f.inputs, <v.inputs>, outputs, theta, z, <flow>, <flow.var>)`

Variables:

`f.inputs` = same as for `f.dea`

`v.inputs` (optional) = same as for `f.dea`

`outputs` = same as for `f.dea`

`theta` = Just provide the value from `f.dea` output

`z` = Same as `theta`

`flow` and `flow.var` = Just ignore these

`output` = Look at how the `f.jim` example runs to see how to use the output from this command. It just links `f.dea` with `f.jim` to save a little repetitious coding.


### `industry.observed`: 
A helper function for `f.jim`.

Syntax: 

`output <- industry.observed(f.inputs, <v.inputs>, outputs, <flow>, <flow.var>)`

Variables:

`f.inputs` = same as for `f.dea`

`v.inputs` (optional) = same as for `f.dea`

`outputs` = same as for `f.dea`

`flow` and `flow.var` = Just ignore these

`output` = Look at how the `f.jim` example runs to see how to use the output from this command. It just links `f.dea` with `f.jim` to save a little repetitious coding.
