# two-stage-cost-trajectory
Code in paper: Wang, S. et al. (2021). Statistical Modeling of Longitudinal Medical Cost Trajectory: Renal Cell Cancer Care Cost Analyses

# Example
```
rm(list=ls())
source("Algorithm.R")
lambda = 1e-4;lambda2 = 1e-2
AR = sn = zn = mn = 0
order = 2;dir = 2
extra.group = T
file.name = "test"
censor = 30;constant.variance=T
n = 500;sigma.e = 2;sigma.tau = 3;rho = .5;rho1 = .2
valid.set <- gen.valid.set(extra.group = extra.group)
dat <- gen.data(n = n,
                sigma.e = sigma.e,sigma.tau = sigma.tau,rho = rho,rho1 = rho1,
                censor = censor,
                seed = 123,
                constant.variance = constant.variance)
plot.data(dat)
param.Naive <- NaiveMethod(dat,extra.group = extra.group)
param.stage1 <- Stage1Method(dat,extra.group = extra.group)
param.stage2 <- Stage2Method(dat,param.stage1,extra.group = extra.group)
```
