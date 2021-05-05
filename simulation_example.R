################################################################################
##  <May 05, 2021> @ Shikun Wang
##  ---------------------------------------------------------------------------
##  A simulation example for manuscript "
##  Statistical Model of Longitudinal Medical Cost Trajectory"
##    
##    
################################################################################
rm(list=ls())
source("Algorithm.R")
lambda <- 1e-4;lambda2= 1e-2
AR = sn = zn = mn = 0
order = 2;dir = 2
LTS.group <- T
file.name = "test"
censor = 30;constant.variance=T
n = 500;sigma.e = 2;sigma.tau = 3;rho = .5;rho1 = .2
valid.set <- gen.valid.set(LTS.group = LTS.group)

dat <- read.csv('simulated_data.csv')
param.Naive <- NaiveMethod(dat,LTS.group = LTS.group)
param.stage1 <- Stage1Method(dat,LTS.group = LTS.group)
param.stage2 <- Stage2Method(dat,param.stage1,LTS.group = LTS.group)

dat.test <- gen.valid.set(LTS.group = LTS.group)
Dtest <- DesignMatrix(dat.test)
traj <- Dtest %*% param.stage2$theta
dat.test$Yest <- traj2
Dtest.LTS <- DesignMatrix(dat.test[dat.test$death>tau,],ks=F)
traj.LTS <- Dtest.LTS %*% param.stage2$eta
dat.test$Yest[dat.test$death>tau] <- traj.LTS
plot(NULL,
     xlim = c(0,tau),# range(dat$time),
     ylim = c(2,12),
     xlab = "quarter",ylab = "Y")
points(dat.test$time[dat.test$death>tau],dat.test$Yest[dat.test$death>tau],col = 'grey')
for(i in 1:length(s.list)){
  points(dat.test$time[dat.test$death==s.list[i]],dat.test$Yest[dat.test$death==s.list[i]],
         col = i)
}


