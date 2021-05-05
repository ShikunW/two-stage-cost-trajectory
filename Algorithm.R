################################################################################
##  <May 05, 2021> @ Shikun Wang
##  ---------------------------------------------------------------------------
##  Functions:
##    DesignMatrix(dat.unc,tau)
##    Convergence(new, old)
##    surv.calc(dat,s.list)
##    weight.calc(dat,para,survp,s.list)
##    FirstStep(dat,para,survp)
##    SecondStep(newdat,para)
##    Stage1Method(dat.unc)
##    Stage2Method(dat, para,knots.t),knots.s)
##    
##    
################################################################################

library(survival)
library(nlme)
library(mvtnorm)
# library(sn)
s.list = seq(4,32,4)
#s.list = 1:40
para.true <- c(2,.5,3,.2)
order <- 1
tau <- max(s.list)
lambda <- 1e-4
set.seed(123)
seed.lst <- sample(10000:99999,2000,replace = F)
gen.data <- function(n = 1000,censor = 30,
                     sigma.e = 2,sigma.tau = 3,rho = .5,rho1 = .2,
                     seed = 123,
                     constant.variance = TRUE){
  #sigma.e = 2;sigma.tau = 3;rho = .5;rho1 = .2
  set.seed(seed)
  if(censor == 30){
    rate = 0.025
  }else if(censor == 50){
    rate = 0.06
  }else{
    print("censor rate not valid!");break}
  
  dat.death <- cbind(ceiling(rexp(n,rate=.05)/4)*4,
                     999)#ceiling(rexp(n,rate)/4)*4)
  dat.surv <- apply(dat.death,1,function(x)min(min(x),tau))
  dat.tau <- 1*(dat.death[,1] > tau)
  dat.delta <- apply(dat.death,1,function(x)1*(x[1] <= x[2] & x[1] <= tau))
  
  dat.tmp <- data.frame(id = rep(1:n,dat.surv),
                        surv = rep(dat.surv,dat.surv),
                        death = rep(dat.death[,1],dat.surv),
                        time = unlist(sapply(dat.surv,function(x)1:x)),
                        delta = rep(dat.delta,dat.surv),
                        tau = rep(dat.tau,dat.surv))
  dat.tmp$Y <- apply(dat.tmp,1,function(x){
    if(x["death"] <= tau){
      if(x["time"]<=4){
        1/3 *(x["time"] - 4)^2 + 7
      }else if(x["time"]<=x["death"]/2){
        7
      }else{
        12 / x["death"]^2 *(x["time"] - x["death"]/2)^2 + 7
      }
    }else{
      if(x["time"]<=6){
        1/5*(x["time"] - 6)^2 + 5
      }else{
        5
      }
    }
  })
  for(id in unique(dat.tmp$id)){
    ni <- sum(dat.tmp$id==id)
    if(dat.tmp[dat.tmp$id == id,"death"][1]<=tau){
      s <- sigma.e
      r <- rho
    }else{
      s <- sigma.tau
      r <- rho1
    }
    Sigma <- matrix(r,ni,ni)
    if(constant.variance){
      diag(Sigma) <- 1 
    }else{
      if(AR == 0){
        #diag(Sigma) <- ((1:ni)^(1/3))
        death <- unique(dat.tmp$death[dat.tmp$id==id])
        diag(Sigma) <- (death/12)^(1/3)
      }else{
        Sigma <- matrix(0,ni,ni)
        for(i in 1:ni){
          for(j in 1:ni){
            Sigma[i,j] <- rho^(abs(i-j))
          }
        }
        diag(Sigma) <- 1
        
      }
    }
    Sigma <- Sigma * s^2
    err <- rmvnorm(1,mean = rep(0,ni),sigma = Sigma)
    if(sn){
      err <- rgamma(ni,10,1) - 10
    }else if(mn){
      err <- rnorm(ni,0,s/2) +4*sapply(1:ni,function(x)rbinom(1,size = 1,prob = .5)) - 2
      #hist(err)
    }
    #hist(rsn(n=1e5, xi=-.7, omega=1, alpha=5, tau=0,  dp=NULL),nclass=50)
    dat.tmp$Y[dat.tmp$id == id] <- dat.tmp$Y[dat.tmp$id == id] + err
    if(zn){
      dat.tmp$Y[dat.tmp$id == id] <- dat.tmp$Y[dat.tmp$id == id] * sapply(1:ni,function(x)rbinom(1,size = 1,prob = .9))
    }
  }
  # randomly set 20% cost > 8 to be missing
  # idx <- which(dat.tmp$Y>8)
  # delidx <- sample(idx, floor(length(idx)/20),replace=F)
  # dat.tmp <- dat.tmp[-delidx,]
  
  # cost at first 12 months<35, 30% censoring; else, 50% censoring
  dat1=dat.tmp[dat.tmp$time<= 4,]
  dat.initial <- aggregate(dat1$Y,list(dat1$id),sum)
  dat.death <- cbind(dat.death[,1],
                     (dat.initial[,2]<35) * ceiling(rexp(n,0.025)/4)*4+ 
                     (dat.initial[,2]>=35) * ceiling(rexp(n,0.06)/4)*4)#ceiling(rexp(n,rate)/4)*4)
  dat.surv <- apply(dat.death,1,function(x)min(min(x),tau))
  dat.tau <- 1*(dat.death[,1] > tau)
  dat.delta <- apply(dat.death,1,function(x)1*(x[1] <= x[2] & x[1] <= tau))
  
  dat.tmp$delta <- dat.delta[match(dat.tmp$id,1:n)]
  dat.tmp$surv <- dat.surv[match(dat.tmp$id,1:n)]
  dat.tmp <- dat.tmp[dat.tmp$time <= dat.tmp$surv,]
  dat.tmp <- dat.tmp[,c("id","surv","death","time","delta","Y")]
  
  return(dat.tmp)
}
plot.data <- function(dat,para = NA,scatter = T,band = T){
  # Plot scatter plot of uncensored data, colored by death group
  #
  # Args:
  #   dat: Data from gen.data() 
  # Returns:
  #   scatter plot
  if(is.na(para)[1]){
    plot(NULL,
         xlim = c(0,tau),# range(dat$time),
         ylim = range(dat$Y),
         xlab = "quarter",ylab = "Y")
    points(dat$time[dat$death>max(s.list)],dat$Y[dat$death>max(s.list)],col = 'grey')
    for(i in 1:length(s.list)){
      points(dat$time[dat$death==s.list[i]],dat$Y[dat$death==s.list[i]],
             col = i)
    }
  }else{
    plot(NULL,
         xlim = c(0,tau),# range(dat$time),
         ylim = c(5,10),# range(dat$Y),
         xlab = "quarter",ylab = "Y")
    if(scatter){
      points(dat$time[dat$death>max(s.list)],dat$Y[dat$death>max(s.list)],col = 'grey')
      for(i in 1:length(s.list)){
        points(dat$time[dat$death==s.list[i]],dat$Y[dat$death==s.list[i]],
               col = "grey")
      }
    }
    valid.set <- para$valid.set
    theta.idx <- 1:(sum(s.list+1))
    Dmat1 <- DesignMatrix(valid.set[valid.set$surv<=max(s.list),])
    p <- ncol(Dmat1)
    para1 <- para$para1
    para2 <- para$para2
    para3 <- para$para3
    
    traj1 <- Dmat1 %*% para1[theta.idx,]
    traj2 <- Dmat1 %*% para2[theta.idx,]
    traj3 <- Dmat1 %*% para3[theta.idx,]
    # 
    if(extra.group){ # with eta
      eta.idx <- max(theta.idx) + 2 + (1: (tau+1))
      Dmat2 <- DesignMatrix(valid.set[valid.set$surv>max(s.list),],
                            ks = F)
      traj1e <- Dmat2 %*% para1[eta.idx,]
      traj2e <- Dmat2 %*% para2[eta.idx,]
      traj3e <- Dmat2 %*% para3[eta.idx,]
      traj1 <- rbind(traj1,traj1e)
      traj2 <- rbind(traj2,traj2e)
      traj3 <- rbind(traj3,traj3e)
    }
    traj1 <- data.frame(time = valid.set$time,
                        surv = valid.set$surv,
                        mean = apply(traj1,1,mean),
                        lb = apply(traj1,1,function(x)quantile(x,.025)),
                        ub = apply(traj1,1,function(x)quantile(x,.975)))
    traj2 <- data.frame(time = valid.set$time,
                        surv = valid.set$surv,
                        mean = apply(traj2,1,mean),
                        lb = apply(traj2,1,function(x)quantile(x,.025)),
                        ub = apply(traj2,1,function(x)quantile(x,.975)))
    traj3 <- data.frame(time = valid.set$time,
                        surv = valid.set$surv,
                        mean = apply(traj3,1,mean),
                        lb = apply(traj3,1,function(x)quantile(x,.025)),
                        ub = apply(traj3,1,function(x)quantile(x,.975)))
    for(i in 1:length(s.list)){
      if(1){ # i %% 2 ==1
        i=i+1;
        sub1 <- traj1[traj1$surv==s.list[i],]
        sub2 <- traj2[traj2$surv==s.list[i],]
        sub3 <- traj3[traj3$surv==s.list[i],]
        if(band){
          polygon(c(sub1$time,rev(sub1$time)),c(sub1$lb,rev(sub1$ub)),col=rgb(0, 0, 0,.1),border = NA)
          polygon(c(sub2$time,rev(sub2$time)),c(sub2$lb,rev(sub2$ub)),col=rgb(1, 0, 0,.1),border = NA)
          polygon(c(sub3$time,rev(sub3$time)),c(sub3$lb,rev(sub3$ub)),col=rgb(0, 1, 0,.1),border = NA)
        }
        lines(mean ~ time,data = sub1,col = 1)
        lines(mean ~ time,data = sub2,col = 2,lty = 2)
        lines(mean ~ time,data = sub3,col = 3,lty = 3)
      }
    }
    sub1 <- traj1[traj1$surv>max(s.list),]
    sub2 <- traj2[traj2$surv>max(s.list),]
    sub3 <- traj3[traj3$surv>max(s.list),]
    if(band){
      polygon(c(sub1$time,rev(sub1$time)),c(sub1$lb,rev(sub1$ub)),col=rgb(0, 0, 0,.1),border = NA)
      polygon(c(sub2$time,rev(sub2$time)),c(sub2$lb,rev(sub2$ub)),col=rgb(1, 0, 0,.1),border = NA)
      polygon(c(sub3$time,rev(sub3$time)),c(sub3$lb,rev(sub3$ub)),col=rgb(0, 1, 0,.1),border = NA)
    }
    lines(mean ~ time,data = sub1,col = 1)
    lines(mean ~ time,data = sub2,col = 2,lty = 2)
    lines(mean ~ time,data = sub3,col = 3,lty = 3)
    legend("topleft",c("1-Stage", "2-Stage", "EM"),
           col = 1:3,lty = 1:3,bty = "n",bg = "white")
    
  }
}
gen.valid.set <- function(extra.group){
  # generate valid set
  #
  # Args:
  #   extra.group: indicator for including the LTS group
  # Returns:
  #   The design matrix of valid set dat.valid.design
  dat.tmp <- data.frame(surv = unlist(sapply(s.list,function(x) rep(x,x))),
                        death = unlist(sapply(s.list,function(x) rep(x,x))),
                        time = unlist(sapply(s.list,function(x)1:x)),
                        delta = 1)
  dat.tmp$Y <- apply(dat.tmp,1,function(x){
    if(x["time"]<=4){
      1/3 *(x["time"] - 4)^2 + 7
    }else if(x["time"]<=x["death"]/2){
      7
    }else{
      12 / x["death"]^2 *(x["time"] - x["death"]/2)^2 + 7
    }
  })
  # extra group
  if(extra.group){
    dat.tmp1 <- data.frame(surv = rep(tau+1,tau),
                           death = rep(tau+1,tau),
                           time = 1:tau,
                           delta = 1)
    dat.tmp1$Y <- apply(dat.tmp1,1,function(x){
      if(x["time"]<=6){
        1/5*(x["time"] - 6)^2 + 5
      }else{
        5
      }
    })
    
    dat.tmp <- rbind(dat.tmp,dat.tmp1)
  }
  
  plot(dat.tmp[,c("time","Y")]) # plot valid dataset
  dat.tmp
}

DesignMatrix <- function(dat.unc,
                         ks = TRUE){
  # Computes the Design Matrix
  #
  # Args:
  #   dat.unc: dat.unca with only uncensored patients
  #   ks: indicator of not creating for LTS gsoup
  # Returns:
  #   The design matrix D
  if(ks){
    tmp <- data.frame(death=rep(s.list,s.list),time = unlist(sapply(s.list,function(x)1:x)))
    D <- matrix(apply(tmp,1,function(x)1*(x['death'] == dat.unc$surv) * (x['time'] == dat.unc$time)),
                nrow(dat.unc),sum(s.list))
  }else{
    D <- matrix(sapply(1:tau,function(x) 1*(x == dat.unc$time)),ncol=tau)
  }
  D
}

PenaltyMatrix <- function(ks = TRUE){
  # Computes the Penalty Matrix
  #
  # Args:
  #   ks: indicator of not creating for LTS group
  # Returns:
  #   The penalty matrix Omega
  if(ks){
    tmp <- data.frame(death=rep(s.list,s.list),time = unlist(sapply(s.list,function(x)1:x)))
    if(order == 1){
      Q <- matrix(0,sum(s.list),sum(s.list))
      for(i in 1:length(s.list)){
        if(s.list[i] >= 2){
          for(j in 1:(s.list[i]-1)){
            Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == j) -1 * (tmp$death == s.list[i])*(tmp$time == j+1),ncol=1)
            Di <- Di %*% t(Di)
            Q <- Q + Di
          }
        }
      }
      if (dir == 2){
        Q2 <- matrix(0,sum(s.list),sum(s.list))
        for(j in 1:tau){
          if(which(s.list >= j)[1] <= (length(s.list) - 1)){
            for(i in which(s.list >= j)[1]:(length(s.list) - 1)){
              Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == j) -1 * (tmp$death == s.list[i+1])*(tmp$time == j),ncol=1)
              Di <- Di %*% t(Di)
              Q2 <- Q2 + Di
            }
          }
        }
        Q <- Q+Q2
        Q3 <- matrix(0,sum(s.list),sum(s.list))
        for(i in 1:(length(s.list)-1)){
          Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == s.list[i])- 1 * (tmp$death == s.list[i+1])*(tmp$time == s.list[i+1]),ncol=1)
          Di <- Di %*% t(Di)
          Q3 <- Q3 + Di
        }
      }
      Q
    }else if(order == 2){
      Q <- matrix(0,sum(s.list),sum(s.list))
      for(i in 1:length(s.list)){
        if(s.list[i] >= 3){
          for(j in 1:(s.list[i]-2)){
            Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == j) + 1*(tmp$death == s.list[i])*(tmp$time == j+2) -2 * (tmp$death == s.list[i])*(tmp$time == j+1),ncol=1)
            Di <- Di %*% t(Di)
            Q <- Q + Di
          }
        }
      }
      if (dir == 2){
        Q2 <- matrix(0,sum(s.list),sum(s.list))
        
        for(j in 1:tau){
          if(which(s.list >= j)[1] <= (length(s.list) - 2)){
            for(i in which(s.list >= j)[1]:(length(s.list) - 2)){
              Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == j) + 1*(tmp$death == s.list[i+2])*(tmp$time == j) -2 * (tmp$death == s.list[i+1])*(tmp$time == j),ncol=1)
              Di <- Di %*% t(Di)
              Q2 <- Q2 + Di
            }
          }
        }
        Q <- Q+Q2
        Q3 <- matrix(0,sum(s.list),sum(s.list))
        for(i in 1:(length(s.list)-1)){
          Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == s.list[i]) -1 * (tmp$death == s.list[i+1])*(tmp$time == s.list[i+1]),ncol=1)
          Di <- Di %*% t(Di)
          Q3 <- Q3 + Di
          Di <- matrix(1*(tmp$death == s.list[i])*(tmp$time == 1) -1 * (tmp$death == s.list[i+1])*(tmp$time == 1),ncol=1)
          Di <- Di %*% t(Di)
          Q3 <- Q3 + Di
        }
        Q <- Q + Q3
      }
      Q
    }
    
  }else{
    Q <- matrix(0,tau,tau)
    if(order == 1){
      for(j in 1:(tau-1)){
        Di <- matrix(0,tau,ncol=1)
        Di[j] <- 1;Di[j+1] <- -1
        Di <- Di %*% t(Di)
        Q <- Q + Di
      }
    }else if (order == 2){
      for(j in 1:(tau-2)){
        Di <- matrix(0,tau,ncol=1)
        Di[j] <- 1;Di[j+1]<- -2;Di[j+2] <- 1
        Di <- Di %*% t(Di)
        Q <- Q + Di
      }
    }
    Q
  }
}

Convergence <- function(new, old){
  # Decide if new and old parameters are close
  #
  # Args:
  #   new: new paramters
  #   old: old paramters
  # Returns:
  #   The convergence TRUE or FALSE
  p <- 4
  n.predictor <- length(new)
  ind <- rep(0, n.predictor)
  for (i in 1:n.predictor) {
    if (new[i] < abs(1e-4)){
      ind[i] <- as.numeric( abs(new[i]-old[i]) < 10^(-p) )
    } else{
      ind[i] <- as.numeric( abs(new[i]-old[i]) / (abs(old[i])+10^(-5)) < 10^(-p) )
    } 
  }
  det <- FALSE
  # cat("converge number = ",  sum(ind), "\n")
  if(sum(ind) == n.predictor ) det <- TRUE
  det
}

surv.calc <- function(dat,extra.group){
  # calculate the KM estimate of death probability in each death group
  #
  # Args:
  #   dat: data
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   vector of death probability
  dat.surv <- dat[!duplicated(dat$id),]
  tt <- summary(survfit(Surv(time = dat.surv$surv, event = dat.surv$delta)~1))
  match <- match(tt$time, s.list)
  tmp.f <- tt$n.event / tt$n.risk * c(1,tt$surv[1:(length(tt$surv)-1)] )
  surv.f <- rep(0, length(s.list))
  surv.f[match] <- tmp.f
  if(extra.group){
    c(surv.f,tail(tt$surv,1))
  }else{
    surv.f
  }
}

weight.calc <- function(dat,survp,para, extra.group){
  # calculate the weight estimate of patients belongs to each death group
  #
  # Args:
  #   dat: data
  #   survp: survival probability
  #   para: parameters
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   matrix of death probability cols are death group, rows are patient id
  
  id.list <- unique(dat$id)
  newdat <- NULL
  for(i in id.list){
    if(!extra.group){wt = rep(0,length(s.list))}else{wt = rep(0,length(s.list)+1)}
    dati <- subset(dat,id == i)
    dti <- dati$delta[1]
    for(j in 1:length(s.list)){
      s <- s.list[j]
      if(dati$delta[1]){
        wt[j] <- 1*(dati$surv[1]==s)
      }else if(dati$surv[1]>=s){
        wt[j] <- 0
      }else{
        ni <- nrow(dati)
        Yi <- matrix(dati$Y,ncol=1)
        dati$surv <- s
        Dmatis <- DesignMatrix(dati,
                               ks = T)
        mu.s <- Dmatis %*% para$theta
        Sigma.s <- matrix(para$rho,nrow = ni, ncol = ni)
        diag(Sigma.s) <- 1
        Sigma.s <- para$sigma^2 * Sigma.s
        p <- 1/sqrt(det(Sigma.s)) * exp(-0.5 * t(Yi - mu.s) %*% solve(Sigma.s) %*% (Yi - mu.s) )  
        wt[j] <- p * survp[j]
      }
    }
    if(dati$delta[1] == 0  & extra.group){
      ni <- nrow(dati)
      Yi <- matrix(dati$Y,ncol=1)
      dati$surv <- tau
      Dmatis <- DesignMatrix(dati,
                             ks = F)
      mu.s <- Dmatis %*% para$eta
      Sigma.s <- matrix(para$rho1,nrow = ni, ncol = ni)
      diag(Sigma.s) <- 1
      Sigma.s <- para$omega^2 * Sigma.s
      p <- 1/sqrt(det(Sigma.s)) * exp(-0.5 * t(Yi - mu.s) %*% solve(Sigma.s) %*% (Yi - mu.s))
      wt[length(wt)] <- p * survp[length(survp)]
    }
    
    if(sum(wt) > 0) wt <- round(wt / sum(wt),3)
    for(j in 1:length(s.list)){
      if(wt[j] > 0){
        dati <- subset(dat,id == i)
        dati$id <- paste0(dati$id,'_',s.list[j])
        dati$surv <- s.list[j]
        dati$delta <- dti
        dati$weight <- wt[j]
        dati$tau <- 0
        newdat <- rbind(newdat, dati)
      }
    }
    if(wt[length(wt)] & extra.group){
      dati <- subset(dat,id == i)
      dati$id <- paste0(dati$id,'_tau')
      dati$surv <- tau
      dati$delta <- 0
      dati$weight <- wt[length(wt)]
      dati$tau <- 1
      newdat <- rbind(newdat, dati)
    }
  }
  newdat
}

para.calc <- function(dat,para = NULL,ks,lbda){
  if(!("weight" %in% names(dat))) dat$weight = 1
  n <- nrow(dat)
  #
  Omega <- PenaltyMatrix(ks = ks)
  p <- nrow(Omega)
  # 
  # initial beta from independence assumption
  if(is.null(para)){
    Dmat <- DesignMatrix(dat,ks = ks)
    Y <- matrix(dat$Y,ncol = 1)
    beta.new <- solve(t(Dmat) %*% Dmat / n + lbda * Omega) %*% (t(Dmat) %*% Y / n)
    rho.new <- 0
    err.dat <- data.frame(id = dat$id,
                          err = Y - Dmat %*% beta.new)
    sigma.new <- 1
  }else{
    beta.new <- para$beta
    rho.new <- para$rho
    sigma.new <- para$sigma
    Dmat <- DesignMatrix(dat,ks = ks)
    Y <- matrix(dat$Y,ncol = 1)
    beta.new <- para$beta
    rho.new <- para$rho
    err.dat <- data.frame(id = dat$id,
                          err = Y - Dmat %*% beta.new)
  }
  itercount <- 0
  repeat{
    beta.old <- beta.new
    rho.old <- rho.new
    sigma.old <- sigma.new
    lm <- gls(err ~ 1, correlation=corCompSymm(form=~1|id), data=err.dat)
    sigma.new <- lm$sigma
    rho.new <- coef(lm$modelStruct$corStruct, unconstrained=FALSE)
    tmp1 <- matrix(0,p,p)
    tmp2 <- matrix(0,p,1)
    
    for(i in unique(dat$id)){
      dati <- subset(dat,id == i)
      ni <- nrow(dati)
      
      Vi <- matrix(rho.new,ni,ni)
      diag(Vi) <- 1
      Vi <- Vi * sigma.new^2
      Vi.i <- solve(Vi)
      Dmati <- DesignMatrix(dati,ks = ks)
      Yi <- matrix(dati$Y,ncol=1)
      tmp1 <- tmp1 + t(Dmati) %*% Vi.i %*% Dmati * dati$weight[1]
      tmp2 <- tmp2 + t(Dmati) %*% Vi.i %*% Yi * dati$weight[1]
    }
    
    beta.new <- solve(tmp1 / n + lbda * Omega) %*% (tmp2 / n)
    old <- c(beta.old,sigma.old,rho.old)
    new <- c(beta.new,sigma.new,rho.new)
    #print(paste0("Itercount = ",itercount))
    if(Convergence(new,old)){
      break
    }else if(itercount == 100){
      print("Reach max iteration.")
      break
    }
    err.dat <- data.frame(id = dat$id,
                          err = Y - Dmat %*% beta.new)
    itercount <- itercount + 1
    # print(head(beta.new))
  }
  return(list(beta = beta.new,sigma = sigma.new,rho = rho.new))
}

FirstStep <- function(dat,survp,para,extra.group){
  #
  # Args:
  #   dat: data
  #   para: parameter
  #   survp: survival probability
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   matrix of death probability cols are death group, rows are patient id
  # Error handling
  weight.calc(dat,survp,para,extra.group)
}

SecondStep <- function(dat,para,extra.group){
  # Compute the M-step: weights for each patients for each death group
  #
  # Args:
  #   newdat: data
  #   para: parameter
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   parameters
  # Error handling
  newdat <- subset(dat,tau == 0)    
  dat.surv <- dat[!duplicated(dat$id),]
  para.init <- list(beta = para$theta, sigma = para$sigma, rho = para$rho)
  para.nonextra <- para.calc(newdat,para.init,ks = T,lbda = lambda)
  if(!extra.group){
    return(list(theta = para.nonextra$beta,
                sigma = para.nonextra$sigma,
                rho = para.nonextra$rho))
  }else{
    dat.tau <- subset(dat, delta==0 & tau == 1)
    para.init <- list(beta = para$eta, sigma = para$omega, rho = para$rho1)
    para.extra <- para.calc(dat.tau,para.init,ks = F,lbda = lambda)
    return(list(theta = para.nonextra$beta,sigma = para.nonextra$sigma,rho = para.nonextra$rho,
                eta = para.extra$beta,omega = para.extra$sigma,rho1 = para.extra$rho))
  }
}

Stage1Method <- function(dat,extra.group){
  # Computes the stage 1 estimate parameters
  #
  # Args:
  #   dat.unc: Data with only uncensored patients
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   The list of theta, sigma and rho.
  # Error handling
  dat.unc <- subset(dat,delta == 1)
  dat.surv <- dat[!duplicated(dat$id),]
  para.nonextra <- para.calc(dat.unc,para = NULL, ks = T,lbda = lambda)
  if(!extra.group){
    return(list(theta = para.nonextra$beta,
                sigma = para.nonextra$sigma,
                rho = para.nonextra$rho))
  }else{
    # extra group
    dat.tau <- subset(dat, delta==0 & surv == tau)
    para.extra <- para.calc(dat.tau, para = NULL,ks = F,lbda = lambda)
    return(list(theta = para.nonextra$beta,sigma = para.nonextra$sigma,rho = para.nonextra$rho,
                eta = para.extra$beta,omega = para.extra$sigma,rho1 = para.extra$rho))
  }
}

Stage2Method <- function(dat,para,extra.group){
  # Computes the stage 2 estimate parameters
  #
  # Args:
  #   dat: data with all patients
  #   para: initial parameters
  #   extra.group: indicator for creating LTS group
  # Returns:
  #   The list of theta, sigma and rho; 
  #   if extra.group is true, also return eta, omega and rho1.
  # Error handling
  survp <- surv.calc(dat,extra.group)
  newdat <- FirstStep(dat,survp,para,extra.group)
  para.new <- SecondStep(newdat,para,extra.group)
  para.new
}


NaiveMethod <- function(dat,extra.group){
  valid.set <- gen.valid.set(extra.group = extra.group)
  start_time <- Sys.time()
  dat.unc <- subset(dat,delta == 1)
  tmp <- data.frame(
    death = rep(s.list,s.list),
    time = unlist(sapply(s.list,function(x)1:x)))
  for(i in 1:nrow(tmp)){
    dati <- subset(dat.unc,surv == tmp$death[i] & time == tmp$time[i])
    tmp$Y[i] <- mean(dati$Y)
  }
  if(!extra.group){
    return(tmp$Y)
  }else{
    
    dat.tau <- subset(dat, delta==0 & surv == tau)
    tmp.extra <- data.frame(
      death = rep(tau,tau),
      time = 1:tau)
    for(i in 1:nrow(tmp.extra)){
      dati <- subset(dat,surv == tmp.extra$death[i] & delta == 0 & time == tmp.extra$time[i])
      tmp.extra$Y[i] <- mean(dati$Y)
    }
    result = (rbind(tmp,tmp.extra)$Y)
  }
  end_time <- Sys.time()
  
  result 
  
}


#extra.group = F
#parametric = F
lambda.lst <- 10^((-7):(-1))
choose.lambda <- function(lambda.lst,extra.group){
  set.seed(123)
  dat <- gen.data(n = 400,n.tau = 100,censor = 30,
                  seed = 123,
                  constant.variance = TRUE)
  n = max(dat$id)
  idx <- sample(1:n,size = n, replace = F)
  mse.vec <- matrix(0,3,length(lambda.lst))
  time.vec <- matrix(0,3,length(lambda.lst))
  for(i in 1:length(lambda.lst)){
    lambda <- lambda.lst[i]
    mse <- matrix(0,3,10)
    time <- matrix(0,3,10)
    for(fold in 1:10){
      dat.test <- subset(dat,id %in% idx[(1:(n/10))+n/10*(fold-1)])
      dat.test <- subset(dat.test,delta == 1)
      dat.train <- subset(dat,!(id %in% idx[(1:(n/10))+n/10*(fold-1)]))
      dat.unc <- subset(dat.train,delta == 1)
      start_time <- Sys.time()
      param.stage1 <- Stage1Method(dat.train,extra.group)
      end_time <- Sys.time()
      time[1,fold] <- end_time - start_time
      start_time <- Sys.time()
      param.stage2 <- Stage2Method(dat.train,param.stage1,extra.group)
      end_time <- Sys.time()
      time[2,fold] <- end_time - start_time
      # start_time <- Sys.time()
      # param.EM <- EMMethod(dat.train,param.stage2,extra.group)
      # end_time <- Sys.time()
      # time[3,fold] <- end_time - start_time
      Dtest <- DesignMatrix(dat.test)
      traj1 <- Dtest %*% param.stage1$theta
      traj2 <- Dtest %*% param.stage2$theta
      # traj3 <- Dtest %*% param.EM$theta
      mse[1,fold] <- mean((dat.test$Y - traj1)^2)
      mse[2,fold] <- mean((dat.test$Y - traj2)^2)
      # mse[3,fold] <- mean((dat.test$Y - traj3)^2)
      print(fold)
      gc()
    }
    mse.vec[,i] <- apply(mse,1,mean)
    time.vec[,i] <- apply(time,1,mean)
    print(mse.vec)
  }
}
