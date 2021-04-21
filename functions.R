# convert the log gamma(large number)
# eg ln(Gamma(900)), ln899!=sum_1_899 ln x
conv.log.gam <- function(a) {
  b <- 0
  c <- a - 1
  if (round(c, 0) == c) {
    while (c > 0) {
      b <- b + log(c)
      c <- c - 1
    }
  }
  else {
    while (c > 0) {
      b <-  b + log(c)
      c <-  c - 1
    }
    b <-  b + log(gamma(c + 1))
  }
  return(b)
}


# MO model, gamma prior

# MOM estimator as starting pt; x is a single colum observation X1, or X2, or X3
mom.lamda <-  function(x) {
  var(x) / mean(x) - 1
}

# the objective to minimize using nlimb
# as the log likelihood is concave, take the negtive
log.likelihood.neg <-  function(par) {
  #par: rho, lam1,lam2, lam3 for tri-NB
  #need n # of observations
  #need prepare row.sum, eg, 46 observations, each 3 varibles, get 46 row sum number
  #x11,x22,x33, vector w.r.t x1,x2,x3, x11<- x1+1, etc, need this so the ln gamma work
  #can add another paramter in the funct
  #to easy the preparation using the following 3 lines
  
  #row.sum<- as.data.frame(rowSums(mydata));nnrow(mydata)
  #x1=mydata[,1];x2=mydata[,2];x3=mydata[,3];
  #x11=as.data.frame(x1+1);x22=as.data.frame(x2+1);x33=as.data.frame(x3+1)
  
  term1  <-  -sum(apply(row.sum + par[1], 1, conv.log.gam))
  term2  <-
    -sum(x1) * log(par[2]) - sum(x2) * log(par[3]) - sum(x3) * log(par[4])
  term3  <-
    n * conv.log.gam(par[1]) + (n * par[1] + sum(x)) * log(sum(par[-1]) +
                                                             1)
  term4  <-
    sum(apply(x11, 1, conv.log.gam)) + sum(apply(x22, 1, conv.log.gam)) +
    sum(apply(x33, 1, conv.log.gam))
  s  <-  term1 + term2 + term3 + term4
  return(s)
}

# simulation using MO model,
# par: rho, lam1,lam2, lam3 for tri-NB
MO.sim  <-  function(nsim, par) {
  theta  <-  rgamma(n = nsim,
                    shape = par[1],
                    scale = 1)
  lambda  <-  par[-1]
  X  <-  matrix(data = NA,
                nrow = nsim,
                ncol = 3)
  for (i in 1:nsim) {
    for (j in 1:3) {
      X[i, j]  <-  rpois(1, lambda = (lambda[j] * theta[i]))
    }
  }
  return(X)
}

## simulation OOC with delta*lambda using MO model
ooc.sim <-
  function(delta = 1,
           delta1 = 1,
           delta2 = 1,
           delta3 = 1,
           delta.rho = 1,
           nsim,
           par) {
    #delta: simutaneously shift on lambda's
    #o.w. shift one by one
    #par: rho, lam1, lam2, lam3
    par.old  <-  par
    par[1]  <-  par.old[1] * delta.rho
    if (delta != 1) {
      par[2:4]  <-  delta * par.old[2:4]
    }
    else{
      par[2]  <-  par.old[2] * delta1
      par[3]  <-  par.old[3] * delta2
      par[4]  <-  par.old[4] * delta3
    }
    X  <-  MO.sim(nsim, par)
    return(X)
  }








########################
########################


#shewhart



# MO model pmf, CDF
MO.logpmf <- function(a, b, c) {
  #x1, lam1, X1=a; x2, lam2, x2=b; x3, lam3, x3=b
  conv.log.gam(a + b + c + rho) - conv.log.gam(rho) - conv.log.gam(a + 1) -
    conv.log.gam(b + 1) - conv.log.gam(c + 1) + a * log(lam1) + b * log(lam2) +
    c * log(lam3) - (a + b + c + rho) * log(lam1 + lam2 + lam3 + 1)
}
MO.pmf <- function(a, b, c) {
  logpmf  <-  MO.logpmf(a, b, c)
  pmf  <-  exp(logpmf)
  return(pmf)
}
MO.cdf <-
  function(a, b, c) {
    #x1, lam1, X1<=a; x2, lam2, x2<=b; x3, lam3, x3<=b
    sum.mo <- 0
    if (a >= 0 & b >= 0 & c >= 0) {
      for (i in 0:a) {
        for (j in 0:b) {
          for (k in 0:c) {
            sum.mo <- sum.mo + MO.pmf(i, j, k)
            #print(sum.mo)
          }
        }
      }
    }
    return(sum.mo)
  }

# bnb pmf cdf
bnb.logpmf <- function(a, b, rho, lam1, lam2) {
  #x1=a,lam1;x2=b,lam2; position is important
  conv.log.gam(a + b + rho) - conv.log.gam(rho) - conv.log.gam(a + 1) -
    conv.log.gam(b + 1) + a * log(lam1) + b * log(lam2) - (a + b + rho) * log(lam1 +
                                                                                lam2 + 1)
}

bnb.pmf <- function(a, b, rho, lam1, lam2) {
  logpmf <- bnb.logpmf(a, b, rho, lam1, lam2)
  pmf <- exp(logpmf)
  return(pmf)
}

bnb.cdf <- function(a, b, rho, lam1, lam2) {
  #x1<=a,lam1;x2<=b,lam2; position is important
  sum = 0
  if (a >= 0 & b >= 0) {
    #iter=0
    for (i in 0:a) {
      for (j in 0:b) {
        sum <- sum + bnb.pmf(i, j, rho, lam1, lam2)
        #iter=iter+1
        #print(sum)
      }
    }
    #print(iter)
  }
  return(sum)
}

# shewheart control limits

# k is a global para you pick
# upper limit
u <- function(rho, lam) {
  # for each Xi i=1,2,3, has a lam, need to specify it
  # ooc is X>u, CDF using X<=u, so only have to return u
  u <- rho * lam + k * sqrt(lam * rho * (1 + lam))
  return(u)
}

# lower limit
l <- function(rho, lam) {
  # as we are using CDF to cal, NB are discrete
  # ooc is X<l, using cdf thus, take care of it here
  # if l is not integer, return l, CDF is using range, will consider it
  # else. return l-1
  l <- rho * lam - k * sqrt(lam * rho * (1 + lam))
  l <- ifelse(floor(l) != l, l, l - 1)
  return(l)
}

# prob(OOC) only upper limit
p <- function(par) {
  #x1, u1, par[2]=lam1, etc
  u1 <- u(par[1], par[2])
  u2 <- u(par[1], par[3])
  u3 <- u(par[1], par[4])
  Q <- MO.cdf(u1, u2, u3)
  p <- 1 - Q
  return(p)
}

# prob(OOC) only lower limit
lower.p <- function(par) {
  #marginal nbinomial k=rho,m=rholam
  rho <- par[1]
  lam1 <- par[2]
  lam2 <- par[3]
  lam3 <- par[4]
  l1 <- l(rho, lam1)
  l2 <- l(rho, lam2)
  l3 <- l(rho, lam3)
  
  p1 <- ifelse(l1 < 0, 0, pnbinom(l1, size = rho, mu = rho * lam1))
  p2 <- ifelse(l2 < 0, 0, pnbinom(l2, size = rho, mu = rho * lam2))
  p3 <- ifelse(l3 < 0, 0, pnbinom(l3, size = rho, mu = rho * lam3))
  term1 <- p1 + p2 + p3
  term2 <-
    bnb.cdf(l1, l2, rho, lam1, lam2) + bnb.cdf(l1, l3, rho, lam1, lam2) +
    bnb.cdf(l3, l2, rho, lam3, lam2)
  term3 <- MO.cdf(l1, l2, l3)
  p <- term1 - term2 + term3
  return(p)
}

# prob(OOC) w/ upper & lower
p.both <- function(par) {
  rho <- par[1]
  lam1 <- par[2]
  lam2 <- par[3]
  lam3 <- par[4]
  l1 <- l(rho, lam1)
  u1 <- u(rho, lam1)
  
  l2 <- l(rho, lam2)
  u2 <- u(rho, lam2)
  
  l3 <- l(rho, lam3)
  u3 <- u(rho, lam3)
  
  
  term1 <-
    pnbinom(u1, size = rho, mu = rho * lam1) + pnbinom(u2, size =
                                                         rho, mu = rho * lam2) + pnbinom(u3, size = rho, mu = rho * lam3)
  
  bnb.u1.u2 <- bnb.cdf(u1, u2, rho, lam1, lam2)
  bnb.u1.u3 <- bnb.cdf(u1, u3, rho, lam1, lam3)
  bnb.u2.u3 <- bnb.cdf(u3, u2, rho, lam3, lam2)
  term2 <- bnb.u1.u2 + bnb.u1.u3 + bnb.u2.u3
  #x1 limit must be in CDF 1st position etc. due to the def of CDF
  term31 <-
    MO.cdf(l1, l2, l3) - MO.cdf(l1, l2, u3) - MO.cdf(l1, u2, l3) -
    MO.cdf(u1, l2, l3)
  term32 <-
    MO.cdf(l1, u2, u3) + MO.cdf(u1, l2, u3) + MO.cdf(u1, u2, l3)
  term33 <- 1 - MO.cdf(u1, u2, u3)
  term3 <- term31 + term32 + term33
  
  s <- -term1 + term2 + term3
  return(s)
}



# ooc shewhart check
ooc.shewhart <- function(X, k, par, upper, lower) {
  #upper=T -> consider upper limit;
  #lower=T -> consider lower limit
  a <- X[1]
  b <- X[2]
  c <- X[3]
  rho <- par[1]
  lam1 <- par[2]
  lam2 <- par[3]
  lam3 <- par[4]
  l1 <- l(rho, lam1)
  u1 <- u(rho, lam1)
  
  l2 <- l(rho, lam2)
  u2 <- u(rho, lam2)
  
  l3 <- l(rho, lam3)
  u3 <- u(rho, lam3)
  
  if (!lower) {
    l1 <- -1
    l2 <- -1
    l3 <- -1
  }
  if (!upper) {
    u1 <- 10 ^ (10 ^ 10)
    u2 <- 10 ^ (10 ^ 10)
    u3 <- 10 ^ (10 ^ 10)
  }
  OOC <- ifelse(a <= l1 | a > u1 | b <= l2 | b > u2 | c <= l3 |
                  c > u3, T, F)
  return(OOC)
}












########################
########################

# SVs

# standardlize the data
std.data <-
  function(mydata,
           stable.mu = apply(mydata, 2, mean),
           stable.vcv = cov(mydata)) {
    #input data is matrix eg, 3-NB, ncol=3; this code works with 2-NB also
    #returns std.data, correlation matrix, stable.vcv,stable.mu
    
    mydata.std <- mydata
    for (i in 1:ncol(mydata)) {
      mydata.std[, i] <-
        (mydata[, i] - stable.mu[i]) / sqrt(stable.vcv[i, i])
    }
    cormat <- stable.vcv
    for (n1 in 1:nrow(stable.vcv)) {
      for (n2 in 1:ncol(stable.vcv)) {
        if (n1 == n2) {
          cormat[n1, n2] <- 1
        }
        else{
          cormat[n1, n2] <-
            stable.vcv[n1, n2] / sqrt(stable.vcv[n1, n1] * stable.vcv[n2, n2])
        }
      }
    }
    return(
      list(
        std.data <- mydata.std,
        stable.mu <- stable.mu,
        stable.vcv <- stable.vcv,
        cormat <- cormat
      )
    )
  }


# 3D plot of SV & kernel distance
plot3D <- function(mydata, sv.index) {
  library(scatterplot3d)
  #par(mfrow=c(1,1))
  test.plot <- scatterplot3d(
    mydata,
    pch = "",
    grid = T,
    box = F,
    angle = 45,
    main = "3D scatter plot"
  )
  source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
  addgrids3d(mydata, grid = c("xy", "xz", "yz"), angle = 45)
  test.plot$points3d(mydata[sv.index, ],
                     pch = 8,
                     type = "h",
                     col = "red")
  test.plot$points3d(mydata[-sv.index, ],
                     pch = 16,
                     type = "h",
                     col = "blue")
  legend(
    "bottom",
    legend = c("Support Vectors", "not Support Vectors"),
    pch = c(8, 16),
    col = c("red", "blue"),
    inset = -0.25,
    xpd = TRUE,
    horiz = TRUE,
    bty = "n"
  )
}


# pairwise plots of SV
plot.pair <- function(mydata, sv.index, wts, nsv) {
  par(mfrow = c(1, 3))
  index.plot <- matrix(data = c(1, 2, 2, 3, 3, 1),
                       ncol = 2,
                       byrow = T)
  for (i in 1:3) {
    #plots L->R: X1 X2, X2 X3, X3 X1
    x <- index.plot[i, 1]
    y <- index.plot[i, 2]
    plot(
      mydata[, x],
      mydata[, y],
      type = "n",
      xlim = c(min(mydata[, x]) - .5, max(mydata[, x] + .5)),
      ylim = c(min(mydata[, y]) - .5, max(mydata[, y] + .5)),
      xlab = c("X", x),
      ylab = c("X", y),
      main = ifelse(i == 2, "pairwise plots of SVs with weights", "")
    )
    text(mydata[, x], mydata[, y], as.character(1:nrow(mydata)), cex = .75)
    points(
      mydata[sv.index, x],
      mydata[sv.index, y],
      cex = wts * 2 * nsv,
      col = 2,
      lwd = 2
    )
  }
}


# get support vectors
getsv <-
  function(mydata,
           stable.mu = apply(mydata, 2, mean),
           stable.vcv = cov(mydata),
           nu = 0.1,
           scaled = F,
           tol = .001,
           kern.wt = 1) {
    library(kernlab)
    # nu upper bound for training error and lower bound on the fraction of data
    #    to become support vectors (0.1 default)
    # incooperate with 2-NB and 3-NB,
    # if need more than that, plots need to be justified whereas others not
    
    list <- std.data(mydata, stable.mu, stable.vcv)
    
    mydata.std <- list$std.data
    stable.mu <- list$stable.mu
    stable.vcv <- list$stable.vcv
    cormat <- list$cormat
    corinv <- solve(cormat)
    
    # Define the kernel distance.
    kernelf <- function(xvec, yvec) {
      dvec <- xvec - yvec
      dist <- (dvec) %*% corinv %*% dvec
      return(exp(-dist / kern.wt))
    }
    
    class(kernelf) = "kernel"
    
    # Get support vectors for the above kernel.
    # No need to scale as mydata.std is already scaled above.
    out1 <- ksvm(
      x = as.matrix(mydata.std),
      type = "one-svc",
      kernel = kernelf,
      nu = nu,
      scaled = scaled,
      tol = tol
    )
    nsv <- out1@nSV
    svmat <- out1@xmatrix
    wts <- out1@alpha
    sv.index <-
      out1@alphaindex # which of the data are the support vectors
    
    #adjust weights so they add to 1
    wts <- wts / sum(wts)
    
    #3d plot of my SV & pairwise plots
    if (ncol(mydata) == 3) {
      par(mfrow = c(1, 2))
      plot3D(mydata.std, sv.index)
      #plot.pair(mydata.std,sv.index,wts,nsv)
      #par(mfrow=c(1,1))
      plot(wts,
           cex = wts * 2 * nsv,
           main = "Support Vectors with weights",
           col = "red")
      text(1:nrow(mydata), wts, as.character(1:nrow(mydata)), cex = .75)
    }
    
    if (ncol(mydata) == 2) {
      par(mfrow = c(1, 2))
      plot(
        mydata.std,
        type = "n",
        xlim = c(min(mydata.std[, 1]) - .5, max(mydata.std[, 1] + .5)),
        ylim = c(min(mydata.std[, 2]) - .5, max(mydata.std[, 2] + .5))
      )
      text(mydata.std, as.character(1:nrow(mydata.std)), cex = .75)
      points(mydata.std[sv.index, ],
             cex = wts * 2 * nsv,
             col = 2,
             lwd = 2)
      plot(wts, cex = wts * 2 * nsv)
      par(mfrow = c(1, 1))
    }
    
    # Compute constant 3rd term in kernel distance radical.
    p2 <- 0
    for (i in 1:nsv) {
      for (j in 1:nsv) {
        p2 <- p2 + wts[i] * wts[j] * kernelf(svmat[i, ], svmat[j, ])
      }
    }
    
    return(
      list(
        mydata = mydata,
        stable.mu = stable.mu,
        stable.vcv = stable.vcv,
        stable.cormat = cormat,
        mydata.std = mydata.std,
        svmatrix = svmat,
        nsv = nsv,
        alpha = wts,
        sv.index = sv.index,
        kernelf = kernelf,
        nu = nu,
        term3 = p2,
        kern.wt = kern.wt
      )
    )
  }


# to check which nu gives the most close kernerl distance
check.kdist <- function(ksvout) {
  # Check if kernel distance is the same for each support vector.
  cormat <- ksvout$stable.cormat
  corinv <- solve(cormat)
  mydata.std <- as.matrix(ksvout$mydata.std)
  ndata <- nrow(mydata.std)
  sv.index <- ksvout$sv.index
  svmat <- mydata.std[sv.index, ]
  nsv <- nrow(svmat)
  wts <- ksvout$alpha
  term3 <- ksvout$term3
  mykern <- ksvout$kernelf
  kern.wt <- ksvout$kern.wt
  
  term1 <- 1 #k(z,z)=1
  kdist <- c()
  for (i in 1:ndata) {
    zmat <- matrix(rep(mydata.std[i, ], nsv),
                   ncol = ncol(mydata.std),
                   byrow = T)
    term2 <-
      -2 * sum(wts * exp(-diag((zmat - svmat) %*% corinv %*% t(zmat -
                                                                 svmat)
      ) / kern.wt))
    kdist[i] <- sqrt(term1 + term2 + term3)
  }
  kernerl.mean <- mean(kdist[sv.index])
  kernerl.stddev <- sd(kdist[sv.index])
  plot(kdist)
  lines(kdist)
  points(sv.index, kdist[sv.index], col = 2, cex = wts * 2 * nsv)
  abline(h = kernerl.mean)
  #print kernerl distance mean and std.dev on plots
  mtext(text = "kernerl mean",
        side = 3,
        line = -1)
  mtext(text = round(kernerl.mean, 4),
        side = 3,
        line = -2)
  mtext(text = "kernerl std.dev",
        side = 3,
        line = 2)
  mtext(text = round(kernerl.stddev, 4),
        side = 3,
        line = 1)
  return(kdist[sv.index])
}


#stable process
kdist.std <- function(zvec, ksvout) {
  library(kernlab)
  # zvec new standardized value
  # vcv of Xvec-Yvec is 2*vcv
  cormat <- ksvout$stable.cormat
  corinv <- solve(cormat)
  mydata.std <- ksvout$mydata.std
  sv.index <- ksvout$sv.index
  svmat <- as.matrix(mydata.std[sv.index, ])
  nsv <- nrow(svmat)
  wts <- ksvout$alpha
  term3 <- ksvout$term3
  mykern <- ksvout$kernelf
  
  term1 <- 1 #k(z,z)=1
  zmat <- matrix(rep(zvec, nsv),
                 ncol = ncol(mydata.std),
                 byrow = T)
  term2 <-
    -2 * sum(wts * exp(-diag((zmat - svmat) %*% corinv %*% t(zmat -
                                                               svmat)
    )))
  
  return(sqrt(term1 + term2 + term3))
}


# simulate kernerl distance
kdist.simu <- function(ksvout, N = 200, par) {
  #ksvout come from getsv
  #par are rho, lam1, lam2, lam3, come from nlimb outout
  kvec <- c()
  #simulate data using MO model & std it
  stable.mu <- ksvout$stable.mu
  stable.vcv <- ksvout$stable.vcv
  newdata <- MO.sim(N, par)
  newdata <- std.data(newdata)$std.data
  
  for (i in 1:N) {
    kvec[i] <- kdist.std(newdata[i, ], ksvout)
  }
  
  par(mfrow = c(1, 2))
  hist(log(kvec), xlab = "log(Kernel Distance)")
  plot(ecdf(kvec), xlab = "Kernel Distance")
  par(mfrow = c(1, 1))
  
  k.ecdf <- ecdf(kvec)
  kdist.vals <- sort(unique(kvec))
  kdist.cdf <- k.ecdf(kdist.vals)
  return(list(
    kdist.vals = kdist.vals,
    kdist.cdf = kdist.cdf,
    k.ecdf = k.ecdf
  ))
}


#check OOC for kernerls, compared with stable process
ooc.kernerl <- function(X, ksvout, k) {
  # X is from ooc simulation;
  # ksvout is from getsv
  # k is ucl you choose from ucl=1-1/arl
  # eg 3-NB ncol(X)=3
  ksvout <- kernerl.svout
  stable.mu <- ksvout$stable.mu
  stable.vcv <- ksvout$stable.vcv
  X <-
    std.data(X, stable.mu = stable.mu, stable.vcv = stable.vcv)$std.data
  kdist <-  kdist.std(X, ksvout)[1, 1]
  ooc <-  ifelse(kdist > k, T, F)
  return(list(ooc = ooc, kdist = kdist))
}









########################
########################

# LRs



#simulate ooc and stable data, get betas
get.betas <- function(n, par, delta.m, pos = 0) {
  # delta.m for built the model
  # pos=0, default, simultaneous shift on all lam
  # o.w, pos= 1,2,3,4. 4 for rho, 1-3 for lam1-lam3
  library(uwIntroStats)
  # n(ooc)=n(stable)=n for simplicity
  X.stable <-  MO.sim(n, par)
  if (pos == 0) {
    X.ooc <-  ooc.sim(delta = delta.m,
                      nsim = n,
                      par = par)
  }
  else if (pos == 1) {
    X.ooc <- ooc.sim(delta1 = delta.m,
                     nsim = n,
                     par = par)
  }
  else if (pos == 2) {
    X.ooc <- ooc.sim(delta2 = delta.m,
                    nsim = n,
                    par = par)
  }
  else if (pos == 3) {
    X.ooc <-  ooc.sim(delta3 = delta.m,
                    nsim = n,
                    par = par)
  }
  else if (pos == 4) {
    X.ooc <-  ooc.sim(delta.rho = delta.m,
                    nsim = n,
                    par = par)
  }
  X.lr <-  rbind(X.stable, X.ooc)
  Y.lr <-  c(rep(0, n), rep(1, n))
  data.lr <-  cbind(X.lr, Y.lr)
  data.lr <-  as.data.frame(data.lr)
  colnames(data.lr) <-  c("X1", "X2", "X3", "Y")
  fit <-  regress("odds", Y ~ X1 + X2 + X3, data = data.lr)
  betas <-  fit$coefficients[, 1]
  #fit=glm(Y~X1+X2+X3,data=data.lr,family=binomial())
  #betas=coef(fit)
  return(betas)
}

get.betas.ave <- function(n, par, delta.m, rep = 100, pos) {
  # ave logistic regression results of betas
  beta.mat <-  get.betas(n, par, delta.m, pos = pos)
  beta.mat <-  as.matrix(beta.mat)
  while (ncol(beta.mat) < rep) {
    beta <-  get.betas(n, par, delta.m, pos = pos)
    beta.mat <-  cbind(beta.mat, beta)
  }
  betas <-  rowMeans(beta.mat)
  return(betas)
}


#logistic regression, get prob(ooc)
p.ooc.log <-  function(betas, X) {
  # X is BNB or TNB input, thus to do matrix cal, prepare it
  betas <-  as.matrix(betas)
  X <-  cbind(rep(1, nrow(X)), X)
  log.odds <-  X %*% betas
  odds <-  exp(log.odds)
  ps <-  odds / (1 + odds)
  return(ps)
}


ooc.log <- function(betas, X, ucl) {
  # X 1 row 3 col TNB/ 2 col BNB
  p <-  as.numeric(p.ooc.log(betas, X))
  ooc <-  ifelse(p > ucl, T, F)
  return(list(ooc = ooc, pooc = p))
}




########################
########################

# LDAs

# get lda ratio
ooc.logpmf <- function(a, b, c, delta.m, pos = 0) {
  #x1, lam1, X1=a; x2, lam2, x2=b; x3, lam3, x3=b
  if (pos == 0) {
    lam1 <-  lam1 * delta.m
    lam2 <-  lam2 * delta.m
    lam3 <-  lam3 * delta.m
  }
  else if (pos == 1) {
    lam1 <-  lam1 * delta.m
  }
  else if (pos == 2) {
    lam2 <-  lam2 * delta.m
  }
  else if (pos == 3) {
    lam3 <-  lam3 * delta.m
  }
  else if (pos == 4) {
    rho <-  rho * delta.m
  }
  else {
    print("Postition is wrong, only 0-4.")
    break
  }
  conv.log.gam(a + b + c + rho) - conv.log.gam(rho) - conv.log.gam(a + 1) -
    conv.log.gam(b + 1) - conv.log.gam(c + 1) + a * log(lam1) + b * log(lam2) +
    c * log(lam3) - (a + b + c + rho) * log(lam1 + lam2 + lam3 + 1)
}

check.ratio <- function(x, delta.m, pos) {
  a <-  x[1]
  b <-  x[2]
  c <-  x[3]
  stable <-  MO.logpmf(a, b, c)
  ooc <-  ooc.logpmf(a, b, c, delta.m, pos = pos)
  diff <-  ooc - stable
  return(diff)
}


# simulate LDA log-likelihood ratio
# only consider n(ooc)=n(stable)
lda.ratio.sim <- function(X, delta.m, pos) {
  nsim <-  nrow(X)
  lda.ratio <-  c()
  for (i in 1:nsim) {
    lda.ratio[i] <-  check.ratio(X[i, ], delta.m, pos = pos)
  }
  return(list(
    nsim = nsim,
    lda.ratio = lda.ratio,
    delta.m = delta.m,
    pos = pos
  ))
}


# check ooc for lda
ooc.lda <- function(x, ucl, delta.m, pos) {
  ratio <-  check.ratio(x, delta.m, pos = pos)
  ooc <-  ifelse(ratio > ucl, T, F)
  return(ooc)
}


# get arl with certain delta.m: model shift need get ucl from delta.m
ooc.arl.lda.m <-
  function(deltas, ntimes, par, ucl, delta.m, pos = 0) {
    ooc.arl.lda <-  matrix(NA, nrow = length(deltas), ncol = ntimes)
    for (i in 1:length(deltas)) {
      delta <-  deltas[i]
      for (j in 1:ntimes) {
        ooc <-  F
        r <-  0
        
        while (!ooc) {
          X <-  ooc.sim(delta = delta,
                      nsim = 1,
                      par = par)
          ooc <-  ooc.lda(X, ucl, delta.m, pos)
          r <-  r + 1
          #print(r)
        }
        ooc.arl.lda[i, j] <-  r
      }
    }
    ooc.arl.lda <-  rowMeans(ooc.arl.lda)
    return(ooc.arl.lda)
  }
