setwd("/Users/jongwookkim/Library/CloudStorage/OneDrive-IndianaUniversity/IU/Dissertation")
source("102822_functions_dissertation.R")
source("irfkrg-fcn01-sphere-harmonics.R")
source("irfkrg-fcn02-legendre-polynomials.R")

sim2 <- function(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1){
  
   # a1=c(0.8,0.8,0.8); a2=c(0.1,0.1,0.1); a3=c(1,1,1)
   # a1=c(0.8,0.6,0.5); a2=c(0.1,0.2,0.3); a3=c(10,5,3)
   # 
   # P <- 300 #location data size
   # t <- 50 #time data size
   # 
   # P <- 100 #location data size
   # t <- 10 #time data size
   # kappa <- 2 #IRF order

   
  # a1 <- c(0.8,0.7,0.6) #should be a1<1 #a parameter for spatial term
  # a2 <- c(0.2,0.3,0.5) #a parameter for time
  # a3 <- c(10,5,3) #scale parameter

  #Check constrains of the parameter values
  if(a1[1] < a1[3] || a1[1] < a1[3]){
    stop("\nwrong parameter values!\n")
  }
  if(a2[1] > a2[3] || a2[1] > a2[3]){
    stop("\nwrong parameter values!\n")
  }
  if(a3[1] < a3[3] || a3[1] < a3[3]){
    stop("\nwrong parameter values!\n")
  }
  
  
  #random sampling for tau
  tau_lat <- runif(1, min=-pi/2, max=pi/2) #tau_latitude
  tau_long <- runif(1, min=0, max=2*pi) #tau_longitude
  
  n <- P*t #total data size
  
  lat <- runif(P, min=-pi/2, max=pi/2) #latitude
  long <- runif(P, min=0, max=2*pi) #longitude
  time <- rep(1:t, each=P)
  
  ####################### Get Spherical distance matrix by cpp code ######################
  dyn.load("C_dissertation.so")
  #is.loaded('SphDist')
  #is.loaded('Cov_mat')
  #is.loaded('G_hat')
  #dyn.unload("C_dissertation.so")
  dat <- data.frame(lat,long)
  
  #save spherical distance matrix
  sDmat <- C_SphDist(dat)
  
  #save the values of Spherical harmonics of the truncated parts
  low_spherical <- data.frame(matrix(0, ncol=4, nrow=nrow(dat)*t))
  colnames(low_spherical) <- c("Y00", "Y1n1", "Y10", "Y11")
  tau_spherical <- rep(0,4)
  
  if(kappa==2){
    
    low_spherical$Y00 <- rep(apply(dat[,1:2],1,Y00),t)
    low_spherical$Y1n1 <- rep(apply(dat[,1:2],1,Y1n1),t)
    low_spherical$Y10 <- rep(apply(dat[,1:2],1,Y10),t)
    low_spherical$Y11 <- rep(apply(dat[,1:2],1,Y11),t)
    
    # low_spherical$Y00 <- rep(apply(dat[,1:2],1,Y.0.0),t)
    # low_spherical$Y1n1 <- rep(apply(dat[,1:2],1,Y.n1.1),t)
    # low_spherical$Y10 <- rep(apply(dat[,1:2],1,Y.0.1),t)
    # low_spherical$Y11 <- rep(apply(dat[,1:2],1,Y.1.1),t)
    
    #values of low order spherical harmonics for tau
    tau_spherical <- c(Y00(c(tau_lat,tau_long)), Y1n1(c(tau_lat,tau_long)), Y10(c(tau_lat,tau_long)), Y11(c(tau_lat,tau_long)))
    #tau_spherical <- c(Y.0.0(c(tau_lat,tau_long)), Y.n1.1(c(tau_lat,tau_long)), Y.0.1(c(tau_lat,tau_long)), Y.1.1(c(tau_lat,tau_long)))
  }
  
  # ####################### Get Spherical distance matrix R version ########################
  # sDmat2 <- matrix(nrow=nrow(dat), ncol=nrow(dat))
  # for(i in 1:nrow(dat)){
  #   for(j in i:nrow(dat)){
  #    sDmat2[i,j] <- sDmat[j,i] <- SphDist(dat[i,1],dat[i,2], dat[j,1],dat[j,2])
  #   }
  # }
  # ########################################################################################

  
  #Spherical distance with tau
  tauDist <- rep(0, P)
  if(kappa != 0){ #we do not need tau if kappa =0
    for(i in 1:nrow(dat)){
      tauDist[i] <- SphDist(dat[i,1], dat[i,2], tau_lat, tau_long) 
    }
  }

  
  index <- 1:P #index of spatial locations (we use it to avoid redundancy of computation of distances)
  dat <- data.frame(dat,time,tauDist,index)
  
  
  ####################### Get covariance function by cpp code ######################
  #dyn.load("C_dissertation.so")
  #dyn.unload("C_dissertation.so")
  
  R <- C_Cov_mat(dat=dat, sDmat, a1, a2, a3, kappa, low_spherical, tau_spherical)

  # #check R is positive semi definite
  #  e <- eigen(R)$values
  #  which(e<=0)
  #  e[which(e<=0)]
  #  summary(e)
  # 
  #  matrixcalc::is.positive.semi.definite(R, tol=1e-8)
  # 
  #  isSymmetric(R) #Check symmetry
  #################################################################################
# 
#   R2 <- Cov_mat(dat, a1, a2, a3, kappa, tau_lat, tau_long, low_spherical, tau_spherical) #covariance matrix with R version
#   all.equal(round(R,4),round(R2,4),tolerance=1e-8) #check R version and cpp version have the same values
#   R[1:5,1:5]; R2[1:5,1:5]
# 
#   e2 <- eigen(R2)$values
#   which(e2<=0)
#   e2[which(e2<=0)]
#   summary(e2)

  #random sampling by Cholesky
  #Cholesky decompostion only works for positive definite functions (not positive-semi definite)
  
  # Rsqrt <- chol(R) #R^(1/2)
  # Rsqrt <- chol(R, pivot=TRUE) #R^(1/2)
  # Z <- t(Rsqrt) %*% matrix(rnorm(n,0,1),n,1)
  
  Z <- mgcv::rmvn(n=1,mu=rep(0,n), V=R)
  
  dat$Z <- Z
  
# 
# #Check which function is faster
#   benchmark("rmvn"= {
#     mgcv::rmvn(n=1,mu=rep(0,n), V=R)
#   },"mvrnorm"={
#     Z <- MASS::mvrnorm(n=1,mu=rep(0,n), Sigma=R) #random sample from Gaussian process
#   }, "cholesky"={
#     Z = t(chol(R)) %*% matrix(rnorm(n,0,1),n,1)
#   },"svd"={
#   # Use svd for R^(1/2)
#   sv <- svd(R)
#   ud <- sv$u %*% diag(sqrt(sv$d))
# 
#   # Sample Z(x)
#   sigma <- 1 #scale parameter
#   Z <- sigma*ud %*% matrix(rnorm(n,0,1),n,1)
#   })
  
  
  ######################Compute MOM estimator#########################
  
  # H <- 10         # Number of distances at which to estimate
  # Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
  # eps <- 0.10     # use to bin distances in estimation
  
  # H <- 10         # Number of distances at which to estimate
  # Hvec <- as.vector(seq(0, 2, by=2/H))[1:H]
  # eps <- Hvec[2]-Hvec[1]     # use to bin distances in estimation
  
  #H <- 20         # Number of distances at which to estimate
  #Hvec <- as.vector(seq(0, 1.5, by = 1/H))[1:H]
  #eps <- 0.05     # use to bin distances in estimation
  
  H <- 20         # Number of distances at which to estimate
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  eps <- 0.1     # use to bin distances in estimation
  
  
  ####################### Get MOM estimate by cpp code ######################
  #dyn.load("C_dissertation.so")
  #dyn.unload("C_dissertation.so")
  
  # t_lag <- ceiling(t/2-1) #time lag
  t_lag <- 10
  
  mom_est <- C_G_hat(dat=dat, sDmat, Hvec=Hvec, eps=eps, P=P, t=t_lag, kappa=kappa, low_spherical)
  R_mom <- matrix(mom_est[[1]]/mom_est[[2]], nrow=length(Hvec)) #MOM estimate
  W <- matrix(mom_est[[2]], nrow=length(Hvec)) #weight
  
  # ############# Get MOM estimate R version##############
  # sDvec <- Dfun(Lmat=dat[,1:2], P, t) #spherical distances
  # 
  # tDmat <- as.matrix(dist(dat$time, method="manhattan")) #time distances
  # tDvec  <-tDmat[upper.tri(tDmat,diag = TRUE)]
  # 
  # mom_est2 <- G_Hat(dat, sDvec, Hvec, tDvec ,eps, P, t_lag, kappa, low_spherical) #MOM estimate
  # R_mom2 <- mom_est2$G
  # W2 <- 1/mom_est2$N #weight
  # ######################################################
  
  #theoretical covariance matrix
  R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
  
  for(i in 1:length(Hvec)){
    for(j in 1:t_lag){
      R_mat[i,j] <- Rvec(dist=Hvec[i], time=j-1, a1[1], a2[1], a3[1], kappa)
    }
  }
  
  
  #Loss function (object function) for optimization
  obj <- function(a){
    Rvec <- function(dist, time, a1=a[1], a2=a[2], a3=a[3]){
      Ct <- a[1]*exp(-a[2]*abs(time))
      R <- a[3]*(1 - Ct^2)/(1-2*cos(dist)*Ct+Ct^2)^(3/2)
      
      if(kappa >= 1) R <- R - a[3] * 1/(4*pi)
      if(kappa >= 2) R <- R - a[3] * Ct * 3/(4*pi) * P.1(cos(dist))
      #if(kappa >= 3) R <- R - a3 * 5*(Ct^2) / (4*pi) * P.2(cos(dist))
      return(R)
    }
    
    #theoretical covariance matrix
    R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
    for(i in 1:length(Hvec)){
      for(j in 1:t_lag){
        R_mat[i,j] <- Rvec(dist=Hvec[i], time=j-1 ,a1=a[1], a2=a[2], a3=a[3])
      }
    }
    sum((R_mom- R_mat)^2) #without the weights
    #sum(W*(R_mom- R_mat)^2)
  }
  
  
  fitted_para <- nlminb(c(0,0,0),obj)$par #fitted parameters by nlminb
  fitted_para2 <- suppressWarnings(nlminb(c(0,0,0),obj,lower=c(0,0,0), upper=c(1,Inf,Inf)))$par #fitted parameters by nlminb
  
  fitted_para3 <- optim(c(0,0,0), obj, method = "BFGS")$par #fitted parameters
  #fitted_para <- optim(c(0,0,0), obj, method = "L-BFGS-B", lower=c(0,0,-Inf), upper=c(1,Inf,Inf))$par #fitted parameters
  return(list(nlminb=fitted_para, nlminb2=fitted_para2, optim=fitted_para3))
}

sim_size <- 1000
para <- data.frame(matrix(ncol=3, nrow=sim_size))
colnames(para) <- c("a1","a2","a3")
para_list <- list(nlminb=para, nlminb2=para, optim=para)

set.seed(111)
for(i in 1:sim_size){
#for(i in i2){
  #para[i,] <- sim(a1=0.8, a2=0.1, a3=1, P=10, t=20)
  #cat("\n Iteration : ", i, "\n")
  #print(para[[1]][i,])

  a1=rep(0.10,3); a2=rep(0.8,3); a3=rep(0.5,3)
  
  para1 <- sim2(a1, a2, a3, P=200, t=20, kappa=2)
  para_list[[1]][i,] <- para1[[1]]
  para_list[[2]][i,] <- para1[[2]]
  para_list[[3]][i,] <- para1[[3]]
  cat("\n Iteration : ", i, "\n")
  print(para1)
}

para_list$nlminb2

#when nlminb fails
i2 <- which(para_list$nlminb2$a2 >= a2*10 | para_list$nlminb2$a3 >= a3*10 | para_list$nlminb2$a3 == 0)
#para_list$nlminb2 <- para_list$nlminb2[-i2,]
length(i2)

c(round(mean(para_list$nlminb2[,1][!is.na(para_list$nlminb2[,1])]),4), #p1
round(mean(para_list$nlminb2[,2][!is.na(para_list$nlminb2[,1])]),4), #p2
round(mean(para_list$nlminb2[,3][!is.na(para_list$nlminb2[,1])]),4)) #p3

# median(sort(para_list$nlminb2[,2])[1:1000]) #p2
# median(sort(para_list$nlminb2[,2])[1:500]) #p2
#plot(density(para_list$nlminb2[,2]),xlim=c(0,3)) #p2
# 
 median(para_list$nlminb2[,1][!is.na(para_list$nlminb2[,1])]) #p1
 median(para_list$nlminb2[,2][!is.na(para_list$nlminb2[,1])]) #p2
 median(para_list$nlminb2[,3][!is.na(para_list$nlminb2[,1])]) #p3
#
#
#
c(round(sd(para_list$nlminb2[,1][!is.na(para_list$nlminb2[,1])]),4), #p1
round(sd(para_list$nlminb2[,2][!is.na(para_list$nlminb2[,1])]),4), #p2
round(sd(para_list$nlminb2[,3][!is.na(para_list$nlminb2[,1])]),4)) #p3

plot(density(para_list$nlminb2[,1][!is.na(para_list$nlminb2[,1])])) #p1
abline(v=a1[1], col="red")
plot(density(para_list$nlminb2[,2][!is.na(para_list$nlminb2[,1])])) #p2
abline(v=a2[1], col="red")
plot(density(para_list$nlminb2[,3][!is.na(para_list$nlminb2[,1])])) #p3
abline(v=a3[1], col="red")

#save(para_list,file="sim2_01_08_05_kappa2.Rda")
#para_list2<- para_list


# mean(para_list2$nlminb2[1:14,1]) #p1
# mean(para_list2$nlminb2[1:14,2]) #p2
# mean(para_list2$nlminb2[1:14,3]) #p3
# 
# sd(para_list2$nlminb2[1:14,1]) #p1
# sd(para_list2$nlminb2[1:14,2]) #p2
# sd(para_list2$nlminb2[1:14,3]) #p3
# 
# plot(density(para_list2$nlminb2[1:14,1])) #p1
# plot(density(para_list2$nlminb2[1:14,2])) #p2
# plot(density(para_list2$nlminb2[1:14,3])) #p3
