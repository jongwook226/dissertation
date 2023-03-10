#setwd("/Users/jongwookkim/Library/CloudStorage/OneDrive-IndianaUniversity/IU/Dissertation")
setwd("/Users/jongwookkim/Documents/dissertation/code")
source("230131_functions_dissertation.R")
source("irfkrg-fcn01-sphere-harmonics.R")
source("irfkrg-fcn02-legendre-polynomials.R")


sim1 <- function(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1){
  
  # a1=c(0.9,0.9,0.9); a2=c(0.001,0.001,0.001); a3=c(1,1,1)
  # a1=c(0.8,0.6,0.5); a2=c(0.1,0.2,0.3); a3=c(10,5,3)
  # a1=c(0.5,0.5,0.5); a2=c(0.2,0.2,0.2); a3=c(1,1,1)
  # 
  # P <- 300 #location data size
  # t <- 50 #time data size
  # 
  # P <- 200 #location data size
  # t <- 20 #time data size
  # kappa <- 1 #IRF order
  
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
  
  #create tau data.frame.
  
  
  if(kappa==0){
    tau <- data.frame(tau_lat=0, tau_long=0)
  }else{
    tau <- data.frame(tau_lat=rep(0,kappa^2), tau_long=rep(0,kappa^2))
    
    # #using the same tau
    #   tau$tau_lat <- runif(1, min=-pi/2, max=pi/2) #tau_latitude
    #   tau$tau_long <- runif(1, min=0, max=2*pi) #tau_longitude
    
    #using different taus
    for(i in 1:kappa^2){
      #random sampling for tau
      tau$tau_lat[i] <- runif(1, min=-pi/2, max=pi/2) #tau_latitude
      tau$tau_long[i] <- runif(1, min=0, max=2*pi) #tau_longitude
    }
  }
  
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
  #tau_spherical <- rep(0,4)
  
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
    #tau_spherical <- c(Y00(c(tau_lat,tau_long)), Y1n1(c(tau_lat,tau_long)), Y10(c(tau_lat,tau_long)), Y11(c(tau_lat,tau_long)))
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
  
  #Spherical distance b/w dataset and taus
  if(kappa==0){
    tauDist <- bwtauDist <- 0
  }else{
    tauDist <- data.frame(matrix(rep(0, P*(kappa^2)), ncol=kappa^2))
    
    if(kappa != 0){ #we do not need tau if kappa =0
      for(i in 1:P){
        for(j in 1:(kappa^2)){
          tauDist[i,j] <- SphDist(dat[i,1], dat[i,2], tau$tau_lat[j], tau$tau_long[j]) 
        }
      }
    }
    
    tauDist <- as.data.frame(lapply(tauDist, rep, t))
    for(i in 1:kappa^2){
      colnames(tauDist)[i] <- paste("tauDist",i, sep="")
    }
    
    #Spherical distance b/w taus
    bwtauDist <- matrix(0, nrow=kappa^2, ncol=kappa^2)
    for(i in 1:kappa^2){
      for(j in i:kappa^2){
        bwtauDist[i,j] <- bwtauDist[j,i] <- SphDist(tau$tau_lat[i], tau$tau_long[i], tau$tau_lat[j], tau$tau_long[j]) 
      }
    }
  }
  
  index <- 1:P #index of spatial locations (we use it to avoid redundancy of computation of distances)
  dat <- data.frame(dat,time,index)
  
  
  ####################### Get covariance function by cpp code ######################
  #dyn.load("C_dissertation.so")
  #dyn.unload("C_dissertation.so")
  
  R <- C_Cov_mat(dat=dat, sDmat, a1, a2, a3, kappa, low_spherical, tauDist, bwtauDist)
  
  # #check R is positive semi definite
  #  e <- eigen(R)$values
  #  which(e<=0)
  #  e[which(e<=0)]
  #  summary(e)
  #  length(e)
  # 
  #  matrixcalc::is.positive.semi.definite(R, tol=1e-8)
  # 
  #  isSymmetric(R) #Check symmetry
  ################################################################################
  
  # R version function
  # R2 <- Cov_mat(dat, a1, a2, a3, kappa, tau$tau_lat, tau$tau_long, low_spherical, tau_spherical) #covariance matrix with R version
  # all.equal(round(R,4),round(R2,4),tolerance=1e-8) #check R version and cpp version have the same values
  # R[1:7,1:7]; R2[1:7,1:7]
  # 
  # e2 <- eigen(R2)$values
  # which(e2<=0)
  # e2[which(e2<=0)]
  # summary(e2)
  
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
  
  # H <- 20         # Number of distances at which to estimate
  # Hvec <- as.vector(seq(0, 1.5, by = 1/H))[1:H]
  # eps <- 0.05     # use to bin distances in estimation
  
  H <- 20         # Number of distances at which to estimate
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  eps <- 0.1     # use to bin distances in estimation

  
  
  ####################### Get MOM estimate by cpp code ######################
  #dyn.load("C_dissertation.so")
  #dyn.unload("C_dissertation.so")
  
  #t_lag <- ceiling(t/2) #time lag
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
        R_mat[i,j] <- Rvec(dist=Hvec[i], time=j-1, a1=a[1], a2=a[2], a3=a[3])
      }
    }
    sum((R_mom- R_mat)^2) #without the weights
    # sum(W*(R_mom- R_mat)^2)
  }
  
  # # Without the scale parameter
  # obj2 <- function(a){
  #   Rvec <- function(dist, time, a1=a[1], a2=a[2]){
  #     Ct <- a[1]*exp(-a[2]*abs(time))
  #     R <- (1 - Ct^2)/(1-2*cos(dist)*Ct+Ct^2)^(3/2)
  #     return(R)
  #   }
  #   
  #   #theoretical covariance matrix
  #   R_mat <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
  #   for(i in 1:length(Hvec)){
  #     for(j in 0:ceiling(t/2-1)){
  #       R_mat[i,j+1] <- Rvec(dist=Hvec[i], time=j, a1=a[1], a2=a[2])
  #     }
  #   }
  #   if(kappa==1){
  #     R_mat <- R_mat - 1/(4*pi)
  #   }
  #   #sum(W*(R_mom- R_mat)^2)
  #   sum((R_mom- R_mat)^2)
  # }
  
  fitted_para <- suppressWarnings(nlminb(c(0,0,0),obj,lower=c(0,0,0), upper=c(1,Inf,Inf)))$par #fitted parameters by nlminb
  
  
  #fitted covariance matrix
  R_fitted <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
  for(i in 1:length(Hvec)){
    for(j in 1:t_lag){
      R_fitted[i,j] <- Rvec(dist=Hvec[i], time=j-1, fitted_para[1], fitted_para[2], fitted_para[3], kappa)
    }
  }
  
  result <- list(R_mat=R_mat, R_fitted=R_fitted, R_mom=R_mom, para_est=fitted_para) #returns true, fitted, Mom covariance matrices, and parameter estimates.
  return(result)
}





#########################  Run the code  ###########################

######### To see the effect of the parameter p_1 ############
result1 <- sim1(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
result2 <- sim1(a1=c(0.5,0.5,0.5), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
par(mfrow=c(1,2))
s_plot(1, result1[[1]], result1[[2]], result1[[3]])
s_plot(1, result2[[1]], result2[[2]], result2[[3]])
par(mfrow=c(1,1))
#############################################################

######### To see the effect of the parameter p_2 ############
result1 <- sim1(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
result2 <- sim1(a1=c(0.8,0.8,0.8), a2=c(0.4,0.4,0.4), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
par(mfrow=c(1,2))
t_plot(1, result1[[1]], result1[[2]], result1[[3]])
t_plot(1, result2[[1]], result2[[2]], result2[[3]])
par(mfrow=c(1,1))
#############################################################


result1 <- sim1(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
result2 <- sim1(a1=c(0.1,0.1,0.1), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #too small parameters. the graphs are flat (not decaying)
result3 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.3,0.3,0.3), a3=c(1,1,1), P=200, t=20, kappa=1) #decay too fast??
result4 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.9,0.9,0.9), a3=c(1,1,1), P=200, t=20, kappa=1) #decay too fast??
result5 <- sim1(a1=c(0.5,0.5,0.5), a2=c(0.2,0.2,0.2), a3=c(1,1,1), P=200, t=20, kappa=1) #a1 too small not converging.
result6 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.001,0.001,0.001), a3=c(1,1,1), P=200, t=20, kappa=1) #check t_plot
result7 <- sim1(a1=c(0.1,0.1,0.1), a2=c(0.8,0.8,0.8), a3=c(1,1,1), P=200, t=20, kappa=1) #check t_plot

#p_1 seems related to overall decay rate.
#p_2 seems more related to temporal terms.

result <- result1

result$para_est #parameter estimates


###################################### Contour Plots ####################################
library(ggplot2)

result <- result1

p1 <- contour_plot(result[[1]], 1) #True Values
p2 <- contour_plot(result[[2]], 2) #Fitted Values
p3 <- contour_plot(result[[3]], 3) #MoM

ggpubr::ggarrange(p1,p2,p3, ncol=2, nrow=2)

# library(cowplot)
# plot_grid(p1,p2,p3, labels=c("A", "B","C"), ncol = 2, nrow = 2)



########################################################

result <- result1

par(mfrow=c(2,3))
  
s_plot(1, result[[1]], result[[2]], result[[3]])
s_plot(2, result[[1]], result[[2]], result[[3]])
s_plot(3, result[[1]], result[[2]], result[[3]])
#s_plot(4, result[[1]], result[[2]], result[[3]])
#s_plot(5, result[[1]], result[[2]], result[[3]])
s_plot(6, result[[1]], result[[2]], result[[3]])
#s_plot(7, result[[1]], result[[2]], result[[3]])
s_plot(8, result[[1]], result[[2]], result[[3]])
#s_plot(9, result[[1]], result[[2]], result[[3]])
s_plot(10, result[[1]], result[[2]], result[[3]])


t_plot(1, result[[1]], result[[2]], result[[3]])
t_plot(2, result[[1]], result[[2]], result[[3]])
t_plot(3, result[[1]], result[[2]], result[[3]])
t_plot(4, result[[1]], result[[2]], result[[3]]) 
t_plot(5, result[[1]], result[[2]], result[[3]])
#t_plot(6, result[[1]], result[[2]], result[[3]])
t_plot(7, result[[1]], result[[2]], result[[3]])
#t_plot(8, result[[1]], result[[2]], result[[3]])
#t_plot(9, result[[1]], result[[2]], result[[3]])
#t_plot(10, result[[1]], result[[2]], result[[3]])

par(mfrow=c(1,1))



########################################################
par(mfrow=c(2,3))
# Our fitted value is better than mom for big lags!!!????

s_plot2(1, result[[1]], result[[2]], result[[3]])
s_plot2(2, result[[1]], result[[2]], result[[3]])
s_plot2(3, result[[1]], result[[2]], result[[3]])
#s_plot2(4, result[[1]], result[[2]], result[[3]])
#s_plot2(5, result[[1]], result[[2]], result[[3]])
s_plot2(6, result[[1]], result[[2]], result[[3]])
#s_plot2(7, result[[1]], result[[2]], result[[3]])
s_plot2(8, result[[1]], result[[2]], result[[3]])
#s_plot2(9, result[[1]], result[[2]], result[[3]])
s_plot2(10, result[[1]], result[[2]], result[[3]])


t_plot2(1, result[[1]], result[[2]], result[[3]])
t_plot2(2, result[[1]], result[[2]], result[[3]])
t_plot2(3, result[[1]], result[[2]], result[[3]])
t_plot2(4, result[[1]], result[[2]], result[[3]]) 
t_plot2(5, result[[1]], result[[2]], result[[3]])
#t_plot2(6, result[[1]], result[[2]], result[[3]])
t_plot2(7, result[[1]], result[[2]], result[[3]])
#t_plot2(8, result[[1]], result[[2]], result[[3]])
#t_plot2(9, result[[1]], result[[2]], result[[3]])
#t_plot2(10, result[[1]], result[[2]], result[[3]])

par(mfrow=c(1,1))
########################################################
par(mfrow=c(2,3))
s_plot(1, result1[[1]], result1[[2]], result1[[3]])
s_plot(1, result2[[1]], result2[[2]], result2[[3]])
s_plot(1, result3[[1]], result3[[2]], result3[[3]])
s_plot(1, result4[[1]], result4[[2]], result4[[3]])
s_plot(1, result5[[1]], result5[[2]], result5[[3]])
s_plot(1, result6[[1]], result6[[2]], result6[[3]])
s_plot(1, result7[[1]], result7[[2]], result7[[3]])

t_plot(1, result1[[1]], result1[[2]], result1[[3]])
t_plot(1, result2[[1]], result2[[2]], result2[[3]])
t_plot(1, result3[[1]], result3[[2]], result3[[3]])
t_plot(1, result4[[1]], result4[[2]], result4[[3]])
t_plot(1, result5[[1]], result5[[2]], result5[[3]])
t_plot(1, result6[[1]], result6[[2]], result6[[3]])
t_plot(1, result7[[1]], result7[[2]], result7[[3]])

# result1 <- sim1(a1=c(0.8,0.8,0.8), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #ideal situation
# result2 <- sim1(a1=c(0.1,0.1,0.1), a2=c(0.1,0.1,0.1), a3=c(1,1,1), P=200, t=20, kappa=1) #too small parameters. the graphs are flat (not decaying)
# result3 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.3,0.3,0.3), a3=c(1,1,1), P=200, t=20, kappa=1) #decay too fast??
# result4 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.9,0.9,0.9), a3=c(1,1,1), P=200, t=20, kappa=1) #decay too fast??
# result5 <- sim1(a1=c(0.5,0.5,0.5), a2=c(0.2,0.2,0.2), a3=c(1,1,1), P=200, t=20, kappa=1) #a1 too small not converging.
# result6 <- sim1(a1=c(0.9,0.9,0.9), a2=c(0.001,0.001,0.001), a3=c(1,1,1), P=200, t=20, kappa=1) #check t_plot
# result7 <- sim1(a1=c(0.1,0.1,0.1), a2=c(0.8,0.8,0.8), a3=c(1,1,1), P=200, t=20, kappa=1) #check t_plot

par(mfrow=c(1,1))
