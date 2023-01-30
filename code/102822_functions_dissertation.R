#####################################################
##############Homogeneous and Stationary#############
#####################################################

#Spherical distance
SphDist <- function(lat1, long1, lat2, long2){
  # Great Circle Formula
  cosd <- sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(long1 - long2)
  suppressWarnings(d <- acos(cosd)) 
  if(is.nan(d)) d <- 0 #acos has an issue that it returns NaN when d=1 although it actually should return 0.
  return(d)
}

#Covariance Function for stationary homogenous
Cov_fn00 <- function(dat1, dat2, a1, a2, a3=1, tau=c(pi/2,pi/2)){
  lat1 <- dat1[[1]]
  lat2 <- dat2[[1]]
  
  long1 <- dat1[[2]]
  long2 <- dat2[[2]]
  
  time1 <- dat1[[3]]
  time2 <- dat2[[3]]
  
  d <- SphDist(lat1, long1, lat2, long2) #spherical distance
  #Ct <- exp(-a2*abs(time1-time2))
  #R <- a1/(4*pi)*(1-Ct^2)/(1-2*d*Ct+Ct^2)^(3/2)
  
  Ct <- a1[1]*exp(-a2[1]*abs(time1-time2))
  R <- a3[1]*(1 - Ct^2)/(1-2*cos(d)*Ct+Ct^2)^(3/2)
  #R <- (1 - Ct^2)/(1-2*cos(d)*Ct+Ct^2)^(3/2) #without scale parameter
  return(R)
}

#Covariance Function for IRF(1)/I(0)
Cov_fn10 <- function(dat1, dat2, a1, a2, a3, tau){
  lat1 <- dat1[[1]]
  lat2 <- dat2[[1]]
  
  long1 <- dat1[[2]]
  long2 <- dat2[[2]]
  
  time1 <- dat1[[3]]
  time2 <- dat2[[3]]
  
  #d <- SphDist(lat1, long1, lat2, long2) #spherical distance
  #Ct <- exp(-a2*abs(time1-time2))
  #R <- a1/(4*pi)*(1-Ct^2)/(1-2*d*Ct+Ct^2)^(3/2)
  
  
  ICF <- function(d,i){ #i is index for the parameters
    Ct <- a1[i]*exp(-a2[i]*abs(time1-time2))
    icf <- (1 - Ct^2)/((1-2*cos(d)*Ct+Ct^2)^(3/2))
    return(icf)
  } 
  
  ##new version
  #R <- a3[1]*ICF(SphDist(lat1,long1, lat2, long2),1) + a3[2]*1/(4*pi)*ICF(0,2) + a3[3]*1/(2*sqrt(pi))*ICF(SphDist(lat1,long1,tau[1],tau[2]),3) + a3[3]*1/(2*sqrt(pi))*ICF(SphDist(lat2,long2,tau[1],tau[2]),3) - 1/(4*pi) - 1/(16*pi^2) - 1/(4*pi^(3/2))
  R <- a3[1]*ICF(SphDist(lat1,long1, lat2, long2),1) + a3[2]/(4*pi)*ICF(0,2) + a3[3]/(2*sqrt(pi))*ICF(SphDist(lat1,long1,tau[1],tau[2]),3) + a3[3]/(2*sqrt(pi))*ICF(SphDist(lat2,long2,tau[1],tau[2]),3) - a3[1]/(4*pi) - a3[2]/(16*pi^2) - a3[3]/(4*pi^(3/2))

  return(R)
}

#Covariance Function for IRF(2)/I(0)
Cov_fn20 <- function(dat1, dat2, a1, a2, a3, tau, low_spherical1, low_spherical2, tau_spherical){
  #i=1;j=50;dat1 <- dat[i,]; dat2 <- dat[j,]; tau=c(tau_lat, tau_long); low_spherical1=low_spherical[i,]; low_spherical2=low_spherical[j,]
  
  lat1 <- dat1[[1]]
  lat2 <- dat2[[1]]
  
  long1 <- dat1[[2]]
  long2 <- dat2[[2]]
  
  time1 <- dat1[[3]]
  time2 <- dat2[[3]]
  
  #d <- SphDist(lat1, long1, lat2, long2) #spherical distance
  #Ct <- exp(-a2*abs(time1-time2))
  #R <- a1/(4*pi)*(1-Ct^2)/(1-2*d*Ct+Ct^2)^(3/2)
  
  Ct <- function(i){
    Ct <- a1[i]*exp(-a2[i]*abs(time1-time2))
    return(Ct)
  }
  
  ICF <- function(d,Ct){ #i is index for the parameters
    icf <- (1 - Ct^2)/((1-2*cos(d)*Ct+Ct^2)^(3/2))
    return(icf)
  }
  
  Ct1 <- Ct(1);  Ct2 <- Ct(2);  Ct3 <- Ct(3)
  
  distPQ <- SphDist(lat1,long1, lat2, long2)
  distPtau <- SphDist(lat1,long1,tau[1],tau[2])
  distQtau <- SphDist(lat2,long2,tau[1],tau[2])
  
  icf1 <- ICF(distPQ, Ct1)
  icf2 <- ICF(0, Ct2)
  icf3 <- ICF(distPtau, Ct3)
  icf4 <- ICF(distQtau, Ct3)
  
  trunc1 <- 1/(4*pi) + 3/(4*pi)*Ct1*distPQ
  trunc2 <- 1/(4*pi) + 3/(4*pi)*Ct2
  trunc3 <- 1/(4*pi) + 3/(4*pi)*Ct3*distPtau
  trunc4 <- 1/(4*pi) + 3/(4*pi)*Ct3*distQtau
  
  Nil2 <- 0
  Nil3 <- sum(low_spherical2)
  Nil4 <- sum(low_spherical1)
  
  for(i in 1:length(low_spherical1)){
    for(j in 1:length(low_spherical2)){
     Nil2 <- Nil2 + low_spherical1[i]*low_spherical2[j]
    }
  }
  
  R <- (a3[1]*(icf1 - trunc1)) + (a3[2] * (icf2 - trunc2) * Nil2) + (a3[3] * (icf3 - trunc3) * Nil3) + (a3[3] * (icf4 - trunc4) * Nil4)
  
  ##new version
  #R <- a3[1]*ICF(SphDist(lat1,long1, lat2, long2),1) + a3[2]*1/(4*pi)*ICF(0,2) + a3[3]*1/(2*sqrt(pi))*ICF(SphDist(lat1,long1,tau[1],tau[2]),3) + a3[3]*1/(2*sqrt(pi))*ICF(SphDist(lat2,long2,tau[1],tau[2]),3) - 1/(4*pi) - 1/(16*pi^2) - 1/(4*pi^(3/2))
  #R <- a3[1]*ICF(SphDist(lat1,long1, lat2, long2),1) + a3[2]/(4*pi)*ICF(0,2) + a3[3]/(2*sqrt(pi))*ICF(SphDist(lat1,long1,tau[1],tau[2]),3) + a3[3]/(2*sqrt(pi))*ICF(SphDist(lat2,long2,tau[1],tau[2]),3) - a3[1]/(4*pi) - a3[2]/(16*pi^2) - a3[3]/(4*pi^(3/2))
  
  return(R)
}

#Function to create covariance matrix
Cov_mat <- function(dat, a1, a2, a3, kappa, tau_lat, tau_long, low_spherical, tau_spherical){
  N <- nrow(dat)
  R <- matrix(nrow=N,ncol=N)
  for(i in 1:N){
    for(j in i:N){
      #####change Cov_fn here###################
      if(kappa == 0){
        R[i,j] <- R[j,i] <- Cov_fn00(dat[i,], dat[j,], a1, a2, a3, tau=c(tau_lat, tau_long))
      }else if(kappa == 1){
        R[i,j] <- R[j,i] <- Cov_fn10(dat[i,], dat[j,], a1, a2, a3, tau=c(tau_lat, tau_long))
      }else if(kappa == 2){
        names(low_spherical) <- NULL
        low_spherical <- as.matrix(low_spherical)
        R[i,j] <- R[j,i] <- Cov_fn20(dat[i,], dat[j,], a1, a2, a3, tau=c(tau_lat, tau_long), low_spherical[i,], low_spherical[j,], tau_spherical)
      }
    }
  }
  return(R)
}

#log-likelihood function
loglike <- function(dat, Z, a1, a2, a3, mu.vec=rep(0,nrow(dat))){
  Sigma <- R <- Cov_mat(dat, a1, a2, a3)
  d <- nrow(dat)
  
  sum_i <- sum(apply(t(Z), 1, function(x) t(x-mu.vec)%*%solve(Sigma)%*%(x-mu.vec)))
  #sum_i <- t(matrix(Z))%*%solve(Sigma)%*%matrix(Z)
  logl <- as.numeric(-.5*d*log(2*pi) - .5*log(det(Sigma)) - .5*sum_i)
  
  return(logl)
}


#Function for matrix of spherical distance
Dfun <- function(Lmat, P, t){
  # returns vector of all distance pairs
  
  Dmat <- matrix(NA,P,P)
  
  for (i in 1:P) {
    for (j in i:P) {
      Dmat[i,j] <- Dmat[j,i] <- SphDist(Lmat[i,1],Lmat[i,2],Lmat[j,1],Lmat[j,2])
    }
  }
  
  Dmat <- do.call("cbind", replicate(t, Dmat, simplify = FALSE)) #all the repeated column values of spherical distances (by time)
  Dmat <- do.call("rbind", replicate(t, Dmat, simplify = FALSE)) #all the repeated row values of spherical distances (by time)
  
  
  #Dmat2 <- Dmat
  #for(k in 1:(t-1)){
  #  Dmat2 <- cbind(Dmat, Dmat2) #all the repeated column values of spherical distances (by time)
  #}
  #Dmat3 <- Dmat2
  #for(l in 1:(t-1)){
  #  Dmat3 <- rbind(Dmat2, Dmat3) #all the repeated row values of spherical distances (by time)
  #}
  #Dmat <- Dmat3
  
  Dvec  <- Dmat[upper.tri(Dmat,diag = TRUE)] # length(which(Dvec == 0))
  
  return(Dvec)
}

## To get MoM of the truncated process
order.approx <- function(data, Lmat, Dvec, Hvec, eps, P){
  
  # data must be one time period of data
  # Lmat is a location matrix
  # must have the Legendre polynomial functions loaded
  # Dvec: vectors of distances between points
  # Hvec: bins for distances
  # eps: used to bin distances (how far apart distances can be to be in same bin)
  
  
  outcome <- "X.0"
  
  x1.vars <- "1"
  x2.vars <- paste(x1.vars, "apply(Lmat,1,Y1n1) + apply(Lmat,1,Y10) + 
                    apply(Lmat,1,Y11)", sep = "+")
  x3.vars <- paste(x2.vars, "apply(Lmat,1,Y2n2) + apply(Lmat,1,Y2n1) + 
                    apply(Lmat,1,Y20) + apply(Lmat,1,Y21) + apply(Lmat,1,Y22)",
                   sep = "+")
  x4.vars <- paste(x3.vars, "apply(Lmat,1,Y3n3) + apply(Lmat,1,Y3n2) + 
                    apply(Lmat,1,Y3n1) + apply(Lmat,1,Y30) + apply(Lmat,1,Y31) + 
                    apply(Lmat,1,Y32) + apply(Lmat,1,Y33)", sep = "+")
  x5.vars <- paste(x4.vars, "apply(Lmat,1,Y4n4) + apply(Lmat,1,Y4n3) + 
                    apply(Lmat,1,Y4n2) + apply(Lmat,1,Y4n1) + 
                    apply(Lmat,1,Y40) + apply(Lmat,1,Y41) + apply(Lmat,1,Y42) + 
                    apply(Lmat,1,Y43) + apply(Lmat,1,Y44)", sep = "+")
  x6.vars <- paste(x5.vars, "apply(Lmat,1,Y5n5) + apply(Lmat,1,Y5n4) + 
                    apply(Lmat,1,Y5n3) + apply(Lmat,1,Y5n2) + 
                    apply(Lmat,1,Y5n1) + apply(Lmat,1,Y50) + 
                    apply(Lmat,1,Y51) + apply(Lmat,1,Y52) + apply(Lmat,1,Y53) + 
                    apply(Lmat,1,Y54) + apply(Lmat,1,Y55)", sep = "+")
  x7.vars <- paste(x6.vars, "apply(Lmat,1,Y6n6) + apply(Lmat,1,Y6n5) + 
                    apply(Lmat,1,Y6n4) + apply(Lmat,1,Y6n3) + 
                    apply(Lmat,1,Y6n2) + apply(Lmat,1,Y6n1) + 
                    apply(Lmat,1,Y60) + apply(Lmat,1,Y61) + apply(Lmat,1,Y62) + 
                    apply(Lmat,1,Y63) + apply(Lmat,1,Y64) + apply(Lmat,1,Y65) + 
                    apply(Lmat,1,Y66)", sep = "+")
  
  
  x1.form <- as.formula(paste(outcome, x1.vars, sep = "~"))
  x2.form <- as.formula(paste(outcome, x2.vars, sep = "~"))
  x3.form <- as.formula(paste(outcome, x3.vars, sep = "~"))
  x4.form <- as.formula(paste(outcome, x4.vars, sep = "~"))
  x5.form <- as.formula(paste(outcome, x5.vars, sep = "~"))
  x6.form <- as.formula(paste(outcome, x6.vars, sep = "~"))
  x7.form <- as.formula(paste(outcome, x7.vars, sep = "~"))
  
  
  
  X.0 <- data
  X.1 <- lm(x1.form)$res # order 1
  X.2 <- lm(x2.form)$res # order 2
  X.3 <- lm(x3.form)$res # order 3
  X.4 <- lm(x4.form)$res # order 4
  X.5 <- lm(x5.form)$res # order 5
  X.6 <- lm(x6.form)$res # order 6
  X.7 <- lm(x7.form)$res # order 7
  
  
  G.0 <- G_Hat(X.0, Dvec, Hvec, eps, P)
  G.1 <- G_Hat(X.1, Dvec, Hvec, eps, P)
  G.2 <- G_Hat(X.2, Dvec, Hvec, eps, P)
  G.3 <- G_Hat(X.3, Dvec, Hvec, eps, P)
  G.4 <- G_Hat(X.4, Dvec, Hvec, eps, P)
  G.5 <- G_Hat(X.5, Dvec, Hvec, eps, P)
  G.6 <- G_Hat(X.6, Dvec, Hvec, eps, P)
  G.7 <- G_Hat(X.7, Dvec, Hvec, eps, P)
  
  G <- t(rbind(G.0$G,G.1$G,G.2$G,G.3$G,G.4$G,G.5$G,G.6$G,G.7$G))
  N <- t(rbind(G.0$N,G.1$N,G.2$N,G.3$N,G.4$N,G.5$N,G.6$N,G.7$N))
  
  return(list(G = G, N = N))
}


#Function for MOM estimator
G_Hat <- function(dat, sDvec, Hvec, tDvec ,eps, P, t_lag, kappa, low_spherical){
  
  if(kappa==0){ #no truncation required
    Zr <- dat$Z
  }else if(kappa==1){ #get truncuated process by residuals
    #order 1
    Zr <- lm(Z~1+time, dat)$res
  }else if(kappa==2){
    Zr <- lm(dat$Z ~ 1 + low_spherical[,2] + low_spherical[,3] + low_spherical[,4] + as.factor(dat$time))$res
  }

  # Create Second Moment Matrix and Vector
  ZZmat <- matrix(0,P*t,P*t)
  
  ZZmat <- Zr%*%t(Zr)
  #  for (i in 1:P) {
  #    for (j in i:P) {
  #      ZZmat[i,j] <- Z[i]*Z[j]
  #    }
  #  } # length(which(ZZmat == 0)) # (P - 1)*P/2
  
  ZZvec <- ZZmat[upper.tri(ZZmat,diag = TRUE)] # (P + 1)*P/2
  
  # Calculate Estimates using MOM Estimator
  G <- N <- matrix(nrow=length(Hvec), ncol=t_lag)
  
  for(i in 1:t_lag){
    
    for (j in 1:(length(Hvec))){
      s <- which(((Hvec[j] - eps) <  sDvec) & (sDvec <  (Hvec[j] + eps)) & tDvec == i-1)
      #s <- which((sDvec > Hvec[j]) & (sDvec <=  Hvec[j]+eps) & tDvec == i)
      N[j,i] <- length(s)
      G[j,i] <- mean(ZZvec[s])# - mean(Z)^2
    }
  }
  
  return(list(G = G, N = N))
}

#Function to compute theoretical covariance values
Rvec <- function(dist, time, a1, a2, a3, kappa){
  Ct <- a1*exp(-a2*abs(time))
  R <- a3*(1 - Ct^2)/(1-2*cos(dist)*Ct+Ct^2)^(3/2)

  if(kappa >= 1) R <- R - a3 * 1/(4*pi)
  if(kappa >= 2) R <- R - a3 * Ct * 3/(4*pi) * P.1(cos(dist))
  #if(kappa >= 3) R <- R - a3 * 5*(Ct^2) / (4*pi) * P.2(cos(dist))
  return(R)
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
    for(j in 0:(t/2-1)){
      R_mat[i,j+1] <- Rvec(dist=Hvec[i], time=j ,a1=a[1], a2=a[2], a3=a[3])
    }
  }
  sum((R_mom- R_mat)^2) #without the weights
  #sum(W*(R_mom- R_mat)^2)
}


############################# Plots for the Simulation Data ###################################
#Plots for "fixed" ylim version
#plot for spherical distances
s_plot <- function(j,R_mat,R_fitted,R_mom){
  #Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  #Hvec <- as.vector(seq(0, 2, by=2/10)[1:10])
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  min_r <- min(R_mat[,1], R_fitted[,1], R_mom[,1])
  max_r <- max(R_mat[,1], R_fitted[,1], R_mom[,1])
  plot1 <- plot(Hvec, R_mat[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("h=",(j-1)), ylim=c(0.8*min_r, 1.1*max_r)); 
  lines(Hvec, R_fitted[,j], col="red", lty=2); lines(Hvec, R_mom[,j], col="blue", lty=4)
  if(j==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}

#plot for time
t_plot <- function(i,R_mat,R_fitted,R_mom){
  #Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  #Hvec <- as.vector(seq(0, 2, by=2/10))[1:10]
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  time <- 0:(ncol(R_mat)-1)
  min_r <- min(R_mat[1,], R_fitted[1,], R_mom[1,])
  max_r <- max(R_mat[1,], R_fitted[1,], R_mom[1,])
  plot1 <- plot(time, R_mat[i,], type="l", lty=1, xlab="Time Lag", ylab="Covariance", main=bquote(psi== .(Hvec[i])), ylim=c(0.8*min_r, 1.1*max_r)); 
  lines(time, R_fitted[i,], col="red", lty=2); lines(time, R_mom[i,], col="blue", lty=4)
  if(i==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}

#Plots for "adjusted" ylim version
#plot for spherical distances
s_plot2 <- function(j,R_mat,R_fitted,R_mom){
  #Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  #Hvec <- as.vector(seq(0, 2, by=2/10))[1:10]
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  min_r <- min(R_mat[,j], R_fitted[,j], R_mom[,j])
  max_r <- max(R_mat[,j], R_fitted[,j], R_mom[,j])
  plot1 <- plot(Hvec, R_mat[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("h=",(j-1)), ylim=c(0.8*min_r, 1.1*max_r)); 
  lines(Hvec, R_fitted[,j], col="red", lty=2); lines(Hvec, R_mom[,j], col="blue", lty=4)
  if(j==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}

#plot for time
t_plot2 <- function(i,R_mat,R_fitted,R_mom){
  #Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  #Hvec <- as.vector(seq(0, 2, by=2/10))[1:10]
  Hvec <- as.vector(seq(0, 2, by = 1/10))
  time <- 0:(ncol(R_mat)-1)
  min_r <- min(R_mat[i,], R_fitted[i,], R_mom[i,])
  max_r <- max(R_mat[i,], R_fitted[i,], R_mom[i,])
  plot1 <- plot(time, R_mat[i,], type="l", lty=1, xlab="Time Lag", ylab="Covariance", main=bquote(psi== .(Hvec[i])), ylim=c(0.8*min_r, 1.1*max_r)); 
  lines(time, R_fitted[i,], col="red", lty=2); lines(time, R_mom[i,], col="blue", lty=4)
  if(i==1){
    legend("topright", legend=c("True", "Fitted", "MoM"), lty=c(1,2,4), col=c("black","red","blue"))
  }
}


################################### Contour Plot #############################################

contour_plot <- function(dat, type){
  colnames(dat) <- 0:9
  row.names(dat) <- as.vector(seq(0, 2, by = 1/10))
  
  df <- reshape2::melt(dat)
  
  # time <- 0:9
  # Hvec <- as.vector(seq(0, 2, by = 1/10))
  # contour_dat <- data.frame(y=as.vector(result$R_fitted), time=rep(time,each=nrow(result$R_fitted)), dist=rep(Hvec,ncol(result$R_fitted)))c
  
  if(type==1){
    p <- ggplot2::ggplot(df, aes(Var1, Var2, z= value, colour=stat(level))) +
      geom_contour() + xlab("Spherical Distance") + ylab("Time Lag") + ggtitle("True Value") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_distiller(palette = "YlGn", direction = 1)
  }else if(type==2){
    p <- ggplot2::ggplot(df, aes(Var1, Var2, z= value, colour=stat(level))) +
      geom_contour() + xlab("Spherical Distance") + ylab("Time Lag") + ggtitle("Fitted Value") +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_distiller(palette = "YlGn", direction = 1)
  }else if(type==3){
    p <- ggplot2::ggplot(df, aes(Var1, Var2, z= value, colour=stat(level))) +
      geom_contour() + xlab("Spherical Distance") + ylab("Time Lag") + ggtitle("MoM") + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_colour_distiller(palette = "YlGn", direction = 1)
  }
  
  
  # p <- ggplot2::ggplot(df, aes(Var1, Var2, z= value)) +
  #   stat_contour(geom="polygon",aes(fill=stat(level))) +
  #   scale_fill_distiller(palette = "Spectral", direction = -1)
  
  # p <- ggplot(df, aes(Var1, Var2, z= value)) +
  #   geom_contour() +
  #   scale_fill_distiller(palette = "Spectral", direction = -1)
  
  #return(plotly::ggplotly(p))
}



############################# Plots for Temperature Anomaly Data ###################################

#plot for spherical distances
temp_s_plot <- function(j, R_mom, R_fitted, ylim_upper=06){
  Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  plot1 <- plot(Hvec, R_mom[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("h=",(j-1)), ylim=c(0, ylim_upper)); 
  lines(Hvec, R_fitted[,j], col="red", lty=2)
  if(j==1){
    legend("topright", legend=c("MoM","fitted"), lty=1:2, col=c("black","red"))
  }
}

#plot for time
temp_t_plot <- function(i, R_mom, R_fitted, ylim_upper=0.6){
  Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  time <- 0:(ncol(R_fitted)-1)
  plot1 <- plot(time, R_mom[i,], type="l", lty=1, xlab="Time lag", ylab="Covariance", main=bquote(psi== .(Hvec[i])), ylim=c(0, ylim_upper));
  lines(time, R_fitted[i,], col="red", lty=2)
  if(i==1){
    legend("topright", legend=c("MoM","fitted"), lty=1:2, col=c("black","red"))
  }
}


#plot for spherical distances
temp_s_plot2 <- function(j, R_fitted, ylim_upper=0.6){
  Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  #plot1 <- plot(Hvec, R_fitted[,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("t=",(j-1)), ylim=c(0, 1.5));
  plot1 <- plot(Hvec, R_fitted[[1]][,j], type="l", lty=1, xlab="Spherical Distance", ylab="Covariance", main=paste("h=",(j-1)), ylim=c(0, ylim_upper)); 
  lines(Hvec, R_fitted[[2]][,j], col="blue" ,lty=2);
  lines(Hvec, R_fitted[[3]][,j], col="red", lty=3);
  lines(Hvec, R_fitted[[4]][,j], col="green", lty=4)
  if(j==1){
    legend("topright", legend=c("1980","1990", "2000","2010"), lty=1:4, col=c("black","blue","red","green"))
  }
}

#plot for time
temp_t_plot2 <- function(i, R_fitted, ylim_upper=0.6){
  Hvec <- as.vector(seq(0, 1, by = 1/10))[1:10]
  time <- 0:(ncol(R_fitted[[1]])-1)
  #plot1 <- plot(time, R_fitted[i,], type="l", lty=1, xlab="Time lag", ylab="Covariance", main=paste(expression(theta),"=",(i-1)), ylim=c(0, 1.5));
  plot1 <- plot(time, R_fitted[[1]][i,], type="l", lty=1, xlab="Time lag", ylab="Covariance", main=bquote(psi== .(Hvec[i])), ylim=c(0, ylim_upper));
  lines(time, R_fitted[[2]][i,], col="blue", lty=2);
  lines(time, R_fitted[[3]][i,], col="red", lty=3);
  lines(time, R_fitted[[4]][i,], col="green", lty=4)
  if(i==1){
    legend("topright", legend=c("1980","1990", "2000","2010"), lty=1:4, col=c("black","blue","red","green"))
  }
}



#################################### Cpp version ####################################
#Wrapper function of C codes
C_Cov_mat <- function(dat=dat, sDmat, a1, a2, a3, kappa, low_spherical, tau_spherical){
  
  dat$index = dat$index-1 #for C
  nrow_dat <- nrow(dat); ncol_dat <- ncol(dat) #data set
  nrow_Dmat <- nrow(sDmat) #spherical distance matrix
  ncol_low_spherical <- ncol(low_spherical) # the number of truncated spherical harmonics
  
  R <- matrix(0, nrow=nrow(dat), ncol=nrow(dat)) #covariance matrix
  
  res <- .C('Cov_mat', as.integer(nrow_dat),  as.integer(ncol_dat), as.double(as.matrix(dat)), as.integer(nrow_Dmat),
            as.double(sDmat), as.double(a1), as.double(a2), as.double(a3), as.integer(kappa), as.double(as.matrix(low_spherical)), as.integer(ncol_low_spherical), as.double(tau_spherical), R=as.double(R))
  
  return(matrix(res$R,nrow=nrow_dat))
}


#Wrapper function of C codes
C_G_hat <- function(dat=dat, sDmat, Hvec=Hvec, eps=eps, P=P, t=t_lag, kappa, low_spherical){
  dat$index = dat$index-1 #for C
  nrow_dat <- nrow(dat); ncol_dat <- ncol(dat) #data set
  nrow_Dmat <- nrow(sDmat) #spherical distance matrix
  
  G <- N <- matrix(0, nrow=length(Hvec), ncol=t)
  
  if(kappa==1){
    #dat$Z <- lm(Z~1,dat)$res #get the truncated process
    dat$Z <- lm(Z~1+time,dat)$res #get the truncated process given t
  }else if(kappa==2){
    dat$Z <- lm(dat$Z ~ 1 + low_spherical[,2] + low_spherical[,3] + low_spherical[,4] + as.factor(dat$time))$res
  }
  
  res <- .C('G_hat', as.integer(nrow_dat),  as.integer(ncol_dat), as.double(as.matrix(dat)), as.integer(nrow_Dmat),
            as.double(sDmat), as.integer(length(Hvec)), as.double(Hvec), as.double(eps), as.integer(t), G=as.double(G), N=as.double(N))
  
  #void G_hat3(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, int* length_Hvec, double* Hvec, double* eps, int* t, double* G, double* N)
  
  return(list(res$G, res$N))
}

#Function to get Spherical Distance Matrix
C_SphDist <- function(dat){
  nrow_dat <- nrow(dat); ncol_dat <- ncol(dat) #data set
  sDmat <- matrix(0, nrow=nrow_dat, ncol=nrow_dat) #spherical distance matrix

  res <- .C('SphDist',as.integer(nrow_dat),  as.integer(ncol_dat), as.double(as.matrix(dat)), sDmat=as.double(sDmat))
  
  return(matrix(res$sDmat,nrow=nrow_dat))
}

######################################################################################