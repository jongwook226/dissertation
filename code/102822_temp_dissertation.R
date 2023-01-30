load("/Users/jongwookkim/Downloads/git2/DA_project/C_code_result.RData")

setwd("/Users/jongwookkim/Library/CloudStorage/OneDrive-IndianaUniversity/IU/Dissertation")
source("102822_functions_dissertation.R")
source("irfkrg-fcn01-sphere-harmonics.R")
source("irfkrg-fcn02-legendre-polynomials.R")

load("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Dmat.Rda") #Spherical distance matrix
sDmat <- Dmat

temp_parameters <- function(year=1980, kappa=1){
  #Load data
  temp1 <- readRDS(paste("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-", year, ".Rds", sep=""))
  temp2 <- readRDS(paste("/Users/jongwookkim/OneDrive - Indiana University/IU/DA/Data and Code for Jongwook/EDA Data/irf-eda01V2-parse-data-", year+1, ".Rds", sep=""))
  
  ############################################
  # #sampling
  # index <- sample(1:nrow(temp1), size=700, replace=FALSE)
  # temp1 <- temp1[index,]
  # temp2 <- temp2[index,]
  ############################################
  
  #remove the redundancy of the column names.
  temp <- cbind(temp1, temp2[,3:14])
  colnames(temp)[15:26] <- c(13:24)
  
  #############################################
  ##create a spherical distance maxtrix
  #P <- nrow(temp)
  #sDmat <- matrix(NA,P,P)
  #for (i in 1:P) {
  #  for (j in i:P) {
  #    sDmat[i,j] <- sDmat[j,i] <- SphDist(temp[i,1],temp[i,2],temp[j,1],temp[j,2])
  #  }
  #}
  ##############################################
  P <- nrow(temp) #location data size
  t <- ncol(temp)-2 #time data size
  n <- P*t #total data size
  
  time <- seq(1, ncol(temp)-2) #time term
  #t <- 1:5
  lat <- rep(temp[,1], length(t))
  long <- rep(temp[,2], length(t))
  time <- rep(time, each=nrow(temp))
  
  dat <- data.frame(lat,long)
  
  Z <- as.vector(temp[,c(-1,-2)])
  index <- rep(1:nrow(temp), length(t)) #Create location index
  
  
  
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
  
  
  ####################### Do we need tau???? ######################
  tau_lat <- runif(1, min=-pi/2, max=pi/2) #tau_latitude
  tau_long <- runif(1, min=0, max=2*pi) #tau_longitude
  
  
  #Spherical distance with tau
  tauDist <- rep(0, P)
  if(kappa != 0){ #we do not need tau if kappa =0
    for(i in 1:nrow(dat)){
      tauDist[i] <- SphDist(dat[i,1], dat[i,2], tau_lat, tau_long) 
    }
  }
  
  #################################################################
  
  dat <- data.frame(lat,long,time,tauDist,index,Z)
  
  
  ####################################################################
  ######################Compute MOM estimator#########################
  ####################################################################
  
  H <- 10         # Number of distances at which to estimate
  Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
  eps <- 0.10     # use to bin distances in estimation
  
  
  ####################Use a C function##########################
  
  dyn.load("C_dissertation.so")
  #is.loaded('SphDist')
  #is.loaded('Cov_mat')
  #is.loaded('G_hat')
  #dyn.unload("C_dissertation.so")
  
  # t_lag <- ceiling(t/2) #time lag
  t_lag <- 10
  
  mom_est <- C_G_hat(dat=dat, sDmat, Hvec=Hvec, eps=eps, P=P, t_lag=t_lag, kappa=kappa, low_spherical)
  R_mom <- matrix(mom_est[[1]]/mom_est[[2]], nrow=length(Hvec)) #MOM estimate
  W <- matrix(mom_est[[2]], nrow=length(Hvec)) #weight
  
  
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
  
  #parameters <- nlminb(c(0,0,0),obj)$par #fitted parameters by nlminb
  fitted_para <- suppressWarnings(nlminb(c(0,0,0),obj,lower=c(0,0,0), upper=c(1,Inf,Inf)))$par #fitted parameters by nlminb
  
  cat(year, fitted_para, "\n\n")
  
  #fitted covariance matrix
  R_fitted <- matrix(ncol=ncol(R_mom), nrow=nrow(R_mom))
  for(i in 1:length(Hvec)){
    for(j in 1:t_lag){
      R_fitted[i,j] <- Rvec(dist=Hvec[i], time=j-1, fitted_para[1], fitted_para[2], fitted_para[3], kappa)
    }
  }
  
  result <- list(R_fitted=R_fitted, R_mom=R_mom, para_est=fitted_para) #retruns fitted, Mom covariance matrices, and parameter estimates.
  
  
  return(result)
}



years <- sapply(c(1:16),function(x) 2*x+1978)
para <- matrix(NA, ncol=3, nrow=length(years))
row.names(para) <- years

result <- list(R_fitted=matrix(rep(NA,100),nrow=10), R_mom= matrix(rep(NA,100),nrow=10), para_est=rep(NA,3))
res_years <- list()
for(i in 1:length(years)) res_years [[i]] <- result
names(res_years) <- paste("y", sapply(c(1:16),function(x) 2*x+1978), sep="")


for(i in 1:length(years)){
  res_years[[i]] <- temp_parameters(year=years[i], kappa=1)
  para[i,] <- res_years[[i]]$para_est
  #para[i,] <- temp_parameters(year=years[i], kappa=1)$para_est
  #save(para, file="parameters.Rda")
}

#save(res_years,para, file="temp_analysis_kappa1.Rda")




################################# RUN FROM HERE!!!!!! ######################################
################################# RUN FROM HERE!!!!!! ######################################
################################# RUN FROM HERE!!!!!! ######################################


setwd("/Users/jongwookkim/Library/CloudStorage/OneDrive-IndianaUniversity/IU/Dissertation")
load("temp_analysis_kappa1.Rda")
source("102822_functions_dissertation.R")
years <- sapply(c(1:16),function(x) 2*x+1978)

##################### Compare the Values of the Parameter estimators #####################

par(mfrow=c(1,1))
plot(years, para[,1], ylab="p1"); lines(years, para[,1]) #p1
plot(years, para[,2], ylab="p2"); lines(years, para[,2]) #p2
plot(years, para[,3], ylab="p3"); lines(years, para[,3]) #p3


plot(years, para[,1], ylim=c(0,2), ylab="Parameters"); lines(years, para[,1]); #p1
lines(years, para[,2], col="red", lty=2, type="l");lines(years, para[,2], col="red", lty=2, type="o");#p2
lines(years, para[,3], col="blue", lty=3, type="l"); lines(years, para[,3], col="blue", lty=3, type="o");#p3
legend("topright", legend=c("p1","p2", "p3"), lty=1:3, col=c("black","red","blue"))




#################### plot to compare MoM and Fitted covariance estimators ###################

par(mfrow=c(2,3))

n=16
matrix(c(1:length(years),years), ncol=2)

temp_s_plot(1, res_years[[n]]$R_mom, res_years[[n]]$R_fitted, ylim_upper=0.6)
temp_s_plot(2, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_s_plot(3, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)

temp_t_plot(1, res_years[[n]]$R_mom, res_years[[n]]$R_fitted, ylim_upper = 0.6)
temp_t_plot(2, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_t_plot(3, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)

# R_mom = res_years[[n]]$R_mom
# R_fitted = res_years[[n]]$R_fitted
# 
# temp_s_plot(1, res_years[[n]]$R_mom, res_years[[n]]$R_fitted, ylim_upper=0.6)
# temp_s_plot(2, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_s_plot(3, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_s_plot(4, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_s_plot(5, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_s_plot(6, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_s_plot(7, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_s_plot(8, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_s_plot(9, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_s_plot(10, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# 
# 
# temp_t_plot(1, res_years[[n]]$R_mom, res_years[[n]]$R_fitted, ylim_upper = 0.6)
# temp_t_plot(2, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_t_plot(3, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_t_plot(4, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_t_plot(5, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_t_plot(6, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# temp_t_plot(7, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_t_plot(8, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_t_plot(9, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
# #temp_t_plot(10, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)



############################# Compare by different years #####################################
par(mfrow=c(2,3))

#Choose years for plot (1980, 1990, 2000, 2010)
R_fitted <- list(res_years$y1980$R_fitted, res_years$y1990$R_fitted, res_years$y2000$R_fitted, res_years$y2010$R_fitted)

n=1
matrix(c(1:length(years),years), ncol=2)

temp_s_plot2(1, R_fitted, ylim_upper = 0.6)
temp_s_plot2(2, R_fitted)
temp_s_plot2(3, R_fitted)
#temp_s_plot2(4, R_fitted)
temp_s_plot2(5, R_fitted)
#temp_s_plot2(6, R_fitted)
temp_s_plot2(7, R_fitted)
#temp_s_plot2(8, R_fitted)
temp_s_plot2(9, R_fitted)
#temp_s_plot2(10, R_fitted)

temp_t_plot2(1, R_fitted)
temp_t_plot2(2, R_fitted)
temp_t_plot2(3, R_fitted)
#temp_t_plot2(4, R_fitted)
temp_t_plot2(5, R_fitted)
#temp_t_plot2(6, R_fitted)
temp_t_plot2(7, R_fitted)
#temp_t_plot2(8, R_fitted)
temp_t_plot2(9, R_fitted)
#temp_t_plot2(10, R_fitted)





temp_s_plot2(2, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_s_plot2(3, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
#temp_s_plot(4, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
#temp_s_plot(5, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_s_plot(6, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
#temp_s_plot(7, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_s_plot(8, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
#temp_s_plot(9, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)
temp_s_plot(10, res_years[[n]]$R_mom, res_years[[n]]$R_fitted)