#Testing Examples - Data generation, short AHMM simulation, permutation matrix generation
require(pracma)
require(Rfast)
require(rpsftm)
require(survival)
source("AHMM/R/Integrals.R")
##########################
#Fixed Initial Parameters#
##########################
lambda0<-4
delta<-1.5
sigma2<-0.5
pi<-1/2
CV<-0.25
m_bar<-10
tau<-1
nsim<-10
N1<-48


#################
#Generating Data#
#################
delta_init<-delta
nsim_init <- nsim
n <- N1
cl <- 1:n
ndelta<-length(delta)
DELTA=matrix(NA,nsim,ndelta)
neg_haz_powerz=0
n_cl_powerz=rep(NA,nsim)
nk <- matrix(NA,n,nsim)
est <- c()

for(s in 1:nsim){
  set.seed(395687+s^2)
  skip<-0

  # Cluster level information
  xi <- rnorm(n,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, n)
  trt[sample(1:n,n/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+ndelta))

  # Check if hazard is negative
  for (k in 1:n){
    if (trt[k]*delta+xi[k]+lambda0<1e-5)
    {neg_haz_powerz=neg_haz_powerz+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Zk <- rep(trt[k],nk[k,s])
      Tk <- -log(runif(nk[k,s]))/c(Zk*delta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Zk,pmin(Tk,Ck),Tk<=Ck))
    }
  }

  ID <- data[,1]
  Z <- data[,1+(1:ndelta),drop=FALSE]
  U <- data[,2+ndelta]
  Delta <- data[,ncol(data)]

  Test_Data<-data.frame(ID,U,Z,Delta)
  beta_est<-AHMM(Surv(U,Delta)~Z, cluster = ID, data = Test_Data) #Testing AHMM function
  est<-rbind(est,beta_est)
}
#delta_bar
#mean(est)




####################################################
#Generating Treatment Allocation Permutation Matrix#
####################################################
if(choose(n, n/2)<=50000){
  # Enumeration of all possible schemes
  emr <- t(combn(n, n / 2))
  R <- dim(emr)[1]
  pmt_test <- matrix(0, R, n)
  for (r in 1:R){
    pmt_test[r, emr[r,]] <- 1
  }
} else{
  S=50000
  # Simulation of permutations
  pmt_test <- matrix(0, S, n)
  for (s in 1:S){
    trt_perm <- sample(cl, n / 2)
    pmt_test[s, trt_perm] <- 1
  }
  pmt_test <- unique(pmt_test) # Indicator matrix
  R <- dim(pmt_test)[1]
}

#Testing PermTest function
PermTest(ID=Test_Data$ID, Z_cl=Test_Data$Zk, U=Test_Data$U, Delta=Test_Data$Delta, pmt = pmt_test)




#Testing Example - AHMM with covariates
rm(list=ls(all=TRUE))
require(pracma)
require(Rfast)
require(rpsftm)
require(survival)
setwd("C:/Users/ob246/Desktop/AHMM_Package/AHMM")
source("R/Integrals.R")
##########################
#Fixed Initial Parameters#
##########################
lambda0<-4
delta<-1.5
beta<--2
theta<-c(delta,beta)  ##treatment effect, one extra cont. variable
sigma2<-0.5
pi<-1/2
CV<-0.25
m_bar<-10
tau<-1
nsim<-100  ##Let's leave the extra dimension here in case we need to generate more datasets at once
N1<-56


#################
#Generating Data#
#################
cl <- 1:N1
npar<-length(theta)
DELTA=matrix(NA,nsim,npar)
neg_haz=0
n_cl=rep(NA,nsim)
nk <- matrix(NA,N1,nsim)
d<-matrix(0,2,1)

for(s in 1:nsim){
  set.seed(39568+s^2)
  skip<-0

  # Cluster level information
  xi <- rnorm(N1,mean=0,sd=sqrt(sigma2))
  trt <- rep(0, N1)
  trt[sample(1:N1,N1/2,replace=F)] <- 1
  data <- matrix(nrow = 0, ncol = (3+npar))

  # Check if hazard is negative
  for (k in 1:N1){
    if (trt[k]*delta+xi[k]+lambda0<1e-5)   ##Does not include Xij
    {neg_haz=neg_haz+1
    skip<-1}
    else{
      if (CV==0) {temp<-m_bar}
      else {
        temp<-ceiling(rgamma(1,shape=CV^-2,rate=m_bar^-1*CV^-2))}
      nk[k,s] <- max(2,temp)
      Zi <- rep(trt[k],nk[k,s])
      #Zi <- as.vector(pracma::rand(nk[k,s],1)<0.5)
      Wij <- as.vector(pracma::rand(nk[k,s],1))
      Xij <- cbind(Zi,Wij)
      Tk <- -log(runif(nk[k,s]))/(Xij%*%theta+xi[k]+lambda0)
      Ck <- runif(nk[k,s])*tau
      data <- rbind(data,cbind(repmat(k,nk[k,s],1),Xij,pmin(Tk,Ck),Tk<=Ck))
    }
  }

  ID <- data[,1]
  X <- data[,1+(1:npar),drop=FALSE]
  U <- data[,2+npar]
  Delta <- data[,ncol(data)]

  Test_Data<-data.frame(ID,U,X,Delta)

  d<-d+AHMM(Surv(U,Delta)~Zi + Wij, cluster = ID, data = Test_Data)
}
d/nsim
#Test_Data<-data.frame(ID,U,X,Delta)
#head(Test_Data)

#AHMM_W(Surv(U,Delta)~Zi + Wij, cluster = ID, data = Test_Data)
