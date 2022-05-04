#' An average cluster size calculation for AHMM with 2-level clustering
#'
#' This function performs a cluster size calculation for an underlying test of
#' treatment effect under the Additive Hazards Mixed Model with two-level
#' clustering and randomization occurring at the higher (cluster) level.
#' The underlying methods are described in (2022 Blaha, Esserman, Li paper).
#'
#'
#' @param lambda0 Anticipated baseline hazard (scalar)
#' @param delta Treatment effect on hazard difference scale (scalar)
#' @param sigma2 Clustering effect - variance of the frailty parameter
#' @param tau Duration of the planned study (administrative censoring)
#' @param n Anticipated sample size
#' @param CV Coefficient of variation among cluster sizes. Default is 0
#' @param pi Randomization ratio (arm B vs. A). Accepted values are between 0
#' and 1. Larger values than 0.5 indicated that more subjects are allocated
#' into the experimental arm. Smaller values than 0.5 indicated that more
#' subjects are allocated into the control arm. Default is 0.5.
#' @param alpha Desired size of the underlying test. Default is alpha=0.05.
#' @param beta Desired power of the test. Default is 80%.
#' @param test Choice between "z" for large sample z-test and "t" for small
#' sample t-test. Default is c("t","z").
#' @return This function returns returns required (average) cluster size of the
#' underlying test of hypothesis H_0:delta=0 vs. two-sided alternative for
#' either t-test when test="t" is selected or z-test if test="z" is selected,
#' given the desired power and size of the test.
#' @keywords methods
#' @export
#' @examples
#' ClusterSize(lambda0=4,delta=1.5,sigma2=0.5,tau=1,n=48)
#'
#'
#'
ClusterSize<-function(lambda0,delta,sigma2,tau,n,CV=0,pi=0.5,alpha=0.05,
                      beta=0.2,test=c("z","t")){

  ####################################################
  #Individual quantities for cluster size calculation#
  ####################################################
  # g1
   g1 = pracma::integral(one,xmin=0,xmax=tau,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi)

  # ICC (Rho tilde)
  rho = (pracma::dblquad(two,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(three,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(four,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(five,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(six,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi)) /
         g1

  # h
    h = (pracma::integral(h,0,tau,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi))

  #################################################################
  #Cluster Size - for given sample size (n) and other input values#
  #     (lambda0, delta, sigma2, tau, n, CV, pi, alpha, beta)     #
  #################################################################
  Q<-(qt(1-beta,df=n-1)-qt(alpha/2,df=n-1))^2
  Z<-(qnorm(1-beta)-qnorm(alpha/2))^2
  m1_init<-Z*(1-rho)/(n*delta^2*(1/g1)*h^2-Z*(1+CV^2)*rho)
  M1<-ceiling(m1_init)
  m2_init<-Q*(1-rho)/(n*delta^2*(1/g1)*h^2-Q*(1+CV^2)*rho)
  M2<-ceiling(m2_init)

  if (sum(test==c("z","t"))==2|sum(test==c("t","z"))==2) {
    return(
      cat("\n",
          "############################","\n",
          "###Required Cluster Size:###","\n",
          "############################","\n",
          "Sample Size for z-test:         ",M1,"\n",
          "Sample Size for t(n-1)-test:    ",M2,"\n","\n"
      )
    )
  } else if (sum(test=="z")==1) {
    return(
      cat("\n",
          "############################","\n",
          "###Required Cluster Size:###","\n",
          "############################","\n",
          "Sample Size for z-test:         ",M1,"\n","\n"
      )
    )
  } else if (sum(test=="t")==1) {
    return(
      cat("\n",
          "############################","\n",
          "###Required Cluster Size:###","\n",
          "############################","\n",
          "Sample Size for t(n-1)-test:    ",M2,"\n","\n"
      )
    )
  } else return("warning: test not properly specified")
}
