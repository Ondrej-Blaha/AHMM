#' A Sample size calculation for AHMM with 2-level clustering
#'
#' This function performs sample size calculation for an underlying test of
#' treatment effect under the Additive Hazards Mixed Model with two-level
#' clustering and randomization occurring at the higher (cluster) level.
#' The underlying methods are described in (2022 Blaha, Esserman, Li paper).
#'
#'
#'
#' @param lambda0 Anticipated baseline hazard (scalar)
#' @param delta Treatment effect on hazard difference scale (scalar)
#' @param sigma2 Clustering effect - variance of the frailty parameter
#' @param tau Duration of the planned study (administrative censoring)
#' @param m_bar Anticipated average cluster size
#' @param CV Coefficient of variation among cluster sizes. Default is 0.
#' @param pi Randomization ratio (arm B vs. A). Accepted values are between 0
#' and 1. Larger values than 0.5 indicated that more subjects are allocated
#' into the experimental arm. Smaller values than 0.5 indicated that more
#' subjects are allocated into the control arm. Default is 0.5.
#' @param alpha Desired size of the underlying test. Default is alpha=0.05.
#' @param beta Desired power of the test. Default is 80%.
#' @param test Choice between "z" for large sample z-test and "t" for small
#' sample t-test. Default is c("t","z").
#' @return This function returns required sample size of the underlying
#' test of hypothesis H_0:delta=0 vs. two-sided alternative for either t-test
#' when test="t" is selected or z-test if test="z" is selected, given the
#' desired power and size of the test.
#' @keywords methods
#' @export
#' @examples
#' SampleSize(lambda0=4,delta=1.5,sigma2=0.5,tau=1,m_bar=10)
#'
#'
#'
SampleSize<-function(lambda0,delta,sigma2,tau,m_bar,CV=0,pi=0.5,alpha=0.05,
                     beta=0.2,test=c("z","t")){

###################################################
#Individual quantities for sample size calculation#
###################################################
# Meat of VBV "B"
B=(pracma::integral(one,xmin=0,xmax=tau,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
     (m_bar*(1+CV^2)-1)*(pracma::dblquad(two,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
                         pracma::dblquad(three,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
                         pracma::dblquad(four,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
                         pracma::dblquad(five,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
                         pracma::dblquad(six,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi)))
# Bread of VBV "V"
V=sqrt(m_bar)*(pracma::integral(h,0,tau,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi))

# Sandwich Variance
sigma2_delta <- solve(V)%*%B%*%solve(V)

########################################################
#Sample Size - for given power, delta, and sigma2_delta#
########################################################
z<-((qnorm(1-beta)+qnorm(1-alpha/2))^2/delta^2)*(sigma2_delta)
N1<-ceiling(as.numeric(z))+pracma::mod(ceiling(as.numeric(z)),2)
powerz<-pnorm(sqrt(N1)*abs(delta)/sqrt(sigma2_delta)-qnorm(1-alpha/2),mean=0,sd=1)

n_init<-N1
while (n_init<((qt(1-beta,df=n_init-1)-qt(alpha/2,df=n_init-1))^2/delta^2)*(sigma2_delta))
{n_init<-n_init+1}
N2<-n_init+pracma::mod(n_init,2)
#powert<-pt(sqrt(N2)*abs(delta)/sqrt(sigma2_delta)-qt(1-alpha/2,df=N2-1),df=N2-1)

if (sum(test==c("z","t"))==2|sum(test==c("t","z"))==2) {
  return(
    cat("\n",
        "###########################","\n",
        "###Required Sample Size:###","\n",
        "###########################","\n",
        "Sample Size for z-test:         ",N1,"\n",
        "Sample Size for t(n-1)-test:    ",N2,"\n","\n"
    )
  )
} else if (sum(test=="z")==1) {
  return(
    cat("\n",
        "###########################","\n",
        "###Required Sample Size:###","\n",
        "###########################","\n",
        "Sample Size for z-test:         ",N1,"\n","\n"
    )
  )
} else if (sum(test=="t")==1) {
  return(
    cat("\n",
        "###########################","\n",
        "###Required Sample Size:###","\n",
        "###########################","\n",
        "Sample Size for t(n-1)-test:    ",N2,"\n","\n"
    )
  )
} else return("warning: test not properly specified")
}
