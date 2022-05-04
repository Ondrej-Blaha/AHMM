#' VIF and ICC calculation for AHMM with 2-level clustering
#'
#' This function performs a VIF and ICC calculation for an underlying test of
#' treatment effect under the Additive Hazards Mixed Model with two-level
#' clustering and randomization occurring at the higher (cluster) level.
#' The underlying methods are described in (2022 Blaha, Esserman, Li paper).
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
#' @return This function returns returns VIF and ICC of the underlying AHMM.
#' @keywords methods
#' @export
#' @examples
#' VIF(lambda0=4,delta=1.5,sigma2=0.5,tau=1,m_bar=10,CV=0.25)
#'
#'
#'
VIF<-function(lambda0,delta,sigma2,tau,m_bar,CV=0,pi=0.5,alpha=0.05,
              beta=0.2){
  # ICC
  rho = (pracma::dblquad(two,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(three,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(four,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(five,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi) +
         pracma::dblquad(six,0,tau,0,tau,dim=2,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi)) /
         pracma::integral(one,xmin=0,xmax=tau,tau=tau,delta=delta,lambda0=lambda0,sigma2=sigma2,pi=pi)

  # VIF
  VIF=(1+((1+CV^2)*m_bar-1)*rho)

  return(
    cat("\n",
        "##########################","\n",
        "###Estimated ICC & VIF:###","\n",
        "##########################","\n",
        "ICC =",rho,"\n",
        "VIF =",VIF,"\n","\n"
        )
  )
}
