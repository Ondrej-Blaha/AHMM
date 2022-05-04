#' Permutation test calculation for AHMM with 2-level clustering
#'
#' This function performs a permutation test calculation for an underlying test of
#' treatment effect under the Additive Hazards Mixed Model with two-level
#' clustering and randomization occurring at the higher (cluster) level.
#' The underlying methods are described in (2022 Blaha, Esserman, Li paper).
#'
#'
#' @param ID Vector specifying in which group/cluster each individual
#' sample belongs to
#' @param Z_cl Vector of actual allocation into treatment arms
#' @param U Vector with observed event times
#' @param Delta Censoring vector
#' @param pmt User specified permutation matrix
#' @return This function returns returns randomization-based test p-value
#' of the underlying AHMM
#' @keywords methods
#' @export
#' @examples
#' PermTest(ID=Test_Data$ID, Z_cl=trt, U=Test_Data$U, Delta=Test_Data$Delta, pmt = pmt)
#'
#'
#'
PermTest=function(ID, Z_cl, U, Delta, pmt = NULL){

  Temp<-data.frame(ID,Z_cl)
  # Sort observations by time
  b <- order(U)
  U <- sort(U)
  ID <- ID[b]
  Delta <- Delta[b]

  # Increment of time
  dU <- U-c(0,U[1:(length(U)-1)])
  # Each row is an individual, each column is a specific time point
  IndUU <- (t(pracma::repmat(U,length(U),1))>=pracma::repmat(U,length(U),1))
  UID <- sort(unique(ID))
  IDind <- pracma::zeros(length(UID), length(U))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }

  # General set up
  n <- length(UID)
  nu <- length(U)
  # Extract trt from the data
  Ranks <- with(Temp, ave(Z_cl, ID, FUN = function(x) rank(x, ties.method="first")))
  trt<-Temp[Ranks == 1,"Z_cl"]
  # The rate of counting process of event, or dN(t)
  NN <- diag(Delta)
  # Sum of U(t) at each t
  temp0 <- t(IndUU)%*%pracma::repmat(1,nu,1)

  # Estimation of marginal baseline hazard
  dHU <- colSums(IndUU*NN)/c(temp0)
  HU <- cumsum(dHU)

  # Permutation distribution of score
  ######################################################################################
  dH <- dHU
  epsilon <- NN-pracma::repmat(t(dH),nu,1)
  score <- as.numeric(IDind%*%(as.matrix(IndUU*epsilon)%*%pracma::repmat(1,nu,1)))
  ######################################################################################


  # Check if pmt matrix has been specified and if not then generate
  if (is.null(pmt)==TRUE){
    if(choose(n, n/2)<=50000){
      # Enumeration of all possible schemes
      emr <- t(combn(n, n / 2))
      R <- dim(emr)[1]
      pmt <- matrix(0, R, n)
      for (r in 1:R){
        pmt[r, emr[r,]] <- 1
      }
    } else{
      S=50000
      # Simulation of permutations
      pmt <- matrix(0, S, n)
      for (s in 1:S){
        trt_perm <- sample(1:n, n / 2)
        pmt[s, trt_perm] <- 1
      }
      pmt <- unique(pmt) # Indicator matrix
      R <- dim(pmt)[1]
    }
  } else{
    #Check validity of the input matrix

  }


  obs = sum(trt*score)
  Pdist = c(pmt%*%score)

  # Printing permutation p-value
  return(
    cat("\n",
        "###############################","\n",
        "###Randomization-based test:###","\n",
        "###############################","\n",
        "p-val =",mean(abs(Pdist) >= abs(obs)),"\n","\n"
    )
  )
}
