#' Treatment effect estimation for AHMM with 2-level clustering (with covariates)
#'
#' This function conducts an estimation procedure of treatment effect
#' in an Additive Hazards Mixed Model with two-level clustering
#' and no additional covariates. The underlying methods are described in
#' (reference our paper).
#'
#' Lots of descriptive details could go here. This line just holds the place.
#'
#' @param formula Model specification in the traditional 'Y ~ X' sense,
#' where Y is in Surv format (function in 'survival' package).
#' @param cluster Variable in the data that contains information about group
#' (cluster) assignment of each individual observation
#' @param data Input data source
#' @param subset If only a subset of the data shall be used
#' @param method Test that should be used in the estimation procedure. Input
#' "z" if z-test or "t" if t-test should be used respectively.
#' @param bcv Input specifying which variance correction should be used in
#' the estimation procedure. "CZ" is the uncorrected variance as in Cai and
#' Zeng, "KC", "MD", "FG", and "MBN" are variance corrections according to
#' Kauermann and Carroll, Mancl and DeRouen, Fay and Graubard, Morel et al.
#' respectively.
#' @param na.action Specify how missing values should be handled.
#' @return This function treatment effect estimate of the underlying AHMM.
#' @keywords methods
#' @export
#' @examples
#' AHMM(Surv(U,Delta)~Z, cluster = ID, data = Test_Data)
#'
#'
#'

AHMM <- function(formula, cluster, data, subset, method = c("z", "t"), bcv = c("CZ", "KC", "MD", "FG", "MBN"), na.action)
{
  Call <- match.call()

  if (missing(formula))
    stop("a formula argument is required")
  ss <- "cluster"
  if (is.list(formula))
    Terms <- if (missing(data))
      terms(formula[[1]], specials = ss)
  else terms(formula[[1]], specials = ss, data = data)
  else Terms <- if (missing(data))
    terms(formula, specials = ss)
  else terms(formula, specials = ss, data = data)
  tcl <- attr(Terms, "specials")$cluster
  if (length(tcl) > 1)
    stop("a formula cannot have multiple cluster terms")
  if (length(tcl) > 0) {
    factors <- attr(Terms, "factors")
    if (any(factors[tcl, ] > 1))
      stop("cluster() cannot be in an interaction")
    if (attr(Terms, "response") == 0)
      stop("formula must have a Surv response")
    if (is.null(Call$cluster))
      Call$cluster <- attr(Terms, "variables")[[1 +
                                                  tcl]][[2]]
    else warning("cluster appears both in a formula and as an argument, formula term ignored")
    Terms <- drop.special(Terms, tcl)
    formula <- Call$formula <- formula(Terms)
  }
  indx <- match(c("formula", "data", "subset", "na.action", "cluster"),
                names(Call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")
  tform <- Call[c(1, indx)]
  tform[[1L]] <- quote(stats::model.frame)
  if (is.list(formula)) {
    multiform <- TRUE
    dformula <- formula[[1]]

    tlab <- unlist(lapply(covlist$rhs, function(x) attr(terms.formula(x$formula),
                                                        "term.labels")))
    tlab <- c(attr(terms.formula(dformula), "term.labels"),
              tlab)
    newform <- reformulate(tlab, dformula[[2]])
    environment(newform) <- environment(dformula)
    formula <- newform
    tform$na.action <- na.pass
  }
  else {
    multiform <- FALSE
    covlist <- NULL
    dformula <- formula
  }

  tform$formula <- if (missing(data))
    terms(formula)
  else terms(formula, data = data)

  mf <- eval(tform, parent.frame())
  Terms <- terms(mf)
  n <- nrow(mf)
  Y <- model.response(mf)
  isSurv2 <- inherits(Y, "Surv2")
  if (isSurv2) {
    if (length(attr(mf, "na.action"))) {
      tform$na.action <- na.pass
      mf <- eval.parent(tform)
    }
    if (!is.null(attr(Terms, "specials")$cluster))
      stop("cluster() cannot appear in the model statement")
    new <- surv2data(mf)
    mf <- new$mf
    Y <- new$y
    n <- nrow(mf)
  }
  else {
    if (!is.Surv(Y))
      stop("Response must be a survival object")
  }
  if (n == 0)
    stop("No (non-missing) observations")
  type <- attr(Y, "type")
  multi <- FALSE
  if (type == "mright" || type == "mcounting")
    multi <- TRUE
  else if (type != "right" && type != "counting")
    stop(paste("Our model doesn't support \"", type,
               "\" survival data", sep = ""))
  data.n <- nrow(Y)
  if (!multi && multiform)
    stop("formula is a list but the response is not multi-state")

  if (length(attr(Terms, "variables")) > 2) {
    ytemp <- rpsftm:::terms.inner(formula[1:2])
    suppressWarnings(z <- as.numeric(ytemp))
    ytemp <- ytemp[is.na(z)]
    xtemp <- rpsftm:::terms.inner(formula[-2])
    if (any(!is.na(match(xtemp, ytemp))))
      warning("a variable appears on both the left and right sides of the formula")
  }

  xlevels <- .getXlevels(Terms, mf)
  cluster <- model.extract(mf, "cluster")
  has.cluster <- !(missing(cluster) || length(cluster) == 0)

  if (!has.cluster) {
    stop("cluster variable not specified")
    ncluster <- 0
    clname <- NULL
  }
  else {
    if (is.factor(cluster)) {
      clname <- levels(cluster)
      cluster <- as.integer(cluster)
    }
    else {
      clname <- sort(unique(cluster))
      cluster <- match(cluster, clname)
    }
    ncluster <- length(clname)
  }


  #  if (length(dropterms)) {
  #    Terms2 <- Terms[-dropterms]
  #    X <- model.matrix(Terms2, mf)
  #    temp <- attr(X, "assign")
  #    shift <- sort(dropterms)
  #    for (i in seq(along.with = shift)) temp <- temp + 1 *
  #      (shift[i] <= temp)
  #    attr(X, "assign") <- temp
  #  }
  #  else
  X <- model.matrix(Terms, mf)
  Xatt <- attributes(X)
  #  if (hasinteractions)
  #    adrop <- c(0, untangle.specials(Terms, "strata")$terms)
  #  else
  adrop <- 0
  xdrop <- Xatt$assign %in% adrop
  X <- X[, !xdrop, drop = FALSE]
  attr(X, "assign") <- Xatt$assign[!xdrop]
  offset <- model.offset(mf)

  if (sum(Y[, ncol(Y)]) == 0) {
    ncoef <- ncol(X)
    ctemp <- rep(NA, ncoef)
    names(ctemp) <- colnames(X)
    concordance = c(concordant = 0, discordant = 0, tied.x = 0,
                    tied.y = 0, tied.xy = 0, concordance = NA, std = NA,
                    timefix = FALSE)
    rval <- list(coefficients = ctemp, var = matrix(0, ncoef,
                                                    ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset,
                 residuals = rep(0, data.n), means = colMeans(X),
                 method = method, n = data.n, nevent = 0, terms = Terms,
                 assign = assign, concordance = concordance, wald.test = 0,
                 y = Y, call = Call)
    class(rval) <- "coxph"
    return(rval)
  }

  if (!all(is.finite(X)))
    stop("data contains an infinite predictor")

  # Just renaming coxph objects into our notation, so we can use intact chunk of our code
  # for estimation while maintaining the coxph objects intact for later use
  ID<-clname
  U<-Y[,1]
  Z<-X
  #Z<-X[,1]
  #Xij<-X[,-1]
  Delta<-Y[,2]

  # General set up
  UID <- sort(unique(ID))
  n <- length(UID)
  nu <- length(U)
  ndelta <- dim(as.matrix(Z))[2]

  # Sort observations by time
  b <- order(U)
  U <- sort(U)
  Z <- Z[b,]
  ID <- ID[b]
  Delta <- Delta[b]

  # Increment of time
  dU <- U-c(0,U[1:(length(U)-1)])
  IndUU <- (t(pracma::repmat(U,length(U),1))>=pracma::repmat(U,length(U),1))
  IDind <- pracma::zeros(length(UID), length(U))
  for (i in 1:length(UID)){
    IDind[i, ID==UID[i]] <- 1
  }
  # the rate of counting process of event, or dN(t)
  NN <- diag(Delta)

  # Sum of U(t) at each t
  temp0 <- t(IndUU)%*%pracma::repmat(1,nu,1)
  # bar{Z}(t) at each t (row) and for each dimension (column)
  barZ <- (t(IndUU)%*%Z)/pracma::repmat(temp0,1,ndelta)
  nom <- pracma::zeros(ndelta,1)
  denom <- pracma::zeros(ndelta,ndelta)

  # Least-squares estimation of delta
  for (k in 1:ndelta){
    tempk <- pracma::repmat(as.matrix(Z[,k]),1,nu)-pracma::repmat(barZ[,k],nu,1)
    nom[k] <- sum(IndUU*tempk*NN)
    for (s in 1:ndelta){
      temps <- pracma::repmat(as.matrix(Z[,s]),1,nu)-pracma::repmat(barZ[,s],nu,1)
      denom[k,s] <- sum((IndUU*tempk*temps)%*%dU)
    }
  }
  delta_est <- pracma::inv(denom)%*%nom
  #Printing parameter estimate
  return(
    cat("\n",
        "###########################","\n",
        "###Estimated Parameters:###","\n",
        "###########################","\n",
        "Treatment effect:",as.numeric(delta_est[1,1]),"\n",
        "Additional covariates:",as.numeric(delta_est[-1,1]),"\n","\n"
    )
  )

  #  if (method == "breslow" || method == "efron") {
  #    if (grepl("right", type))
  #      fit <- coxph.fit(X, Y, istrat, offset, init,
  #                       control, weights = weights, method = method,
  #                       rname, nocenter = nocenter)
  #    else fit <- agreg.fit(X, Y, istrat, offset, init,
  #                          control, weights = weights, method = method,
  #                          rname, nocenter = nocenter)
  #  }
  #  else if (method == "exact") {
  #    if (type == "right")
  #      fit <- coxexact.fit(X, Y, istrat, offset, init,
  #                          control, weights = weights, method = method,
  #                          rname, nocenter = nocenter)
  #    else fit <- agexact.fit(X, Y, istrat, offset, init,
  #                            control, weights = weights, method = method,
  #                            rname, nocenter = nocenter)
  #  }
  #  else stop(paste("Unknown method", method))

  # Transfer of the estimates into coxph objects for further processing
  # fit$coefficients<-delta_est


  #  if (is.character(fit)) {
  #    fit <- list(fail = fit)
  #    class(fit) <- "coxph"
  #  }
  #  else {
  #    if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
  #      vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
  #      msg <- paste("X matrix deemed to be singular; variable",
  #                   paste(vars, collapse = " "))
  #      if (!singular.ok)
  #        stop(msg)
  #    }
  #    fit$n <- data.n
  #    fit$nevent <- sum(Y[, ncol(Y)])
  #    fit$terms <- Terms
  #    fit$assign <- assign
  #    class(fit) <- fit$class
  #    fit$class <- NULL
  #    if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
  #      fit$naive.var <- fit$var
  #      fit2 <- c(fit, list(x = X, y = Y, weights = weights))
  #      if (length(istrat))
  #        fit2$strata <- istrat
  #      if (length(cluster)) {
  #        temp <- residuals.coxph(fit2, type = "dfbeta",
  #                                collapse = cluster, weighted = TRUE)
  #        if (is.null(init))
  #          fit2$linear.predictors <- 0 * fit$linear.predictors
  #        else fit2$linear.predictors <- c(X %*% init)
  #        temp0 <- residuals.coxph(fit2, type = "score",
  #                                 collapse = cluster, weighted = TRUE)
  #      }
  #      else {
  #        temp <- residuals.coxph(fit2, type = "dfbeta",
  #                                weighted = TRUE)
  #        fit2$linear.predictors <- 0 * fit$linear.predictors
  #        temp0 <- residuals.coxph(fit2, type = "score",
  #                                 weighted = TRUE)
  #      }
  #      fit$var <- t(temp) %*% temp
  #      u <- apply(as.matrix(temp0), 2, sum)
  #      fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u,
  #                                control$toler.chol)$test
  #    }
  #    if (length(fit$coefficients) && is.null(fit$wald.test)) {
  #      nabeta <- !is.na(fit$coefficients)
  #      if (is.null(init))
  #        temp <- fit$coefficients[nabeta]
  #      else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
  #      fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta],
  #                                   temp, control$toler.chol)$test
  #    }
  #    if (length(cluster))
  #      temp <- concordancefit(Y, fit$linear.predictors,
  #                             istrat, weights, cluster = cluster, reverse = TRUE,
  #                             timefix = FALSE)
  #    else temp <- concordancefit(Y, fit$linear.predictors,
  #                                istrat, weights, reverse = TRUE, timefix = FALSE)
  #    if (is.matrix(temp$count))
  #      fit$concordance <- c(colSums(temp$count), concordance = temp$concordance,
  #                           std = sqrt(temp$var))
  #    else fit$concordance <- c(temp$count, concordance = temp$concordance,
  #                              std = sqrt(temp$var))
  #    na.action <- attr(mf, "na.action")
  #    if (length(na.action))
  #      fit$na.action <- na.action


  #    if (y)
  #      fit$y <- Y
  #    fit$timefix <- control$timefix
  #  }

  #    names(fit$coefficients) <- paste0(names(fit$coefficients),
  #                                      suffix)
  #    if (x)
  #      fit$strata <- istrat
  #    class(fit) <- c("coxphms", class(fit))
  #  }
  #  names(fit$means) <- names(fit$coefficients)
  #  fit$formula <- formula(Terms)
  #  if (length(xlevels) > 0)
  #    fit$xlevels <- xlevels

  #  fit$call <- Call
  #  fit
}
