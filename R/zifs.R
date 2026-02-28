#' Zero-Inflated Flory-Schulz Probability Mass Function
#'
#' Computes the probability mass function (PMF) of the
#' Zero-Inflated Flory-Schulz (ZIFS) distribution.
#'
#' @param x Non-negative integer value(s).
#' @param p Parameter satisfying 0 < p < 1.
#' @param beta Zero-inflation parameter satisfying 0 <= beta < 1.
#'
#' @return Probability mass at x.
#'
#' @examples
#' dzifs(1, 0.8, 0.5)
#' dzifs(0, 0.3, 0.75)
#'
#' @export
dzifs <- function(x, p, beta) {
  ifelse(
    x < 0,
    0,
    ifelse(
      x == 0,
      beta,
      (1 - beta) * x * (1-p)^(x-1) * p^2
    )
  )
}


#' Zero-Inflated Flory-Schulz Cumulative Distribution Function
#'
#' Computes the cumulative distribution function (CDF)
#' of the ZIFS distribution.
#'
#' @param x Numeric value.
#' @param p Parameter satisfying 0 < p < 1.
#' @param beta Zero-inflation parameter satisfying 0 <= beta < 1.
#'
#' @return Cumulative probability P(X <= x).
#'
#' @examples
#' pzifs(1, 0.8, 0.5)
#' pzifs(0.4, 0.3, 0.75)
#'
#' @export
pzifs <- function(x, p, beta) {
  if (x < 0) return(0)
  x <- floor(x)
  sum(dzifs(0:x, p, beta))
}


#' Zero-Inflated Flory-Schulz Survival Function
#'
#' Computes the survival function P(X > x) of the ZIFS distribution.
#'
#' @inheritParams pzifs
#'
#' @return Survival probability.
#'
#' @examples
#' szifs(1, 0.8, 0.5)
#'
#' @export
szifs <- function(x, p, beta) {
  1 - pzifs(x, p, beta)
}


#' Mean of the ZIFS Distribution
#'
#' Computes the theoretical mean of the ZIFS distribution.
#'
#' @param beta Zero-inflation parameter (0 <= beta < 1).
#' @param p Parameter (0 < p < 1).
#'
#' @return Mean value.
#'
#' @examples
#' zifs_mean(0.3, 0.4)
#'
#' @export
zifs_mean <- function(beta, p) {


  if (beta < 0 || beta >= 1) {
    stop("beta must satisfy 0 <= beta < 1")
  }

  if (p <= 0 || p >= 1) {
    stop("p must satisfy 0 < p < 1")
  }


  mean_value <- (1 - beta) * (2 - p) / p

  return(mean_value)
}


#' Variance of the ZIFS Distribution
#'
#' Computes the theoretical variance of the ZIFS distribution.
#'
#' @inheritParams zifs_mean
#'
#' @return Variance value.
#'
#' @examples
#' zifs_variance(0.3, 0.4)
#'
#' @export
zifs_variance <- function(beta, p) {


  if (beta < 0 || beta >= 1) {
    stop("beta must satisfy 0 <= beta < 1")
  }

  if (p <= 0 || p >= 1) {
    stop("p must satisfy 0 < p < 1")
  }


  variance_value <- (1 - beta) *
    ((beta * p^2 - (2 + 4 * beta) * p + 2 + 4 * beta) / p^2)

  return(variance_value)
}


#' Index of Dispersion of the ZIFS Distribution
#'
#' Computes the index of dispersion (variance/mean).
#'
#' @inheritParams zifs_mean
#'
#' @return Index of dispersion value.
#'
#' @examples
#' zifs_iod(0.3, 0.4)
#'
#' @export
zifs_iod <- function(beta, p) {


  if (beta < 0 || beta >= 1) {
    stop("beta must satisfy 0 <= beta < 1")
  }

  if (p <= 0 || p >= 1) {
    stop("p must satisfy 0 < p < 1")
  }


  iod_value <- (beta * (p - 2)^2 + 2 * (1 - p)) / (p * (2 - p))

  return(iod_value)
}


#' Quantile Function of the ZIFS Distribution
#'
#' Computes the quantile function of the ZIFS distribution
#' using the Lambert W function.
#'
#' @param r Probability value (0 <= r < 1).
#' @param p Parameter (0 < p < 1).
#' @param beta Zero-inflation parameter (0 <= beta < 1).
#'
#' @return Quantile corresponding to probability r.
#'
#' @examples
#' qzifs(0.7, 0.4, 0.3)
#'
#' @export
qzifs <- function(r, p, beta) {


  if (r < 0 || r >= 1) stop("r must satisfy 0 <= r < 1")
  if (p <= 0 || p >= 1) stop("p must satisfy 0 < p < 1")
  if (beta < 0 || beta >= 1) stop("beta must satisfy 0 <= beta < 1")


  if (r <= beta) {
    return(0)
  }


  log_term <- log(1 - p)

  A <- (1 - (r - beta)/(1 - beta)) *
    (log_term / p) *
    (1 - p)^(1/p)

  W1 <- lamW::lambertWm1(A)

  c_val <- (1 / log_term) * W1 - (1 / p)


  return(ceiling(c_val))
}


#' Mode of the ZIFS Distribution
#'
#' Computes the mode of the ZIFS distribution.
#'
#' @param p Parameter (0 < p < 1).
#' @param beta Zero-inflation parameter (0 <= beta < 1).
#'
#' @return Mode value. May return two values in bimodal case.
#'
#' @examples
#' modezifs(0.4, 0.3)
#'
#' @export
modezifs <- function(p, beta) {


  if (p <= 0 || p >= 1) stop("p must satisfy 0 < p < 1")
  if (beta < 0 || beta >= 1) stop("beta must satisfy 0 <= beta < 1")


  if (beta == 0) {
    return(floor(1 / p))
  }

  k_star <- floor(1 / p)
  max_term <- k_star * (1 - p)^(k_star - 1)


  if (p > 0.5) {
    if (beta > (1 - beta) * p^2 * max_term) {
      return(0)
    } else {
      return(k_star)
    }
  }


  if (p <= 0.5) {
    if (beta > (1 - beta) * p^2) {
      return(c(0, k_star))  # bimodal case
    } else {
      return(k_star)
    }
  }
}


#' Random Generation from the ZIFS Distribution
#'
#' Generates random observations from the ZIFS distribution.
#'
#' @param k Sample size.
#' @param p Parameter (0 < p < 1).
#' @param beta Zero-inflation parameter (0 <= beta < 1).
#'
#' @return Numeric vector of random values.
#' @importFrom stats runif
#' @export
rzifs <- function(k, p, beta) {
  u <- runif(k)
  y <- numeric(k)

  for (i in seq_len(k)) {
    y[i] <- qzifs(u[i], p, beta)
  }
  y
}


#' Maximum Likelihood Estimation for ZIFS
#'
#' Computes MLE estimators for p and beta.
#'
#' @param input_data Numeric vector of observations.
#'
#' @return Named vector with estimates.
#'
#' @export
MLE_ZIFS <- function(input_data){
  n = length(input_data)
  n_0 = sum(input_data == 0)
  m_1 = mean(input_data)
  beta_hat = n_0/n
  p_hat = 2/((n*m_1)/(n-n_0) + 1)
  print(paste("Estimated value of p through MLE:", p_hat))
  print(paste("Estimated value of Beta through MLE:", beta_hat))
}


#' Method of Moments Estimation for ZIFS
#'
#' Computes method-of-moments estimators for p and beta.
#'
#' @param input_data Numeric vector of observations.
#'
#' @return Named vector with estimates.
#'
#' @export
mom_ZIFS <- function(input_data){
  m_1 = mean(input_data)
  m_2 = mean(input_data^2)
  beta_hat = 1 - (
    m_1 * ((3*m_1 + m_2) - sqrt(3*m_1^2 + m_2^2))
  ) / (
    -m_1 + m_2 + sqrt(3*m_1^2 + m_2^2)
  )
  p_hat = ((3*m_1 + m_2) - sqrt(3*m_1^2 + m_2^2)) / (m_1 + m_2)
  print(paste("Estimated value of p through MOM:", p_hat))
  print(paste("Estimated value of Beta through MOM:", beta_hat))
}


#' Log-Likelihood of the ZIFS Distribution
#'
#' Computes the log-likelihood for observed data.
#'
#' @param n_0 Number of zeros in the sample.
#' @param x Vector of non-zero observations.
#' @param p Parameter (0 < p < 1).
#' @param beta Zero-inflation parameter (0 <= beta < 1).
#'
#' @return Log-likelihood value.
#'
#' @export
ll_zifs=function(n_0,x,p,beta){
  n_0*log(beta) + sum(log((1 - beta) * x * (1-p)^(x-1) * p^2))
}

#' AIC for ZIFS Model
#'
#' Computes Akaike Information Criterion.
#'
#' @inheritParams ll_zifs
#'
#' @return AIC value.
#'
#' @export
aic_zifs = function(n_0,x,p,beta){2*2 - 2*ll_zifs(n_0,x,p,beta)}

#' BIC for ZIFS Model
#'
#' Computes Bayesian Information Criterion.
#'
#' @param sample_size Total sample size.
#' @inheritParams ll_zifs
#'
#' @return BIC value.
#'
#' @export
bic_zifs = function(n_0,x,sample_size,p,beta){2*log(sample_size) - 2*ll_zifs(n_0,x,p,beta)}

