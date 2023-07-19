#' Probability Density Functions for Probabilistic
#' Uncertainty Analysis
#' 
#' Define a distribution for PSA parameters.
#' 
#' These functions are not exported, but only used
#' in [define_psa()]. To specify a user-made 
#' function use [define_distribution()].
#' 
#' [use_distribution()] uses gaussian kernel 
#' smoothing with a bandwith parameter calculated 
#' by [stats::density()]. Values for unobserved
#' quantiles are calculated by linear
#' interpolation.
#' 
#' @name distributions
#' @param mean Distribution mean.
#' @param sd Distribution standard deviation.
#' @param ... Dirichlet distribution parameters.
#' @param prob Proportion.
#' @param size Size of sample used to estimate 
#'   proportion.
#' @param meanlog Mean on the log scale.
#' @param sdlog SD on the log scale.
#' @param mu Mean on the logit scale.
#' @param sigma SD on the logit scale.
#' @param shape1 for beta distribution
#' @param shape2 for beta distribution
#' @param lower lower bound of triangular 
#'   distribution.
#' @param upper upper bound of triangular 
#'   distribution.
#' @param peak peak of triangular distribution.
#' @param x A distribution function, see details.
#' @param distribution A numeric vector of 
#'   observations defining a distribution, usually
#'   the output from an MCMC fit.
#' @param smooth Use gaussian kernel smoothing?
#'   
#' @examples 
#' 
#' # a mixture of 2 gaussians
#' x <- c(rnorm(100), rnorm(100, 6))
#' plot(density(x))
#' 
#' use_distribution(x)
#' 
normal <- function(mean, sd) {
  function(x) stats::qnorm(p = x, mean = mean, sd = sd)
}

#' @rdname distributions
lognormal <- function(mean, sd, meanlog, sdlog) {
  if (missing(sdlog)) sdlog <- sqrt(log(1 + sd^2/mean^2))
  if (missing(meanlog)) meanlog <- log(mean) - sdlog^2/2
  function(x) stats::qlnorm(p = x, meanlog = meanlog, sdlog = sdlog)
}

#' @rdname distributions
gamma <- function(mean, sd) {
  shape <- mean^2/sd^2
  scale <- sd^2/mean
  function(x) stats::qgamma(p = x, shape = shape, scale = scale)
}

#' @rdname distributions
binomial <- function(prob, size) {
  function(x) stats::qbinom(p = x, size = size, prob = prob) / size
}

#' @rdname distributions
multinomial <- function(...) {
  list_param <- list(...)
  
  structure(
    lapply(list_param, function(x) r_multinomial(x)),
    class = "multinom_param"
  )
}
r_multinomial <- function(n) {
  function(x) stats::qgamma(x, shape = n, scale = 1)
}

#' @rdname distributions
logitnormal <- function(mu, sigma) {
  if (! requireNamespace("logitnorm")) {
    stop("'logitnorm' package required for logitnormal distributions.")
  }
  
  function(x) logitnorm::qlogitnorm(p = x, mu = mu, sigma = sigma)
}

#' @rdname distributions
beta <- function(shape1, shape2){
  function(x){stats::qbeta(p = x, shape1 = shape1, shape2 = shape2)}
}

#' @rdname distributions
triangle <- function(lower, upper, peak = (lower + upper)/2) {
  if (! requireNamespace("triangle")) {
    stop("'triangle' package required for logitnormal distributions.")
  }
  stopifnot(peak >= lower,
            upper >= peak,
            upper > lower
            )
  function(x) triangle::qtriangle(p = x, a = lower, b = upper, c = peak)
}


#' @rdname distributions
poisson <- function(mean) {
  function(x) stats::qpois(p = x, lambda = mean)
}

#' @rdname distributions
#' @export
use_distribution <- function(distribution, smooth = TRUE) {
  distribution <- sort(distribution)
  
    function(x) {
      if (smooth) {
        noise <- stats::rnorm(
          n = length(x),
          mean = 0,
          sd = stats::density(distribution)$bw
        )
      } else {
        noise <- 0
      }
      
      stats::approxfun(
        x = seq(0, 1, length = length(distribution)),
        y = distribution)(x) + noise
    }
}
