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
#' smoothing with a bandwidth parameter calculated 
#' by [stats::density()]. Values for unobserved
#' quantiles are calculated by linear
#' interpolation.
#' 
#' [define_distribution()] takes as argument a 
#' function with a single argument, `x`, 
#' corresponding to a vector of quantiles. It 
#' returns the distribution values for the given 
#' quantiles. See examples.
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
#' define_distribution(
#'   function(x) stats::qexp(p = x, rate = 0.5)
#' )
#' 
#' # a mixture of 2 gaussians
#' x <- c(rnorm(100), rnorm(100, 6))
#' plot(density(x))
#' 
#' use_distribution(x)
normal <- function(mean, sd) {
  list(r_normal(mean, sd))
}
r_normal <- function(mean, sd) {
  function(x) stats::qnorm(p = x, mean = mean, sd = sd)
}

#' @rdname distributions
lognormal <- function(mean, sd, meanlog, sdlog) {
  if (missing(sdlog)) sdlog <- sqrt(log(1 + sd^2/mean^2))
  if (missing(meanlog)) meanlog <- log(mean) - sdlog^2/2
  
  list(r_lognormal(meanlog, sdlog))
}
r_lognormal <- function(meanlog, sdlog) {
  function(x) stats::qlnorm(p = x, meanlog = meanlog, sdlog = sdlog)
}

#' @rdname distributions
gamma <- function(mean, sd) {
  list(r_gamma(mean^2/sd^2, sd^2/mean))
}

r_gamma <- function(shape, scale) {
  function(x) stats::qgamma(p = x, shape = shape, scale = scale)
}

#' @rdname distributions
binomial <- function(prob, size) {
  list(r_binomial(prob, size))
}

r_binomial <- function(prob, size) {
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
  list(r_logitnormal(mu, sigma))
}

r_logitnormal <- function(mu, sigma) {
  function(x) logitnorm::qlogitnorm(p = x, mu = mu, sigma = sigma)
}

#' @rdname distributions
beta <- function(shape1, shape2){
  list(r_beta(shape1, shape2))
}

r_beta <- function(shape1, shape2){
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
  list(r_triangle(lower, upper, peak))
}
r_triangle <- function(lower, upper, peak) {
  function(x) triangle::qtriangle(p = x, a = lower, b = upper, c = peak)
}


#' @rdname distributions
poisson <- function(mean) {
  list(r_poisson(mean))
}
r_poisson <- function(mean) {
  function(x) stats::qpois(p = x, lambda = mean)
}

#' @rdname distributions
#' @export
define_distribution <- function(x) {
  list(x)
}

#' @rdname distributions
beta <- function(shape1, shape2){
  list(r_beta(shape1, shape2))
}
r_beta <- function(shape1, shape2){
  function(x){stats::qbeta(p = x, shape1 = shape1, shape2 = shape2)}
}

#' @rdname distributions
triangle <- function(lower, upper, peak = (lower + upper)/2) {
  if (! requireNamespace("triangle")) {
    stop("'triangle' package required for triangle distributions.")
  }
  stopifnot(peak >= lower,
            upper >= peak,
            upper > lower
  )
  list(r_triangle(lower, upper, peak))
}
r_triangle <- function(lower, upper, peak) {
  function(x) triangle::qtriangle(p = x, a = lower, b = upper, c = peak)
}

#' @rdname distributions
#' @export
use_distribution <- function(distribution, smooth = TRUE) {
  distribution <- sort(distribution)
  
  define_distribution(
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
  )
}

#' Resample survival distribution
#' 
#' 
#' @param x a `surv_object`
#' @param n the number of observations to generate if dist is specified or x is a `surv_dist`
#' object
#' 
#' The lower n is, the higher is the variability
#' @export
resample_surv <- function(n){
  if (missing(n)) return(
    structure(expr(resample_surv_boot()),
              class = c("surv_psa"))
  ) 
  structure(expr(resample_surv_dist(!!n)),
            class = c("surv_psa"))
}


#' @rdname resample_surv
#' @keywords internal
resample_surv_boot <- function(x){
  structure(list(r_boot_survfit(x)),
            class = c("surv_psa"))
}


#' @rdname resample_surv
#' @keywords internal
resample_surv_dist <- function(x, n){
  if (inherits(x, "surv_fit")){
    cli::cli_warn("{.fn resample_surv} should not contain the {.arg n} argument \
                  for {.cls surv_fit} object",
                  .frequency = "regularly",
                  .frequency_id = "resample_surv_dist")
    return(resample_surv_boot(x))
  }
  if (! requireNamespace("flexsurv")) {
    stop("'flexsurv' package required.")
  }
  pf <- get(paste0("r", x$distribution),
            envir = asNamespace("flexsurv"))
  
  args <- x[- match("distribution", names(x))]
  ret <- do.call(pf, c(n, args))
  structure(list(r_resample_surv_dist(ret, x$distribution, args)),
            class = c("surv_psa"))
}

#' @rdname resample_surv
#' @keywords internal
r_resample_surv_dist <- function(distribution, type, args){
  if (!missing(args)){
    y <- ecdf(distribution)(distribution)
    df <- list(x = distribution, y = y)
    args2 <- setNames(syms(names(args)), names(args))
    rhs <- rlang::call2(paste0("p", type), quote(x), !!!args2)
    formula <- rlang::new_formula(quote(y), rhs, env = asNamespace("flexsurv"))
    fit <- try(nls(
      formula, 
      data = df, start = args
    ), silent = TRUE)
    if (inherits(fit, "try-error")){
      fit <- nls(
        formula, 
        data = df, start = args,
        algorithm = "port",
        lower = 1E-1
      )
    }
    return(do.call(define_surv_dist, c(type, as.list(coef(fit)))))
  }
}

#' @rdname resample_surv
#' @keywords internal
r_boot_survfit <- function(x){
  init_surv_object <- x
  data <- rlang::call_args(x) %>% 
    `[[`("data") 
  e_data <- data %>% 
    eval_tidy()
  if (is.null(attr(x, "strata"))){
    new_data <- e_data[sample.int(nrow(e_data), 
                                replace = TRUE),]
  } else {
    strata <- 
      names(attr(x, "strata")) %>% strsplit(., ", ") %>% 
      unlist(use.names = F) %>% 
      gsub("=.*", "", .) %>%
      unique()
    new_data <-  e_data %>%
      split(.[[strata]]) %>% 
      lapply(function(x){
        x[sample.int(nrow(x), replace = TRUE),]
      }) %>% 
      bind_rows()
  }
 # assign(deparse(data), new_data, envir = getOption("heemod.env"))
  
   new_env <- rlang::env()
   assign(deparse(data), new_data, envir = new_env)
  # res <- rlang::call_modify(x, data = quote(new_data))
  # res <- new_quosure(res, new_env)
  
  res <- new_quosure(data, env = new_env)
  
  # if (is.list(init_surv_object) && "dist" %in% names(init_surv_object)){
  #   structure(c(list(dist = res),
  #               init_surv_object[setdiff(names(init_surv_object), "dist")]),
  #             class = class(init_surv_object))
  # } else res
  res
}
