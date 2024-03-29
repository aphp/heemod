% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{calibrate_model}
\alias{calibrate_model}
\title{Calibrate Model Parameters}
\usage{
calibrate_model(
  x,
  parameter_names,
  fn_values,
  target_values,
  initial_values = NULL,
  method = c("Nelder-Mead", "BFGS", "L-BFGS-B"),
  ...
)
}
\arguments{
\item{x}{Result from \code{\link[=run_model]{run_model()}} or \code{\link[=update]{update()}}.}

\item{parameter_names}{Names of the parameters to
calibrate.}

\item{fn_values}{Function applied to the model that
returns the values of interest as a numeric vector.}

\item{target_values}{Values to match, same length as the
output from \code{fn_values}.}

\item{initial_values}{Optional starting values. See
details.}

\item{method}{Optimisation method (\code{Nelder-Mead},
\code{BFGS}, or \code{L-BFGS-B}).}

\item{...}{Optional arguments passed to
\code{\link[optimx:optimx]{optimx::optimx()}}.}
}
\value{
A data frame in which each row has the calibrated
values of parameters given in \code{parameter_names}, for
the corresponding row of \code{initial_values}, along with
the convergence code for each run.
}
\description{
Search for the appropriate value of unknown parameters to
obtain specific model results.
}
\details{
Parameters not being optimized are unchanged from the
values in the model run. If \code{initial_values} is \code{NULL},
the initial parameter values will also be taken from the
model run.

\code{initial_values} can be a vector or a table. In the
second case each row corresponds to a set of initial
parameter values: the calibration will be run once per
set.

Passing in multiple initial values allows (among other
things) the user to check whether the calibration gets
the same results from different starting points.

Multi-dimensional problems are optimized with
\code{\link[optimx:optimx]{optimx::optimx()}}, 1-dimensional problems with
\code{\link[stats:optimize]{stats::optimise()}} (except when a \code{method} is given).
\code{convcode} is always \code{NA} with \code{\link[stats:optimize]{stats::optimise()}}.

Running \code{\link[=calibrate_model]{calibrate_model()}} does not change the model
parameters; the user must create a new model and run it
if desired.

See also \code{vignette("k-calibration")}.
}
\examples{

param <- define_parameters(p = 0.8)

mat <- define_transition(
  p, C,
  0, 1
)
mod <- define_strategy(
  transition = mat,
  A = define_state(cost=10, effect = 0.5), 
  B = define_state(cost = 5, effect = 0.8)
)

res_mod <- run_model(
  mod = mod,
  parameters = param,
  init = c(1000L, 0L),
  cycles = 10,
  cost = cost,
  effect = effect,
  method = "end"
)

f <- function(x) {
  dplyr::filter(
    get_counts(x),
    state_names == "A" & model_time == 10
  )$count
}
f(res_mod)

#'\dontrun{
#'calibrate_model(
#'  res_mod,
#'  parameter_names = "p",
#'  fn_values = f,
#'  target_values = 130,
#'  initial_values = data.frame(p = c(0.5, 0.9)),
#'  lower = 0, upper = 1
#')
#'}
}
