vswitch <- function(x, ...)
  UseMethod("vswitch")

vswitch.factor <- function(x, ...) {
  x <- levels(x)[x]
  
  vswitch(x, ...)
}

vswitch.character <- function(x, ...) {
  listRes <- list(...)
  nRes <- names(listRes)
  if (is.null(nRes))
    stop("'...' should be named when 'x' is a character vector.")
  if (any(!is.na(x) & ! x %in% nRes))
    stop("Some values of 'x' do not correspond to any name in '...'.")
  x <- match(x, nRes)
  names(listRes) <- NULL
  
  do.call(vswitch, c(list(x = x), listRes))
}

vswitch.default <- function(x, ...) {
  listRes <- lapply(list(...),
                    function(x) rep(x, length.out = length(x)))
  if (! is.null(names(listRes)))
    warning("Named '...' with non character/factor 'x'. Names will be ignored.")
  if (any(!is.na(x) & (x > length(listRes) | x < 1)))
    stop("Some values of 'x' are out of range.")
  
  tabRes <- as.data.frame(listRes)
  iRows <- seq(length.out = nrow(tabRes))
  
  as.matrix(tabRes)[cbind(iRows, x)]
}

#' Dispatch Values According to Strategy
#' 
#' Returns different values depending on the strategy.
#' 
#' @param .strategy Optional strategy name. If not specified
#'   it is implicitely added.
#' @param ... Values of the parameter named depending on the
#'   strategy.
#'   
#' @return A vector of values.
#' @export
#' 
#' @examples
#' 
#' define_parameters(
#'   val = 456,
#'   x = dispatch_strategy(
#'     strat_1 = 1234,
#'     strat_2 = 9876,
#'     strat_3 = val * 2 + model_time
#'   )
#' )
dispatch_strategy <- function(.strategy, ...) {
  .dots <- list(...)
  if (is.null(names(.dots)) || any(is.na(names(.dots)))) {
    stop("All arguments to 'dispatch_strategy()' must be named.")
  }
  if (! is.character(.strategy)) {
    stop("'.strategy' must be a character vector.")
  }
  if (any(is.na(.strategy))) {
    stop("Missing data in '.strategy'.")
  }
  vswitch(.strategy, ...)
}

#' Hack to Automate Use of Strategy Name
#' 
#' This function is a hack to automate the definition of the
#' argument `.strategy` in
#' [dispatch_strategy()].
#' 
#' The hack consists in replacing calls to 
#' `dispatch_strategy(...)` by 
#' `dispatch_strategy(.strategy = strategy, ...)` if
#' `.strategy_name` is not already defined.
#' 
#' @param .dots An `expressions` object.
#'   
#' @return A modified `expressions` object.
#'   
#' @keywords internal
dispatch_strategy_hack <- function(.dots) {
  f <- function (x) {
    if (is.atomic(x) || is.name(x)) {
      x
    } else if (is.call(x)) {
      if (dispatch_strategy_check(x[[1]])) {
        x <- call_standardise(x)
        if (is.null(x$.strategy)) {
          x$.strategy <- substitute(strategy)
        }
      }
      as.call(lapply(x, f))
    } else if (is.pairlist(x)) {
      as.pairlist(lapply(x, f))
    } else {
      stop(sprintf(
        "Don't know how to handle type %s.",
        typeof(x)))
    }
  }
  
  do.call(
    structure,
    c(list(
      .Data = lapply(
        .dots, f
      )),
      attributes(.dots)
    )
  )
}

# Ensure only heemod version of dispatch_strategy gets used
dispatch_strategy_check <- function(x) {
  if (identical(x, quote(dispatch_strategy))) {
      TRUE
  } else {
    FALSE
  }
}
