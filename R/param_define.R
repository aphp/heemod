#' Define Markov Model Parameters
#' 
#' Define parameters called to compute the transition matrix
#' or state values for a Markov model. Parameters can be 
#' time dependent by using the `model_time` 
#' parameter.
#' 
#' Parameters are defined sequentially, parameters defined 
#' earlier can be called in later expressions.
#' 
#' Vector length should not be explicitly set, but should 
#' instead be stated relatively to `model_time` 
#' (whose length depends on the number of simulation 
#' cycles). Alternatively, `dplyr` functions such as 
#' [dplyr::n()] can be used.
#' 
#' This function relies heavily on the `dplyr` package.
#' Parameter definitions should thus mimic the use of 
#' functions such as [dplyr::mutate()].
#' 
#' Variable names are searched first in the parameter 
#' definition (only parameters defined earlier are visible) 
#' then in the environment where `define_parameters` 
#' was called.
#' 
#' For the `modify` function, existing parameters are 
#' modified, but no new parameter can be added. Parameter 
#' order matters since only parameters defined earlier can 
#' be referenced in later expressions.
#' 
#' @param ... Name-value pairs of expressions defining 
#'   parameters.
#' @param .OBJECT An object of class 
#'   `uneval_parameters`.
#' @param .dots Used to work around non-standard evaluation.
#'   
#' @return An object of class `uneval_parameters` 
#'   (actually a named list of quosures).
#' @export
#' 
#' @importFrom dplyr n row_number
#' @example inst/examples/example_define_parameters.R
#'   
define_parameters <- function(...) {
  .dots <- exprs(...)
  lapply(.dots, function(x){
    if (is.symbol(x)) return(TRUE)
    if(identical(x[[1]], quote(`<-`))) {
      cli::cli_abort("use of assignment symbol (`<-`) for variable {.var {x[[2]]}}. 
                     Replace with `=` symbol.")
    }
  })
  .dots <- as_quosures(.dots, getOption("heemod.env"))
  deprecated_x_cycle(.dots)
  define_parameters_(.dots)
}

#' @rdname define_parameters
#' @export
define_parameters_ <- function(.dots) {
  
  if (length(.dots)){
    check_names(names(.dots))
  }
  structure(.dots,
            class = c("uneval_parameters", class(.dots)))
}

#' Return parameters names
#' 
#' Extract parameters names.
#' 
#' @param x An object with parameters.
#'   
#' @return A character vector of parameter names.
#'   
#' @keywords internal
get_parameter_names <- function(x) {
  UseMethod("get_parameter_names")
}
#' @export
get_parameter_names.updated_model <- function(x) {
  get_parameter_names(get_model(x))
}

#' @export
get_parameter_names.uneval_parameters <- function(x) {
  names(x)[! names(x) %in% c("markov_cycle", "strategy",
                             "model_time")]
}
#' @export
get_parameter_names.eval_parameters <- function(x) {
  get_parameter_names.uneval_parameters(x)
}
#' @export
get_parameter_names.run_model <- function(x) {
  get_parameter_names(get_parameters(x))
}

#' Modify Object
#' 
#' This generic function allows the modification of various 
#' objects such as parameters, transitions matrix or states.
#' 
#' More details are available on the respective help page of
#' each object definition.
#' 
#' @param .OBJECT Various objects.
#' @param ... Modifications.
#'   
#' @return Same class as `x`.
#' @export
#' 
modify <- function(.OBJECT, ...) {
  UseMethod("modify")
}

modify_ <- function(.OBJECT, .dots) {
  UseMethod("modify_")
}

#' @export
#' @rdname define_parameters
modify.uneval_parameters <- function(.OBJECT, ...) {
  .dots <- quos(...)
  
  modify_(.OBJECT = .OBJECT, .dots = .dots)
}

#' @export
modify_.uneval_parameters <- function(.OBJECT, .dots) {
  if (length(.dots)) {
    check_names(names(.dots))
    utils::modifyList(.OBJECT, .dots)
  } else {
    .OBJECT
  }
}

#' Define Inflow for a BIA
#' 
#' @param ... Name-value pairs of expressions defining
#'   inflow counts.
#' @param .dots Used to work around non-standard evaluation.
#'   
#' @return An object similar to the return value of
#'   [define_parameters()].
#' @export
define_inflow <- function(...) {
  .dots <- quos(...)
  define_inflow_(.dots)
}

#' @export
#' @rdname define_inflow
define_inflow_ <- function(.dots) {
  structure(
    .dots,
    class = c("uneval_inflow", class(.dots)))
}

#' Define Initial Counts
#' 
#' @param ... Name-value pairs of expressions defining
#'   initial counts.
#' @param .dots Used to work around non-standard evaluation.
#'   
#' @return An object similar to the return value of
#'   [define_parameters()].
#' @export
define_init <- function(...) {
  .dots <- quos(...)
  define_init_(.dots)
}

#' @export
#' @rdname define_init
define_init_ <- function(.dots) {
  structure(
    .dots,
    class = c("uneval_init", class(.dots)))
}

#' Define Starting State Values
#' 
#' This function is meant to be used inside [define_strategy()] and 
#' [define_state()]. 
#' 
#' @param ... Name-value pairs of expressions defining
#'   starting values. The names must correspond to an existing state value.
#' @param .dots Used to work around non-standard evaluation.
#'
#' @details The behaviour is different following the function using [define_starting_values()]
#' as an argument.
#' \itemize{
#' \item When used inside [define_strategy()], the state values are modified for the 
#' first cycle in each state
#' \item When used inside [define_state()], the state values are modified for counts
#' entering the state
#' }
#' @return An object similar to the return value of
#'   [define_parameters()].
#' @export
define_starting_values <- function(...) {
  .dots <- quos(...)
  define_starting_values_(.dots)
}

#' @export
#' @rdname define_starting_values
define_starting_values_ <- function(.dots) {
  structure(
    .dots,
    class = c("uneval_starting_values", class(.dots)))
}

to_check <- "'define_init()', 'define_inflow()' or 'define_starting_values()'"

check_init <- function(x, ref) {
  UseMethod("check_init")
}

#' @export
check_init.quosures <- function(x, ref) {
  original_class <- class(x)
  
  if (length(x)) {
    if (is.null(names(x)) || all(names(x) == "")) {
      stop(to_check, " values must be named.")
    }
    
    if (! all(names(x) %in% ref)) {
      stop("Some ", to_check, " names are incorrect.")
    }
    
    if (any(duplicated(names(x)))) {
      stop("Duplicated names in ", to_check, ".")
    }
  }
  
  # res <- stats::setNames(
  #   object = as_quosures(
  #     lapply(ref, function(x) 0),
  #     env = globalenv()),
  #   nm = ref)
  
  res <- setNames(rep(quos(0), length(ref)), nm = ref)
  res <- utils::modifyList(
    res, x
  )
  
  structure(res, class = original_class)
}

#' @export
check_init.default <- function(x, ref) {
  
  if (! length(x) == length(ref)) {
    stop("Incorrect length in ", to_check, ".")
  }
  
  if (is.null(names(x))) {
    names(x) <- ref
  }
  
  if (! all(sort(names(x)) == sort(ref))) {
    stop("Some ", to_check, " names are incorrect.")
  }
  
  define_init_(as_quosures(
    lapply(x[ref], function(x) x),
    env = globalenv()))
}

check_inflow <- function(x, ...) {
  res <- check_init(x, ...)
  structure(
    res,
    class = c("uneval_inflow", class(res)))
}

check_starting_values <- function(x, ...) {
  res <- check_init(x, ...)
  structure(
    res,
    class = c("starting_values", class(res)))
}
