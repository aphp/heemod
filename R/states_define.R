#' Define a Markov Model State
#' 
#' Define the values characterising a Markov Model state for
#' 1 cycle.
#' 
#' As with [define_parameters()], state values are
#' defined sequentially. Later state definition can thus
#' only refer to values defined earlier.
#' 
#' For the `modify` function, existing values are 
#' modified, no new values can be added. Values order 
#' matters since only values defined earlier can be 
#' referenced in later expressions.
#' 
#' @param ... Name-value pairs of expressions defining state
#'   values.
#' @param starting_values Optional starting values defined
#'   with [define_starting_values()].
#' @param .OBJECT An object of class `state`.
#' @param x Used to work around non-standard evaluation.
#'   
#' @return An object of class `state` (actually a named
#'   list of quosures).
#' @export
#' 
#' @example inst/examples/example_define_state.R
#'   
define_state <- function(..., starting_values = define_starting_values()) {
  .dots <- quos(...)
  deprecated_x_cycle(.dots)
  define_state_(list(.dots = .dots, starting_values = starting_values))
}

#' @export
#' @rdname define_state
define_state_ <- function(x) {
  .dots <- x$.dots
  check_names(names(.dots))
  starting_values <- check_starting_values(
    x = x$starting_values,
    ref = names(.dots)
  )
  structure(list(
            .dots = .dots, 
            starting_values = starting_values),
            class = c("state", "list"))
}

#' @export
#' @rdname define_state
modify.state <- function(.OBJECT, ...) {
  .dots <- quos(...)
  
  modify_(.OBJECT = .OBJECT, .dots = .dots)
}

#' @export
modify_.state <- function(.OBJECT, .dots) {
  check_names(names(.dots))
  
  if ("starting_values" %in% names(.dots)) {
    starting_values <- check_starting_values(
      x = eval_tidy(.dots$starting_values),
      ref = names(.OBJECT$.dots)
    )
    .OBJECT <- utils::modifyList(.OBJECT, list(starting_values = starting_values))
    .dots <- utils::modifyList(.dots, list(starting_values = NULL))
  }
  
  if (!length(.dots)){
    return(.OBJECT)
  }
  if (! all(names(.dots) %in% names(.OBJECT$.dots))) {
    stop(sprintf(
      "The following state values are not defined: %s.",
      names(.dots)[names(.dots) %in% names(.OBJECT$.dots)]
    ))
  }
  
  utils::modifyList(.OBJECT$.dots, .dots)
}

#' Define Markov Model State List
#' 
#' Define the states of a Markov model by combining 
#' `state` objects.
#' 
#' State names have to correspond to those specified through
#' [define_transition()].
#' 
#' All states should have the same value names.
#' 
#' The `modify` function can modify existing states or 
#' add new ones.
#' 
#' @param ... Name-value pairs of expressions defining model
#'   states.
#' @param .OBJECT An `uneval_states` object.
#' @param .dots List of states, only used by 
#'   `define_state_list_` to avoid using `...`.
#'   
#' @return An object of class `uneval_state_list` (a 
#'   list of `state` objects).
#'   
#' @examples
#' \dontrun{
#' s1 <- define_state(cost = 1, util = 1)
#' s2 <- define_state(cost = 3, util = .4)
#' 
#' states_mod <- define_state_list(
#'   healthy = s1,
#'   sick = s2
#' )
#' 
#' states_mod
#' 
#' s1_bis <- define_state(cost = 0, util = 1)
#' s3 <- define_state(cost = 10, util = .1)
#' 
#' modify(
#'   states_mod,
#'   healthy = s1_bis,
#'   sicker = s3
#' )
#' }
#'   
#' @keywords internal
define_state_list <- function(...) {
  .dots <- list(...)
  define_state_list_(.dots)
}

#' @rdname define_state_list
define_state_list_ <- function(.dots) {
  # states <- lapply(.dots, function(x){
  #   structure(
  #     x$.dots,
  #     class = class(x)
  #   )
  # })
  states <- .dots
  state_names <- names(states)
  
  if (is.null(state_names)) {
    if (!identical(Sys.getenv("TESTTHAT"), "true"))
      message("No named state -> generating names.")
    state_names <- LETTERS[seq_along(states)]
    names(states) <- state_names
  }
  
  if (any(state_names == "")) {
    warning("Not all states are named -> generating names.")
    state_names <- LETTERS[seq_along(states)]
    names(states) <- state_names
  }
  
  if (any(duplicated(names(states)))) {
    stop("Some state names are duplicated.")
  }
  
  if (! all(unlist(lapply(states,
                          function(x) "state" %in% class(x))))) {
    
    .x <- names(states)[! unlist(lapply(
      states,
      function(x) "state" %in% class(x)))]
    
    stop(sprintf(
      "Incorrect state object%s: %s",
      plur(length(.x)),
      paste(.x, collapse = ", ")
    ))
  }
  
  check_states(states)
  structure(
    states,
    class = c("uneval_state_list", class(states))
  )
}

#' @rdname define_state_list
#' @export
modify.uneval_state_list <- function(.OBJECT, ...) {
  .dots <- list(...)
  
  modify_(.OBJECT = .OBJECT, .dots = .dots)
}

#' @export
modify_.uneval_state_list <- function(.OBJECT, .dots) {
  res <- utils::modifyList(.OBJECT, .dots)
  check_states(res)
  
  res
}

#' Check Model States for Consistency
#' 
#' For internal use.
#' 
#' All states should have the same value names.
#' 
#' @param x An object of class `uneval_states`.
#'   
#' @return `NULL`
#'   
#' @keywords internal
check_states <- function(x){
  if (! list_all_same(lapply(x, function(y) length(y$.dots)))) {
    stop("Number of state values differ between states.")
  }
  
  if (! list_all_same(lapply(x, function(y) sort(names(y$.dots))))) {
    stop("State value names differ between states.")
  }
  NULL
}

#' Return Number of State
#' 
#' For internal use.
#' 
#' Work with both `uneval_states` and
#' `eval_states`.
#' 
#' @param x An object containing states.
#'   
#' @return An integer: number of states.
#'   
#' @keywords internal
get_state_number <- function(x){
  # !mod!
  # rename get_state_count
  length(get_state_names(x))
}

#' Return Names of State Values
#' 
#' @param x An object containing states.
#'   
#' @return A character vector of state value names.
#'   
#' @keywords internal
get_state_value_names <- function(x){
  UseMethod("get_state_value_names")
}

#' @export
get_state_value_names.uneval_state_list <- function(x) {
  names(x[[1]][[1]])
}

#' @export
get_state_value_names.state <- function(x){
  names(x[[1]])
}

#' Get State Names
#' 
#' Retrieve state names from an object containing states.
#' 
#' @param x An object containing states.
#' @param ... Additional arguments passed to methods.
#'   
#' @return A character vector of state names.
#'   
#' @keywords internal
get_state_names <- function(x, ...){
  UseMethod("get_state_names")
}

#' @export
get_state_names.eval_state_list <- function(x, ...){
  names(x$.dots)
}

#' @export
get_state_names.default <- function(x, ...){
  names(x)
}
