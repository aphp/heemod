#' Iteratively Evaluate a Markov Model With New Parameter 
#' Values
#' 
#' Given a data.frame with on set of new parameters values 
#' per row, iteratively evaluate the model over the set of 
#' new values.
#' 
#' New parameters with a missing value (`NA`) do not 
#' replace existing parameters.
#' 
#' @param x Result from [run_model()].
#' @param strategy Name or index of model to recompute.
#' @param newdata a data.frame whose names match parameters 
#'   names. `strategy` will be evaluated iteratively, 
#'   taking successive values from each row.
#'   
#' @return A data.frame containing the values of 
#'   `newdata` and each Markov Model evaluation in 
#'   `res`.
#'   
#' @example inst/examples/example_eval_strategy_newdata.R
#'   
#' @keywords internal
eval_strategy_newdata <- function(x, strategy = 1, newdata) {
  strategy <- check_strategy_index(x = x, i = strategy)
  cycles <- get_cycles(x)
  init <- get_uneval_init(x)
  inflow <- get_inflow(x)
  method <- get_method(x)
  old_parameters <- get_parameters(x)
  uneval_strategy <- x$uneval_strategy_list[[strategy]]
  expand_limit <- get_expand_limit(x, strategy)
  
  if (status_cluster(verbose = FALSE)) {
    cl <- get_cluster()
    
    num_cores <- length(cl)
    if (!identical(Sys.getenv("TESTTHAT"), "true"))
    message(paste("Using a cluster with", num_cores, "cores."))
    
    split_vec <- rep(1:num_cores, each = nrow(newdata) %/% num_cores)
    split_vec <- c(split_vec, rep(num_cores, nrow(newdata) %% num_cores))
    
    pnewdata <- split(newdata, split_vec)
    parallel::clusterExport(
      cl, 
      c("uneval_strategy", "old_parameters", "pnewdata", 
        "cycles", "init", "method", "strategy"),
      envir = environment()
    )
    suppressMessages(
      pieces <- parallel::parLapply(cl, pnewdata, function(newdata) {
        newdata %>% 
          dplyr::rowwise() %>% 
          dplyr::do(
            .mod = eval_newdata(
              .data,
              strategy = uneval_strategy,
              old_parameters = old_parameters,
              cycles = cycles,
              init = init,
              inflow = inflow,
              method = method,
              strategy_name = strategy,
              expand_limit = expand_limit
            )
          ) %>% 
          dplyr::ungroup() %>% 
          dplyr::bind_cols(
            newdata
          )
      })
    )
    res <- dplyr::bind_rows(pieces)
    rownames(res) <- NULL
    
  } else {
    suppressMessages(
      res <- newdata %>% 
        dplyr::rowwise() %>% 
        dplyr::do(
          .mod = eval_newdata(
            .data,
            strategy = uneval_strategy,
            old_parameters = old_parameters,
            cycles = cycles,
            init = init,
            method = method,
            inflow = inflow,
            strategy_name = strategy,
            expand_limit = expand_limit
          )
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::bind_cols(
          newdata
        )
    )
  }
  res
}

get_new_surv_parameters <- function(new_parameters, env = getOption("heemod.env")){
  surv_new_parameters <- Filter(function(x) inherits(x, "surv_psa"), new_parameters)
  
  new_list <- map(surv_new_parameters,1)
  
  quo_surv <- Filter(is_quosure, new_list)
  non_quo_surv <- new_list[setdiff(names(new_list), names(quo_surv))]
  
  nm <- map(quo_surv, function(x){
    deparse(get_expr(x))
  })
  data_surv <- map(seq_along(quo_surv), function(i){
    get_env(quo_surv[[i]])[[nm[[i]]]]
  }) %>% setNames(nm)
  
  list2env(c(non_quo_surv, data_surv), envir = env)
  
  surv_new_parameters
}

eval_newdata <- function(new_parameters, strategy, old_parameters,
                         cycles, init, method, inflow,
                         strategy_name, expand_limit) {
  
  new_parameters <- Filter(
    function(x) all(rlang::is_call(x) || ! is.na(x)),
    new_parameters
  )
  
  new_parameters <- setdiff(new_parameters, get_new_surv_parameters(new_parameters))
  
  tidy_new_param <- to_dots(new_parameters)
  
  parameters <- utils::modifyList(
    old_parameters,
    tidy_new_param
  )
  
  eval_strategy(
    strategy = strategy,
    parameters = parameters,
    cycles = cycles,
    init = init,
    method = method,
    inflow = inflow,
    strategy_name = strategy_name,
    expand_limit = expand_limit
  )
}
