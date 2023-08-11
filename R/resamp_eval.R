alert_psa_surv <- function(psa){
  list_vars <- map(psa$surv, function(x){
    var <- getOption("heemod.env")[[x]]
    if(!inherits(var, c("surv_dist", "surv_fit"))){
      cli::cli_abort("{.var {x}} must be a {.cls surv_fit} or {.cls surv_dist} object, i.e the initial survival \\
      distribution before the modification leading to the creation of a {.cls {class(var)[[1]]}} object")
    }
    })
}

remove_undefined_psa <- function(psa, param){
  list_vars <- map(param, all.vars) 
  all_vars <- unique(c(unlist(list_vars, use.names = F), 
                       names(list_vars), 
                       ls(getOption("heemod.env"))))
  not_found <- setdiff(names(psa$list_qdist), all_vars)
  if (!length(not_found)) return(psa)
    cli::cli_warn(glue::glue('{paste(not_found, collapse = ", ")} \\
                                      neither used by define_parameters nor found in the environment. \\
                                      Will be skipped.'))
    psa$list_qdist <- psa$list_qdist[-which(names(psa$list_qdist) %in% not_found)]
    if (!length(psa$list_qdist)) {
      cli::cli_abort("No defined parameter to resample")
    }
    psa$surv <- psa$surv[-which(psa$surv %in% not_found)]
    psa$multinom <- psa$multinom[-which(psa$multinom %in% not_found)]
    psa
}

#' Run Probabilistic Uncertainty Analysis
#' 
#' @param model The result of [run_model()].
#' @param psa Resampling distribution for parameters defined
#'   by [define_psa()].
#' @param N &gt; 0. Number of simulation to run.
#' @param keep logical; if TRUE, all models will be returned
#'   
#' @return A list with the following elements
#' * psa:  a `data.frame` with one row per model.
#' * run_model: a `data.frame` with mean cost and utility for each strategy
#' * model: the initial model object 
#' * N: the number of simulations ran
#' * resamp_par: the resampled parameters
#' * full: if `keep` is TRUE, a list of each model objects created at each iteration
#' 
#' @export
#' 
#' @example inst/examples/example_run_psa.R
#'   
run_psa <- function(model, psa, N, keep = FALSE) {
  stopifnot(
    N > 0,
    ! is.null(N)
  )
  if (! all(c(".cost", ".effect") %in% names(get_model_results(model)))) {
    stop("No cost and/or effect defined, probabilistic analysis unavailable.")
  }
  copy_param_env(model$parameters)
  on.exit(copy_param_env(model$parameters))
  psa <- remove_undefined_psa(psa, model$parameters)
  alert_psa_surv(psa)
  newdata <- eval_resample(psa, N)
  
  list_res <- list_full <- list()
  
  for (n in get_strategy_names(model)) {
    if (!identical(Sys.getenv("TESTTHAT"), "true"))
      message(sprintf("Resampling strategy '%s'...", n))
    e_newdata <- eval_strategy_newdata(
      x = model,
      strategy = n,
      newdata = newdata
    )
    message(sprintf("Resampling strategy '%s'...", n))
    list_res <- c(
      list_res,
      list(
        e_newdata
        %>% 
          dplyr::rowwise() %>% 
          dplyr::do(get_total_state_values(.data$.mod)) %>% 
          dplyr::bind_cols(newdata) %>% 
          dplyr::ungroup()
      )
    )
    if (keep){
      list_full <- c(
        list_full, 
        list(
          e_newdata$.mod
        )
      )
    }
  }
  
  names(list_res) <- get_strategy_names(model)
  if (keep) names(list_full) <- get_strategy_names(model)
  index <- seq_len(N)
  
  for (n in names(list_res)) {
    list_res[[n]]$.strategy_names <- n
    list_res[[n]]$.index <- index
  }
  
  x_tidy <- get_ce(model)
  
  res <- 
    dplyr::bind_rows(list_res)
  res <- dplyr::mutate(res, !!!x_tidy)
  
  run_model <- res %>% 
    dplyr::select(-.index) %>% 
    dplyr::group_by(.data$.strategy_names) %>%
    dplyr::select_if(is.numeric) %>%
    dplyr::summarise_all(mean) %>% 
    as.data.frame()
  
   # restore environment
  
  structure(
    c(list(
      psa = res,
      run_model = run_model[names(get_model_results(model))],
      model = model,
      N = N,
      resamp_par = names(newdata)
    ),
    if (keep)
      list(full = list_full)
    ),
    class = c("psa", class(list()))
  )
}

get_model <- function(x) {
  UseMethod("get_model")
}

#' @export
get_model.psa <- function(x) {
  x$model
}

#' @export
get_model_results.psa <- function(x) {
  x$run_model
}

#' @export
get_cycles.psa <- function(x) {
  get_cycles(get_model(x))
}

#' @export
get_uneval_init.psa <- function(x) {
  get_uneval_init(get_model(x))
}

#' @export
get_method.psa <- function(x) {
  get_method(get_model(x))
}

#' @export
get_central_strategy.psa <- function(x, ...) {
  get_central_strategy(get_model(x))
}

#' @export
get_noncomparable_strategy.psa <- function(x, ...) {
  get_noncomparable_strategy(summary(x, ...)$res_comp)
}

#' @export
get_root_strategy.psa <- function(x, ...) {
  get_root_strategy(summary(x, ...)$res_comp)
}

eval_correlation <- function(x, var_names) {
  res <- diag(length(var_names))
  colnames(res) <- var_names
  rownames(res) <- var_names
  
  for (i in seq_len(nrow(x))) {
    res[x$v1[i], x$v2[i]] <- x$cor[i]
    res[x$v2[i], x$v1[i]] <- x$cor[i]
  }
  res
}

#' Evaluate Resampling Definition
#' 
#' @param psa A [define_psa()] object.
#' @param N &gt; 0. Number of simulation to run.
#'   
#' @return A `data.frame` of resampled values with on 
#'   column per parameter and `N` rows.
#'   
#' @keywords internal
eval_resample <- function(psa, N) {
  list_qdist_filter <- psa$list_qdist[setdiff(names(psa$list_qdist), psa$surv)]
  res <- data.frame(._start = logical(N))
  if (length(list_qdist_filter)){
    mat_p <- stats::pnorm(
      mvnfast::rmvn(
        n = N,
        mu = rep(0, length(list_qdist_filter)),
        sigma = psa$correlation
      )
    )
  
    list_res <- mapply(
      function(i, f) f(mat_p[, i]),
      seq_len(ncol(mat_p)),
      list_qdist_filter
    )
    
    if (length(dim(list_res)) < 2) {
      list_res <- matrix(list_res, ncol = length(list_res))
    }
    
    colnames(list_res) <- names(list_qdist_filter)
    res <- as.data.frame(list_res)
    
    for (m in psa$multinom) {
      call_denom <- make_call(m, "+")
      list_expr <- as_quosures(
        c(list(
          .denom = call_denom),
        stats::setNames(
          lapply(
            m,
            function(x) as.call(list(as.name("/"), as.name(x), as.name(".denom")))),
          m)), env = parent.frame()) 
      res <- dplyr::mutate(res, !!!list_expr) %>% 
        dplyr::select(-.denom)
    }
    res$.denom <- NULL
  }
  
  loc <- list()
  for (s in psa$surv){
    x <- psa$list_qdist[[s]]
    lhs <- as.name(s)
    rhs <- rlang::call2(x[[1]], lhs, !!!rlang::call_args(x))
    res[[s]] <- replicate(N, eval(rhs, envir = getOption("heemod.env")), simplify = F) %>%
      map(function(.x){
        if (inherits(.x[[1]], "boot_surv")){
          tmp <- .x[[1]]()
          structure(list(function(x) tmp), class = c(class(.x), "surv_psa"))
        } else .x
      })
    #loc[[s]] <- 
  }
  res$._start <- NULL
  res
}


#' @keywords internal
compute_surv_ci <- function(x, times, type, psa, Nrep){
  resamples <- eval_resample(psa, Nrep)
  env <- rlang::env()
  res <- map(seq_len(nrow(resamples)), function(i){
    get_new_surv_parameters(map(resamples,i), env = env)
    suppressMessages(compute_surv_(eval_tidy(x, env = env), times, type = type, env = env))
  })
  unlist(res, use.names = FALSE) %>% 
    matrix(ncol = Nrep) %>% 
    apply(MARGIN = 1, FUN = function(x) c(quantile(x,c(0.025, 0.975), na.rm = TRUE)), simplify = T)%>%
    t() %>% 
    as.data.frame%>%
    dplyr::mutate(times = times)
}