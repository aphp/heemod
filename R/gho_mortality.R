#' Use WHO Mortality Rate
#' 
#' 
#' Returns age and sex-specific mortality probabilities for 
#' a given region
#' 
#' Only locally cached data are available. 
#' For memory space reasons
#' local data is only available for WHO high-income 
#' countries (pooled), and only for the latest year.
#' 
#' 
#' @name who_mortality
#' @param age age as a continuous variable.
#' @param sex sex as `"FMLE"`-`"MLE"`, `0`-`1` (male = 0,
#'   female = 1) or `1`-`2` (male = 1, female = 2).
#' @param region Region code. 
#' @param year Use data from that year. Defaults to 
#'   `"latest"`.

#'   
#' @return This function should be used within 
#'   [define_transition()] or [define_parameters()].
#'   
#' @examples 
#' 
#' define_transition(
#'   C, get_who_mr(age = 50 + model_time, sex = "FMLE", region = "EUR"),
#'   0, 1
#' )
#' 
get_who_mr <- function(age, sex = NULL, region = NULL, 
                            year = "latest") {
    message("Fetching mortality data from package cached data.")
    mr_data <- get_package_mr(
      year = as.character(year)
    )
  ages_gho <- trans_age_gho(age)
  age_gho <- if ("AGE100+" %in% mr_data$AGEGROUP){
    ages_gho[[1]]
  } else {
    ages_gho[[2]]  
  }
  ref_data <- tibble(
    AGEGROUP = age_gho
  )  
  
  if (!is.null(sex)) {
    ref_data$SEX <- trans_sex_gho(sex)
  } else {
    ref_data$SEX <- "BTSX"
  }
  
  ref_data$SEX <- paste0("SEX_", ref_data$SEX)
  ref_data$AGEGROUP <- paste0("AGEGROUP_", ref_data$AGEGROUP)
  
  
  suppressMessages({
    dplyr::left_join(
      ref_data,
      mr_data
    )$NumericValue
  })
}

get_package_mr <- function(year) {
  if (year != "latest" && year != list_morta$data$year) {
    stop(sprintf(
      "No local data available for specified year (specified: %s, available: %s).",
      year,
      list_morta$data$year
    ))
  }
  message(sprintf(
    "Using cached data from year %s.",
    list_morta$year
  ))
  
  list_morta$data
}

trans_age_gho <- function(age) {
  stopifnot(! is.null(age))
  stopifnot(
    age >= 0,
    is.numeric(age),
    ! any(is.na(age)),
    length(age) > 0
  )
  labs <- list(
    c(
      "AGELT1", "AGE1-4", "AGE5-9",
      "AGE10-14", "AGE15-19", "AGE20-24",
      "AGE25-29", "AGE30-34", "AGE35-39",
      "AGE40-44", "AGE45-49", "AGE50-54",
      "AGE55-59", "AGE60-64", "AGE65-69",
      "AGE70-74", "AGE75-79", "AGE80-84",
      "AGE85-89", "AGE90-94", "AGE95-99",
      "AGE100+"
    )
  )
  labs[[2]] <- c(labs[[1]][1:18], "AGE85PLUS")
  
  breaks <- list(
    c(0, 1, seq(5, 100, 5), +Inf),
    c(0, 1, seq(5, 85, 5), +Inf)
  )
  
  lapply(1:2, function(i){
    cut(
      age,
      breaks[[i]],
      right = FALSE,
      labels = labs[[i]]
    ) %>% 
      as.character()
  })
}

trans_sex_gho <- function(sex) {
  stopifnot(! is.null(sex))
  u_sex <- sort(unique(sex))
  
  if (length(u_sex) > 2)
    stop("More than 2 sex modalities.")
  
  stopifnot(
    all(u_sex %in% c("FMLE", "MLE")) ||
      all(u_sex %in% 0:1) ||
      all(u_sex %in% 1:2),
    ! any(is.na(sex)),
    length(sex) > 0
  )
  
  if (all(u_sex == 1)) {
    stop("All values of 'sex' are equal to 1, this format is ambiguous.")
  }
  
  if (all(u_sex %in% c("FMLE", "MLE"))) {
    return(as.character(sex))
    
  } else if (all(u_sex %in% 0:1) || all(u_sex %in% 1:2)) {
    message("Converting sex values (male = 0, female = 1) or (male = 1, female = 2).")
    sex %>% 
      factor(labels = c("MLE", "FMLE")) %>% 
      as.character() %>% 
      return()
    
  } else {
    stop("Error during conversion of labels for 'sex'.")
  }
}
