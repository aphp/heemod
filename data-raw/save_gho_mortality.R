library(rgho)
library(heemod)
library(dplyr)

get_latest_morta <- function() {
  mr_data <- rgho::get_gho_data(
    code = "LIFE_0000000029",
    filter = list(
      WORLDBANKINCOMEGROUP = "WB_HI"
    )
  )
  
  if (nrow(mr_data) == 0) return(NULL)
  
  study_year <- max(mr_data$YEAR)
  mr_data_year <- mr_data[mr_data$YEAR == study_year, ]
  
  if (nrow(mr_data_year) != 57) {
    stop("Strange GHO mortality data.")
  }
  
  list(
    data = mr_data_year,
    year = study_year,
    pool = FALSE
  )
}

list_morta <- get_latest_morta()

usethis::use_data(
  list_morta,
  internal = TRUE,
  overwrite = TRUE
)
