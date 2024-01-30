#' Markov Models for Health Economic Evaluations
#' 
#' An implementation of the modelling and
#' reporting features described in reference
#' textbooks and guidelines: deterministic and
#' probabilistic sensitivity analysis, 
#' heterogeneity analysis, time dependency
#' on state-time and model-time (semi-Markov
#' and non-homogeneous Markov models), etc.
#' 
#' @keywords internal
"_PACKAGE"

## usethis namespace: start

#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr group_by
#' @importFrom dplyr as_tibble
#' @importFrom dplyr tibble
#' @importFrom dplyr as_data_frame
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr %>%
#' @importFrom dplyr desc
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate_if
#' @importFrom dplyr funs
#' @importFrom dplyr anti_join
#' 
#' @importFrom utils head
#' @importFrom utils modifyList
#' @importFrom utils globalVariables
#' @importFrom utils as.roman
#'   
#' @importFrom stats pnorm
#' @importFrom stats dist
#' @importFrom stats qbeta
#' @importFrom stats qbinom
#' @importFrom stats qgamma
#' @importFrom stats qlnorm
#' @importFrom stats qnorm
#' @importFrom stats terms
#' @importFrom stats setNames
#' @importFrom stats reorder
#' @importFrom stats na.omit
#' @importFrom stats update
#' @importFrom stats as.formula
#' @importFrom stats var
#' @importFrom stats coef
#' @importFrom stats model.matrix
#' @importFrom stats formula
#' @importFrom stats ecdf
#' @importFrom stats nls
#' @importFrom stats quantile
#' @importFrom stats plogis
#' @importFrom stats qlogis
#'   
#' @importFrom mvnfast rmvn
#'   
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 %+replace%
#'   
#' @importFrom memoise memoise
#' @importFrom memoise timeout
#'   
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom utils packageVersion
#'   
#' @importFrom tools file_ext
#' 
#' @importFrom grDevices dev.off
#' @importFrom grDevices cairo_pdf
#' @importFrom grDevices png
#' 
#' @importFrom graphics plot
#' @importFrom graphics par
#'   
#' @importFrom tibble tibble
#' 
#' @importFrom rlang sym syms quo .data quos quos as_quosures enexprs 
#' set_expr as_label as_quosure get_expr is_quosure expr eval_tidy new_quosure
#' parse_exprs get_env call_match caller_env enquo exprs enexpr %||%
#'
#' @importFrom purrr map map_chr map_dbl map_lgl
## usethis namespace: end
NULL

#' @export
dplyr::`%>%`

utils::globalVariables(".")