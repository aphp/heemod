.onLoad <- function(libname, pkgname) {
  op <- options()
  op.heemod <- list(
    heemod.verbose = FALSE,
    heemod.memotime = 3600,
    heemod.inf_parameter = "stop"
  )
  toset <- !(names(op.heemod) %in% names(op))
  if(any(toset)) options(op.heemod[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  options("heemod.env" =  rlang::new_environment(parent = asNamespace("heemod")))
  assign("start_tibble", data.frame(), getOption("heemod.env"))
  invisible()
}

