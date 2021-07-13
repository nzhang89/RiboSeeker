# utilities

#' Get current time
#'
#' @return A character variable of the current time.
#'
.now = function() {
  return(sprintf('[%s]', Sys.time()))
}
