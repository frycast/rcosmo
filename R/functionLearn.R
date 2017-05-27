#' Add together two numbers
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @family blerbs.
#' @examples
#' add(1,1)
#' add(10,1)
#' \dontrun{add(1,2,3)}
add <- function(x,y) {
 x + y
}


#' Sum of vector elements.
#'
#' \code{sum} returns the sum of all the values present in its arguments
#'
#' This is a generic function: methods can be defined for it directly or via the
#' \code{\link{Summary}} group generic. For this to work properly, the arguments
#' \code{...} should be unnamed, and dispatch is on the first argument.
#'
#' @family blerbs
#' @aliases plus summation
#' @keywords methods misc
#' @param ... Numeric, complex or logical vectors.
#' @param na.rm A logical scalar. Should missing values (including NaN) be removed?
#' @return If all inputs are integer and logical, then the output
#'   will be an integer. If integer overflow
#'   \url{http://en.wikipedia.org.wiki/Integer_overflow}  occurs, the output
#'   will be NA with a warning. Otherwise it will be length-one numeric or complex vector.
#'
#'   Zero-length vectors have sum 0 by definition.
#'
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F,F,F,T,T)
#' \dontrun{sum("a")}
#' @section Practice Section:
#'   Just adding this text to see what roxygen does.
#' @seealso Look at this link in my package: \code{\link{add}}
#'
#'   and this link in the \code{\link{stats}} package: \code{\link[stats]{var}}
#'
#'   and my cool URL \url{www.sciencestems.com}.
sum <- function(..., na.rm = TRUE) {}







