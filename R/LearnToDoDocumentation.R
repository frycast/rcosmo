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
#' @param na.rm A logical scalar. Should missing values (including NaN) be
#'   removed?
#' @return If all inputs are integer and logical, then the output will be an
#'   integer. If integer overflow
#'   \url{http://en.wikipedia.org.wiki/Integer_overflow} occurs, the output will
#'   be NA with a warning. Otherwise it will be length-one numeric or complex
#'   vector.
#'
#'   Zero-length vectors have sum 0 by definition.
#'
#' @examples
#' sum(1:10)
#' sum(1:5, 6:10)
#' sum(F,F,F,T,T)
#' \dontrun{sum("a")}
#' @section Practice Section: Just adding this text to see what roxygen does.
#' @seealso Look at this link in my package: \code{\link{add}}
#'
#'   and this link in the \code{\link{stats}} package: \code{\link[stats]{var}}
#'
#'   and my cool URL \url{www.sciencestems.com}.
sum <- function(..., na.rm = TRUE) {}


#' The Foo Function
#'
#' @param a This is the first argument
foo <- function(a) a + 10


#' The Bar Function
#'
#' @param b This is the second argument
#' @inheritParams foo
bar <- function(a,b) {
  foo(a)*b
}


#' Foo Bar generic
#'
#' The generic will \strong{foo} certain things, and then bar for
#' \emph{reasons}. Contact \email{daniel-fryer@@live.com.au} for more functions
#' like this if you
#' \enumerate{
#' \item Have money.
#' \item Are willing to spend it.
#' }
#' because
#' \deqn{a + b = \U222B &#x222b x dx}
#'
#' @param x Object to foo
foobar <- function(x) UseMethod("foobar")

#' @describeIn foobar Difference between the mean and the median.
foobar.numeric <- function(x) abs(mean(x) - median(x))

#' @describeIn foobar First and last values pasted together in a string.
foobar.character <- function(x) paste0(x[1], "-", x[length(x)])

#' @rdname foobar
#' @param y Numeric for thisIsOutOfPlace
thisIsOutOfPlace <- function(x,y) x
