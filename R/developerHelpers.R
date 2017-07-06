#' Tabular Function For Rcosmo Developers
#'
#' Turns an R data frame into the correct format for a roxygen2 table.
#'
#' @param df A data frame
#' @param ... Arguments passed to \code{\link{format}}
#' @examples
#' cat(tabular(mtcars[1:5,1:5]))
#' @seealso
#' Courtesy of Hadley Wickham, \href{http://r-pkgs.had.co.nz/man.html}{here}.
tabular <- function(df, ...) {
 stopifnot(is.data.frame(df))

  align <- function(x) if (is.numeric(x)) "r" else "1"
  col_align <- vapply(df, align, character(1))

  cols <- lapply(df, format, ...)
  contents <- do.call("paste",
    c(cols, list(sep = " \\tab", collapse = "\\cr\n  ")))

  paste("\\tabular{", paste(col_align, collapse = ""), "}{\n ",
    contents, "\n}\n", sep = "")
}
