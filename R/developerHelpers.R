
contains <- function(word,string)
{
  return(isTRUE(grep(word, string) == 1))
}

# Check that a character parameter is one of the
# strings from ... and if not then output an error message
checkParam <- function(param, ...)
{
  # to be completed
}
