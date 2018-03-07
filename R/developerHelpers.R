
contains <- function(word,string)
{
  return(isTRUE(grep(word, string) == 1))
}
