
contains <- function(word,string)
{
  return(isTRUE(grep(word, string) == 1))
}


copyCMBAttributes <- function(cmbdf1, cmbdf2, exclude)
{

  atts <- c("ordering", "nside", "coords", "window", "row.names",
            "resolution", "header1", "header2")

  if ( !missing(exclude) )
  {
    atts <- atts[-which(atts %in% exclude)]
  }

  for ( a in atts )
  {
    attr(cmbdf1, a) <- attr(cmbdf2, a)
  }

  return(cmbdf1)
}
