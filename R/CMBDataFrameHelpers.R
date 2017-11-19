# Function returns HEALPix pixels from a CMBDF -------------------------------------------------
pix <- function(cmbdf)
{
  as.numeric( row.names(cmbdf) )
}



# THIS FUNCTION NEEDS FINISHING
# Function returns ordering scheme from a CMBDF and can convert ordering scheme ------------------
ordering <- function( cmbdf, newOrdering = NA )
{
  if ( is.na(newOrdering) )
  {
    attr( cmbdf, "ordering" )
  }
  else
  {
    if ( attr(cmbdf, "ordering") == newOrdering )
    {
      stop("Ordering is already " + newOrdering )
    }
    else if ( newOrdering == "nested" )
    {
      # Convert 'ring' to 'nested' ordering
    }
    else if ( newOrdering == "ring" )
    {
      # Convert 'nested' to 'ring' ordering
    }
  }
}




# Function returns Nside from a CMBDF -----------------------------------------------------------
Nside <- function( cmbdf )
{
  attr( cmbdf, "Nside" )
}




# THIS FUNCTION NEEDS FINISHING
# Function returns coordinate system of a CMBDF and can convert coordinates ----------------------
# if newCoords argument is passed
coords <- function( cmbdf, newCoords = NA )
{
  # Check that argument is a CMBDF
  if ( !("CMBDataFrame" %in% class(cmbdf)) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  
  # If newCoords argument not given then return the coordinate type
  if ( is.na(newCoords) )
  {
    attr(cmbdf, "coords")
  } 
  else 
  {
    # Make sure that newCoords doesn't match current coords
    if ( attr(cmbdf, "coords") == newCoords )
    {
      stop("Coordinates are already " + newCoords)
    }
    else if ( newCoords == "spherical" )
    {
      # Convert to spherical
    }
    else if ( newCoords == "cartesian" )
    {
      # Convert to cartesian
    }
  }
}



