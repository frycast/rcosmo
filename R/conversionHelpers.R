

# THE REMAINING FUNCTIONS SHOULD NOT BE MADE AVAILABLE TO PACKAGE USERS

# Function takes longitude in [0,2pi] and transforms it
# to longitude in [-pi,pi]
lonWrap <- function(lon) {
  (lon + 180) %% 360 - 180
}