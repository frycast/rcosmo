library(devtools)
#library(pathological)
#os_path()
#example(startup)
#Sys.getenv() # show current environment variables
#help("environment variables") # details on important environment variables

# set Personal Access Token environment variable for the session
Sys.setenv("GITHUB_PAT" = "<Insert PIT here>")
# Or alternatively make it permanent by adding it to .Renviron:
# cat("GITHUB_PAT=<Insert PIT here>\n",
#     file = file.path(normalizePath("~/"), ".Renviron"), append = TRUE)
## Sys.getenv("GITHUB_PAT") # check that the PAT is now available.

# Install the package from repo
devtools::install_github("VidaliLama/cmbstat")

# Remove the PAT environment variable
Sys.unsetenv("GITHUB_PAT")
