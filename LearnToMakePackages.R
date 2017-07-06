
# USEFUL HOTKEYS ----------------------------------------------------------

# ALT + SHIFT + K                                                                hotkeys
# CTRL + SHIFT + L                                                               load all, run all, and save all open files
# click func in code, press F2, press 'CTRL + .' then start typing the name.     find function in packages and R/ directory, show code
# CTRL + F9                                                                      back - and open closed tab
# CTRL + SHIFT + D                                                               convert roxygen comments to .Rd files
# CTRL + SHIFT + B                                                               build & reload package (slow but thorough) (doesn't 'load all')
# CTRL + SHIFT + /                                                               wrap comments so they are less than 80 char per line
# CTRL + SHIFT + K                                                               knit vignette
# CTRL + SHIFT + T                                                               test code with testthat

# USEFUL COMMANDS (STRUCTURE AND R/) ----------------------------------------
library(devtools)                    # loads and attaches devtools to search path (devtools::<method> can automatically load without attaching)
library(testthat)
library(roxygen2)
#library("package")                  # throws an error if package isn't found in lapply(.libPaths(), dir)
#require("package")                  # prints warning and returns false instead.
#devtools::session_info()            # R version, loaded packages, etc
#devtools::has_devel()               # check if you're ready to develop
#file.exists("~/.ssh/id_rsa.pub")    # check for RSA public key (git)
#devtools::create("Rcosmo")
#dir(full.names = TRUE)              # like 'ls' for the working directory
#install.packages()                  # downloads and installs a binary from CRAN.
#devtools::install_github()          # e.g. "hadley/devtools" (downloads a source (not binary) package from GitHub, runs build() to make vignettes, and then uses R CMD INSTALL to do the install.)
#devtools::build()                   # build bundle package (wrapper for R CMD build)
#devtools::use_build_ignore("notes") # exclude specific file or directory "notes" from source -> bundle (adds Perl regular expression to .Rbuildignore (.Rinstignore is equiv for bundle -> binary))
#devtools::build(binary = TRUE)      # build binary package
#devtools::install()                 # wrapper for R CMD INSTALL (install from binary or bundle)
#devtools::install_url(), devtools::install_gitorious(), devtools::install_bitbucket() # work similar to packages found elsewhere on net.
#devtools::load_all()                # reload your code: skip install and load package directly to memory. Use instead of source() to avoid changing landscape.
#.libPaths()                         # see which libraries are currently active
#lapply(.libPaths(), dir)            # see all packages installed in all active libraries (not necessarily loaded or on search path)

#formatR::tidy_dir("R")              # clean up R code (doesn't match my comment style on this file though).
#lintr::lint_package()               # highlight style issues rather than fixing!

#old <- options(stringsAsFactors = FALSE) # save old global options and graphics par()
#on.exit(options(old), add = TRUE)        # restore old global options and graphics par() when function exits

#old <- setwd(tempdir())                  # save old working directory
#on.exit(setwd(old), add = TRUE)          # restore old working directory when function exits

#default.stringsAsFactors()               # return global option default for stringsAsFactors (careful using functions that assume the user has some default e.g. read.csv())


# .onAttach <- function(libname, pkgname) {         #.onAttach can alter the landscape when a package is attached
#   packageStartupMessage("Welcome to my package")
# }
#
#
# .onLoad <- function(libname, pkgname) {           #.onLoad alters landscape on load
# op <- options()
# op.devtools <- list(
#   devtools.path = "~/R-dev",                      # replace "devtools" with "myPackage"
#   devtools.install.args = "",
#   devtools.name = "Your name goes here",          # so getOption("devtools.name") can get the name of package author.
#   devtools.desc.author = '"First Last <first.last@example.com> [aut, cre]"',
#   devtools.desc.license = "What license is it under?",
#   devtools.desc.suggests = NULL,
#   devtools.desc = list()
# )
# toset <- !(names(op.devtools) %in% names(op))    # only set those options which not yet there (dont conflict with user options)
# if(any(toset)) options(op.devtools[toset])
#
# invisible()
# }
#
#
# MORE USES OF .onLoad AND .onAttach:
# rJava::.jpackage()      to use rJava to talk to a .jar file.
# Rcpp::loadRcppModules() to make C++ classes available as reference classes in R with Rcpp modules.
# tools::vignetteEngine() register vignette engines.
# REMEMBER TO CONSIDER .onUnload()ing

#stringi::stri_escape_unicode(x)                  # convert unicode in string x to correct escape characters for CRAN (otherwise ASCII is fine)




# USEFUL COMMANDS (PACKAGE METADATA) --------------------------------------

# devtools::use_package("dplyr")              # Places package in DESCRIPTION under imports
# devtools::use_package("dplyr", "Suggests")  # Places under suggests.

# TWO OPTIONS IF USER DOESN'T HAVE SUGGESTED PACKAGE:
# # You need the suggested package for this function
# my_fun <- function(a, b) {
#   if (!requireNamespace("pkg", quietly = TRUE)) {
#     stop("Pkg needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
# }
#
# # There's a fallback method if the package isn't available
# my_fun <- function(a, b) {
#   if (requireNamespace("pkg", quietly = TRUE)) {
#     pkg::f()
#   } else {
#     g()
#   }
# }




# USEFUL COMMANDS (DOCUMENTATION) -----------------------------------------
# devtools::document()        convert roxygen comments to .Rd files.



# USEFUL COMMANDS (VIGNETTES) ---------------------------------------------

#devtools::build()                  #rebuilds package including vignettes

#devtools::install_github(build_vignettes = TRUE)    #installs and builds all vignettes and install all suggested packages

#browseVignettes()
#browseVignettes("packagename")

#devtools::use_vignette("my-vignette")

#knitr::opts_chunk$set(              # run this from a knitr block in vignette to set options for all blocks.
#  opt1 = val1,
#  opt2 = val2
#)
#
#OPTIONS
#eval = FALSE
#echo = FALSE
#results = "hide"
#warning = FALSE
#error = TRUE, purl = FALSE
#collapse = TRUE
#results = 'asis'
#fig.show = "hold"
#fig.width = 5, fig.height = 5




# USEFUL COMMANDS (TESTING) -----------------------------------------------

# devtools::use_testthat()               set up package for using testthat

#  ------------------------------------------------------------------------
devtools::use_build_ignore("cmbstatLearn")
