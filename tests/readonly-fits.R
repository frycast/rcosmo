# Standalone regression test: the legacy testthat runner is disabled.
library(rcosmo)

local({
  source <- system.file("extdata", "CMB_testmap.fits", package = "rcosmo")
  filename <- tempfile(fileext = ".fits")
  stopifnot(file.copy(source, filename))
  on.exit({
    # Release the temporary mappings created internally by CMBDataFrame().
    gc()
    Sys.chmod(filename, "0600")
    unlink(filename)
  }, add = TRUE)

  before <- tools::md5sum(filename)
  expected <- CMBDat(filename)
  expected_df <- CMBDataFrame(filename, spix = 1:10)

  # mmap's Windows implementation still requests native write access even
  # with PROT_READ. Exercise permission denial on Unix, matching CRAN Debian.
  if (.Platform$OS.type == "unix") {
    stopifnot(Sys.chmod(filename, "0444"))
    # Root can bypass file modes; the mapping protection is checked below too.
    if (Sys.info()[["effective_user"]] != "root") {
      stopifnot(file.access(filename, 2L) == -1L)
    }
  }

  stopifnot(identical(CMBDat(filename), expected))
  stopifnot(identical(CMBDat(filename, spix = 1:10)$data,
                      expected$data[1:10, , drop = FALSE]))

  mapped <- CMBDat(filename, mmap = TRUE)
  on.exit(mmap::munmap(mapped$data), add = TRUE, after = FALSE)
  stopifnot(mmap::is.mmap(mapped$data))
  stopifnot(identical(mapped$data$prot, mmap::mmapFlags("PROT_READ")))
  stopifnot(identical(mapped$data[1:12], expected$data))
  stopifnot(identical(CMBDataFrame(filename, spix = 1:10), expected_df))
  stopifnot(identical(tools::md5sum(filename), before))
})
