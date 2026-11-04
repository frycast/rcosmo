## rcosmo 1.1.6

This update addresses CRAN's report of Debian failures when the installed
user library is read-only, with a requested correction before 2026-11-22.
CMBDat() now explicitly requests PROT_READ instead of mmap's default
read/write protection. This also fixes FITS reads through CMBDataFrame().

A standalone regression test reads a temporary FITS file with Unix write
permissions removed, compares full and sampled reads, checks mapping
protection, exercises CMBDataFrame(), and verifies the file checksum.
The legacy testthat runner remains disabled; the new test runs independently.

## Check results

R CMD build followed by R CMD check --as-cran, GitHub Actions:
https://github.com/frycast/rcosmo/actions/runs/41784252574

All five jobs passed, including examples and the regression test:
* Ubuntu R-release and R-oldrel-1: Status OK.
* Ubuntu R-devel: one incoming NOTE about the short interval since the last
  release. This update is in response to CRAN's correction request.
* macOS R-release: one runner NOTE about outdated HTML Tidy.
* Windows R-release: passed.

Coverage CI also passed. The CI font-installation hang on Windows was fixed
by selecting TeX Live's bundled utilities.
