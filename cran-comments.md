## Resubmission of archived package: rcosmo 1.1.5

rcosmo was archived on 2022-05-04 because its dependency geoR was archived.
geoR is available on CRAN again. This release retains that dependency and
updates documentation for current R and dependency versions.

## Test environments

GitHub Actions, 2026-09-30:

* Ubuntu 24.04.5 LTS, R-devel (2026-09-29 r90598)
* Ubuntu 24.04.5 LTS, R 4.6.1 and R 4.5.3
* Windows Server 2022, R 4.6.1
* macOS Tahoe 26.6.2 (arm64), R 4.6.1

## R CMD check results

Full --as-cran checks include the PDF manual and HTML/math validation.
There are no errors or warnings. Linux release/oldrel and Windows report
Status: OK. The R-devel incoming and macOS HTML checks report the NOTES
explained below.

The R-devel incoming NOTE identifies this as a new submission of an archived
package and reports HTTP 503 (Service Unavailable) responses from the ESA
Planck Legacy Archive links in covPwSp.Rd and downloadCMBPS.Rd. These are
external service responses; the original archiving reason is resolved by
the return of geoR to CRAN.

The macOS HTML manual NOTE reports that the runner's system HTML Tidy is
not recent enough for HTML validation. PDF generation and math rendering pass.
