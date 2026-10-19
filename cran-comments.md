## Resubmission of archived package: rcosmo 1.1.5

rcosmo was archived on 2022-05-04 because its dependency geoR was archived.
geoR is available on CRAN again.

## Response to CRAN review

Thank you for reviewing the previous submission.

* Expanded Description into one paragraph covering the package's purpose,
  functionality and implemented methods; removed excess whitespace.
* Removed the redundant LICENSE file and its DESCRIPTION reference.
  The package remains licensed under GPL-3.
* Added the 43 missing return-value sections, including output types and
  meanings or side effects, and the missing replacement-method arguments.
* Replaced commented-out example code with small runnable examples.
  FITS examples use a bundled synthetic twelve-pixel map. Download examples
  use explicit temporary destinations and are guarded by if (interactive()).
* Replaced unconditional diagnostic printing with suppressible warnings or
  errors. Remaining console output is in print/summary methods or controlled
  by the existing messages argument.
* Removed implicit download destinations. downloadCMBMap requires destfile;
  downloadCMBPS requires destfile when saving. Examples use temporary files.
* plotAngDis saves graphics settings and immediately registers on.exit()
  restoration, including when an error occurs.

Added the package paper reference to DESCRIPTION in author-year/DOI format:
Fryer, Li and Olenko (2020) <doi:10.32614/RJ-2020-012>.

## Check results

GitHub Actions, 2026-10-19 (UTC), commit 7d52cc7:
https://github.com/frycast/rcosmo/actions/runs/39198943908

R CMD check --as-cran reports no errors or warnings:

* Ubuntu 24.04.5 LTS, R 4.6.1 and R 4.5.3: Status OK.
* Windows Server 2022, R 4.6.1: Status OK.
* Ubuntu 24.04.5 LTS, R-devel (2026-10-17 r90697): one incoming NOTE
  identifying a new submission of an archived package.
* macOS Tahoe 26.6.2 (arm64), R 4.6.1: one NOTE because the runner's
  HTML Tidy is too old for HTML validation; PDF generation passes.

These results include example checks. The existing tests/testthat.R runner
is commented out, so they do not establish that the testthat suite passes.
