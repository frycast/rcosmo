testthat::test_that("downloaders require explicit destinations before network access", {
  testthat::expect_error(downloadCMBMap(), "destfile must be an explicit")
  testthat::expect_error(downloadCMBPS(), "destfile must be an explicit")
  for (dest in list("", NA_character_, character(), c("a", "b"))) {
    testthat::expect_error(downloadCMBMap(destfile = dest), "destfile must be an explicit")
    testthat::expect_error(downloadCMBPS(destfile = dest), "destfile must be an explicit")
  }
})

testthat::test_that("invalid ring pixel indices fail instead of printing and continuing", {
  testthat::expect_error(ring2nest(1, 0), "pix must contain HEALPix indices")
  testthat::expect_error(ring2nest(1, 13), "pix must contain HEALPix indices")
})

testthat::test_that("plotAngDis restores graphical parameters on success and error", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(3, 3, 2, 1))
  oldpar <- graphics::par(c("mfrow", "mar"))
  cmbdf <- CMBDataFrame(nside = 8, I = sin(seq_len(12 * 8^2)), ordering = "nested")
  plotAngDis(cmbdf)
  testthat::expect_equal(graphics::par(c("mfrow", "mar")), oldpar)

  ## Non-finite intensities force a plotting error after the par() change.
  cmbdf$I <- rep(Inf, nrow(cmbdf))
  testthat::expect_error(suppressWarnings(plotAngDis(cmbdf)))
  testthat::expect_equal(graphics::par(c("mfrow", "mar")), oldpar)
})
