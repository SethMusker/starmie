# test-loadAdmixture.R
context("I/O checks for loadAdmixture")

# try some q and p and logfiles files
q3 <- system.file("extdata/hapmap3_files", "hapmap3.3.Q", package="starmie")
p3 <- system.file("extdata/hapmap3_files", "hapmap3.3.Q", package="starmie")
log3_fin <- system.file("extdata/hapmap3_files", "log3.out", package="starmie")
p2 <- system.file("extdata/hapmap3_files", "hapmap3.2.P", package="starmie")
log1_fin <-system.file("extdata/hapmap3_files", "log1.out", package="starmie")

q3_bias <- system.file("extdata/hapmap3_files", "hapmap3.3.Q_bias", package="starmie")
q3_se <- system.file("extdata/hapmap3_files", "hapmap3.3.Q_se", package="starmie")

test_that("loadAdmixture I/O checks", {
  expect_error(loadAdmixture(c("a", "b"), "c"), "qfile must be a character vector of length 1.")
  expect_error(loadAdmixture(p2, c("a", "b")), "pfile must be a character vector of length 1.")
  expect_error(loadAdmixture(q3, p3, c("a", "b")), "logfile must be a character vector of length 1.")
  expect_error(loadAdmixture(q3, p2),
               "Number of populations does not match between Q and P files")
  expect_error(loadAdmixture(q3, p3, log1_fin), "logfile does not match input P and Q files")
  expect_s3_class(loadAdmixture(q3,p3, log3_fin), "admix")
  expect_false(is.null(loadAdmixture(q3, p3, bootstrap_bias = q3_bias)$Q_bias_df))
  expect_false(is.null(loadAdmixture(q3, p3, bootstrap_se = q3_se)$Q_se_df))

})


