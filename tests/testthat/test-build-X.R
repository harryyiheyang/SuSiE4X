test_that("build_X validates its user-facing inputs", {
  expect_error(
    build_X(data.frame(SNP = "rs1"), "chr_{CHR}"),
    "CHR and SNP"
  )
  expect_error(
    build_X(data.frame(CHR = "1", SNP = "rs1"), "chr_1"),
    "must contain \\{CHR\\}"
  )
  expect_error(
    build_X(data.frame(CHR = NA, SNP = "rs1"), "chr_{CHR}"),
    "no usable chromosome and SNP records"
  )
  expect_error(
    build_X(data.frame(CHR = "1", SNP = "rs1"), "chr_{CHR}",
            clumping = NA),
    "clumping must be TRUE or FALSE"
  )
  expect_error(
    build_X(data.frame(CHR = "1", SNP = "rs1"), "chr_{CHR}",
            clumping = TRUE, r2 = 1.01),
    "r2 must be a single number"
  )
})

test_that("build_X clumping defaults are stable", {
  expect_identical(formals(build_X)$clumping, FALSE)
  expect_identical(formals(build_X)$r2, 0.81)
})
