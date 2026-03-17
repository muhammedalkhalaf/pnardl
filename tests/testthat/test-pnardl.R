test_that("pnardl runs and returns expected structure", {
  set.seed(1)
  N <- 4; TT <- 25
  df <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y    = rnorm(N * TT),
    x    = rnorm(N * TT)
  )

  res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
                sr = character(0), asymmetric = "x",
                id = "id", time = "time", p = 1, asy_test = TRUE)

  expect_s3_class(res, "pnardl")
  expect_true(is.numeric(res$coefficients))
  expect_true(length(res$panel_coefs) > 0)
  expect_true("x_pos" %in% names(res$partial_sums))
  expect_true("x_neg" %in% names(res$partial_sums))
})

test_that("partial sums have correct signs", {
  set.seed(2)
  N <- 2; TT <- 20
  df <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y    = rnorm(N * TT),
    x    = rnorm(N * TT)
  )
  res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
                sr = character(0), asymmetric = "x",
                id = "id", time = "time", p = 1, asy_test = FALSE)

  # positive partial sums should be non-decreasing
  for (pid in unique(df$id)) {
    ps <- res$partial_sums[res$partial_sums$id == pid, "x_pos"]
    ps <- ps[!is.na(ps)]
    expect_true(all(diff(ps) >= -1e-10))
  }
})

test_that("print.pnardl does not error", {
  set.seed(3)
  N <- 3; TT <- 20
  df <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y    = rnorm(N * TT),
    x    = rnorm(N * TT)
  )
  res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
                sr = character(0), asymmetric = "x",
                id = "id", time = "time", p = 1, asy_test = FALSE)
  expect_output(print(res))
})

test_that("summary.pnardl returns list invisibly", {
  set.seed(4)
  N <- 3; TT <- 20
  df <- data.frame(
    id   = rep(1:N, each = TT),
    time = rep(1:TT, N),
    y    = rnorm(N * TT),
    x    = rnorm(N * TT)
  )
  res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
                sr = character(0), asymmetric = "x",
                id = "id", time = "time", p = 1, asy_test = FALSE)
  s <- summary(res)
  expect_true(is.list(s))
})
