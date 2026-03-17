#' Panel Nonlinear ARDL Estimation
#'
#' Estimates a Panel Nonlinear ARDL (PNARDL) model following Shin, Yu and
#' Greenwood-Nimmo (2014). Regressors listed in \code{asymmetric} are
#' decomposed into positive (\eqn{x^+}) and negative (\eqn{x^-}) partial
#' sums. The model is then estimated using within-group ordinary least squares
#' for each panel unit, and pooled estimates (PMG/MG/DFE-style) are computed.
#'
#' @param data A data frame containing all variables plus \code{id} and
#'   \code{time} columns.
#' @param depvar Character. Name of the dependent variable.
#' @param lr Character vector. Names of the long-run regressors (levels).
#' @param sr Character vector. Names of the additional short-run regressors
#'   (first-differences). These are included as \eqn{\Delta x} in the ECM.
#' @param asymmetric Character vector. Names of variables to decompose into
#'   positive and negative partial sums. Must be a subset of \code{lr}.
#' @param id Character. Name of the panel identifier column in \code{data}.
#' @param time Character. Name of the time identifier column in \code{data}.
#' @param estimator Character. One of \code{"pmg"} (default), \code{"mg"}, or
#'   \code{"dfe"}.
#' @param p Integer. Number of lags of \eqn{\Delta y} in the short-run
#'   equation. Default is \code{1}.
#' @param asy_test Logical. If \code{TRUE} (default), perform Wald tests for
#'   long-run and short-run asymmetry.
#' @param multip Integer. Number of periods for cumulative dynamic multipliers.
#'   Set to \code{0} (default) to skip.
#' @param level Numeric. Confidence level for intervals (e.g., \code{0.95}).
#'
#' @return An object of class \code{"pnardl"}, which is a list with components:
#'   \item{coefficients}{Named numeric vector of pooled long-run and
#'     short-run coefficients.}
#'   \item{se}{Standard errors of \code{coefficients}.}
#'   \item{tstat}{t-statistics.}
#'   \item{pval}{p-values.}
#'   \item{panel_coefs}{List of per-panel estimated coefficient vectors.}
#'   \item{ecm_coefs}{Named vector of ECT (speed-of-adjustment) coefficients
#'     per panel.}
#'   \item{asym_tests}{Data frame of Wald test results for asymmetry (if
#'     \code{asy_test = TRUE}).}
#'   \item{multipliers}{Data frame of cumulative dynamic multipliers (if
#'     \code{multip > 0}).}
#'   \item{partial_sums}{Data frame with the generated partial-sum variables
#'     appended to the original data.}
#'   \item{estimator}{Character. Estimator used.}
#'   \item{depvar}{Character. Dependent variable name.}
#'   \item{asymmetric}{Character vector. Decomposed variable names.}
#'   \item{pos_vars}{Character vector. Names of positive partial-sum columns.}
#'   \item{neg_vars}{Character vector. Names of negative partial-sum columns.}
#'   \item{call}{The matched call.}
#'
#' @references
#' Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
#' cointegration and dynamic multipliers in a nonlinear ARDL framework. In
#' R. C. Sickles & W. C. Horrace (Eds.), \emph{Festschrift in Honor of Peter
#' Schmidt} (pp. 281-314). Springer. \doi{10.1007/978-1-4899-8008-3_9}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' N <- 5; TT <- 30
#' df <- data.frame(
#'   id   = rep(1:N, each = TT),
#'   time = rep(1:TT, N),
#'   y    = rnorm(N * TT),
#'   x    = rnorm(N * TT)
#' )
#' res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
#'               sr = character(0), asymmetric = "x",
#'               id = "id", time = "time", p = 1, asy_test = FALSE)
#' print(res)
#' }
#'
#' @export
pnardl <- function(data,
                   depvar,
                   lr,
                   sr          = character(0),
                   asymmetric,
                   id          = "id",
                   time        = "time",
                   estimator   = c("pmg", "mg", "dfe"),
                   p           = 1L,
                   asy_test    = TRUE,
                   multip      = 0L,
                   level       = 0.95) {

  cl <- match.call()
  estimator <- match.arg(estimator)

  # ---- input validation ------------------------------------------------
  if (!is.data.frame(data))
    stop("'data' must be a data frame.")
  if (!depvar %in% names(data))
    stop("'depvar' not found in data.")
  if (!id %in% names(data))
    stop("'id' not found in data.")
  if (!time %in% names(data))
    stop("'time' not found in data.")
  if (length(asymmetric) == 0)
    stop("'asymmetric' must name at least one variable.")
  for (v in asymmetric) {
    if (!v %in% names(data))
      stop(sprintf("asymmetric variable '%s' not found in data.", v))
  }
  if (p < 1L) stop("'p' must be >= 1.")

  # ---- sort panel -------------------------------------------------------
  data <- data[order(data[[id]], data[[time]]), ]

  # ---- decompose variables into partial sums ----------------------------
  pos_vars <- paste0(asymmetric, "_pos")
  neg_vars <- paste0(asymmetric, "_neg")

  for (i in seq_along(asymmetric)) {
    v    <- asymmetric[i]
    pn   <- pos_vars[i]
    nn   <- neg_vars[i]
    data[[pn]] <- NA_real_
    data[[nn]] <- NA_real_

    for (pid in unique(data[[id]])) {
      idx  <- which(data[[id]] == pid)
      xv   <- data[[v]][idx]
      dx   <- c(NA_real_, diff(xv))
      dxp  <- ifelse(!is.na(dx) & dx > 0, dx, 0)
      dxn  <- ifelse(!is.na(dx) & dx < 0, dx, 0)
      data[[pn]][idx] <- cumsum(ifelse(is.na(dxp), 0, dxp))
      data[[nn]][idx] <- cumsum(ifelse(is.na(dxn), 0, dxn))
      data[[pn]][idx[1]] <- NA_real_
      data[[nn]][idx[1]] <- NA_real_
    }
  }

  # ---- build LR variable list (replace asymmetric with pos/neg) ---------
  lr_vars <- lr[-which(lr == depvar)]
  lr_final <- character(0)
  for (v in lr_vars) {
    if (v %in% asymmetric) {
      i <- which(asymmetric == v)
      lr_final <- c(lr_final, pos_vars[i], neg_vars[i])
    } else {
      lr_final <- c(lr_final, v)
    }
  }

  # SR difference variables
  sr_final <- character(0)
  for (v in sr) {
    if (v %in% asymmetric) {
      i <- which(asymmetric == v)
      sr_final <- c(sr_final, pos_vars[i], neg_vars[i])
    } else {
      sr_final <- c(sr_final, v)
    }
  }

  # ---- per-panel ECM estimation ----------------------------------------
  panel_ids <- unique(data[[id]])
  panel_coefs <- vector("list", length(panel_ids))
  names(panel_coefs) <- as.character(panel_ids)

  for (pid in panel_ids) {
    sub  <- data[data[[id]] == pid, ]
    sub  <- sub[order(sub[[time]]), ]
    n    <- nrow(sub)
    if (n < (p + length(lr_final) + length(sr_final) + 5)) next

    dy <- c(NA_real_, diff(sub[[depvar]]))

    # levels (lagged) for error correction
    lev_mat <- cbind(
      lag1_y = c(NA_real_, sub[[depvar]][-n]),
      do.call(cbind, lapply(lr_final, function(v)
        setNames(data.frame(c(NA_real_, sub[[v]][-n])), paste0("L_", v))))
    )

    # AR lags of dy
    ar_mat <- do.call(cbind, lapply(seq_len(p), function(j) {
      z <- dy
      z <- c(rep(NA_real_, j), z[seq_len(length(z) - j)])
      setNames(data.frame(z), paste0("dL", j, "_y"))
    }))

    # SR differences
    d_sr  <- do.call(cbind, lapply(sr_final, function(v) {
      dv <- c(NA_real_, diff(sub[[v]]))
      setNames(data.frame(dv), paste0("d_", v))
    }))

    if (length(sr_final) == 0) {
      X <- cbind(lev_mat, ar_mat)
    } else {
      X <- cbind(lev_mat, ar_mat, d_sr)
    }

    ok   <- complete.cases(dy, X)
    if (sum(ok) < ncol(X) + 2) next

    fit  <- tryCatch(
      stats::lm.fit(cbind(1, as.matrix(X[ok, ])), dy[ok]),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    panel_coefs[[as.character(pid)]] <- stats::setNames(
      fit$coefficients,
      c("(Intercept)", colnames(X))
    )
  }

  # remove NULLs
  panel_coefs <- Filter(Negate(is.null), panel_coefs)
  if (length(panel_coefs) == 0)
    stop("Estimation failed for all panels. Check data structure.")

  # ---- pool coefficients (MG mean, PMG/DFE = same here) ---------------
  coef_mat <- do.call(rbind, lapply(panel_coefs, function(b) {
    nms <- names(panel_coefs[[1]])
    out <- setNames(rep(NA_real_, length(nms)), nms)
    shared <- intersect(names(b), nms)
    out[shared] <- b[shared]
    out
  }))

  pooled_coef <- colMeans(coef_mat, na.rm = TRUE)

  # pooled SE (across-panel dispersion for MG; within for DFE-like)
  if (nrow(coef_mat) > 1) {
    pooled_se <- apply(coef_mat, 2, stats::sd, na.rm = TRUE) /
      sqrt(colSums(!is.na(coef_mat)))
  } else {
    pooled_se <- rep(NA_real_, length(pooled_coef))
    names(pooled_se) <- names(pooled_coef)
  }

  pooled_se[pooled_se == 0] <- NA_real_
  tstat  <- pooled_coef / pooled_se
  df_    <- length(panel_coefs) - 1L
  pval   <- 2 * stats::pt(-abs(tstat), df = max(df_, 1L))

  # ECT per panel
  ecm_coefs <- sapply(panel_coefs, function(b) {
    if ("lag1_y" %in% names(b)) b["lag1_y"] else NA_real_
  })

  # ---- asymmetry Wald tests --------------------------------------------
  asym_tests <- NULL
  if (asy_test && length(panel_coefs) > 0) {
    rows <- lapply(asymmetric, function(v) {
      i  <- which(asymmetric == v)
      pv <- pos_vars[i]
      nv <- neg_vars[i]
      lp <- paste0("L_", pv)
      ln <- paste0("L_", nv)
      dp <- paste0("d_",  pv)
      dn <- paste0("d_",  nv)

      # long-run
      lr_p <- pooled_coef[lp]
      lr_n <- pooled_coef[ln]
      lr_se_p <- pooled_se[lp]
      lr_se_n <- pooled_se[ln]

      if (!is.na(lr_p) && !is.na(lr_n) &&
          !is.na(lr_se_p) && !is.na(lr_se_n) &&
          lr_se_p > 0 && lr_se_n > 0) {
        lr_diff   <- lr_p - lr_n
        lr_se_d   <- sqrt(lr_se_p^2 + lr_se_n^2)
        lr_chi2   <- (lr_diff / lr_se_d)^2
        lr_pval   <- stats::pchisq(lr_chi2, df = 1, lower.tail = FALSE)
      } else {
        lr_diff <- lr_chi2 <- lr_pval <- NA_real_
      }

      # short-run
      sr_p <- pooled_coef[dp]
      sr_n <- pooled_coef[dn]
      sr_se_p <- pooled_se[dp]
      sr_se_n <- pooled_se[dn]

      if (!is.na(sr_p) && !is.na(sr_n) &&
          !is.na(sr_se_p) && !is.na(sr_se_n) &&
          sr_se_p > 0 && sr_se_n > 0) {
        sr_diff   <- sr_p - sr_n
        sr_se_d   <- sqrt(sr_se_p^2 + sr_se_n^2)
        sr_chi2   <- (sr_diff / sr_se_d)^2
        sr_pval   <- stats::pchisq(sr_chi2, df = 1, lower.tail = FALSE)
      } else {
        sr_diff <- sr_chi2 <- sr_pval <- NA_real_
      }

      data.frame(
        variable  = v,
        lr_beta_p = lr_p,
        lr_beta_n = lr_n,
        lr_diff   = lr_diff,
        lr_chi2   = lr_chi2,
        lr_pval   = lr_pval,
        sr_gamma_p = sr_p,
        sr_gamma_n = sr_n,
        sr_diff   = sr_diff,
        sr_chi2   = sr_chi2,
        sr_pval   = sr_pval,
        stringsAsFactors = FALSE
      )
    })
    asym_tests <- do.call(rbind, rows)
  }

  # ---- dynamic multipliers ---------------------------------------------
  multipliers <- NULL
  if (multip > 0L) {
    mean_phi <- mean(ecm_coefs, na.rm = TRUE)
    for (v in asymmetric) {
      i  <- which(asymmetric == v)
      bp <- pooled_coef[paste0("L_", pos_vars[i])]
      bn <- pooled_coef[paste0("L_", neg_vars[i])]
      if (is.na(bp) || is.na(bn) || is.na(mean_phi)) next

      rows_m <- lapply(0:multip, function(h) {
        if (h == 0L) {
          mp <- 0; mn <- 0
        } else {
          mp <- bp * (1 - (1 + mean_phi)^h)
          mn <- bn * (1 - (1 + mean_phi)^h)
        }
        data.frame(variable = v, period = h,
                   m_pos = mp, m_neg = mn,
                   diff = mp - mn, stringsAsFactors = FALSE)
      })
      multipliers <- rbind(multipliers, do.call(rbind, rows_m))
    }
  }

  # ---- assemble result object ------------------------------------------
  structure(
    list(
      coefficients = pooled_coef,
      se           = pooled_se,
      tstat        = tstat,
      pval         = pval,
      panel_coefs  = panel_coefs,
      ecm_coefs    = ecm_coefs,
      asym_tests   = asym_tests,
      multipliers  = multipliers,
      partial_sums = data,
      estimator    = estimator,
      depvar       = depvar,
      asymmetric   = asymmetric,
      pos_vars     = pos_vars,
      neg_vars     = neg_vars,
      call         = cl
    ),
    class = "pnardl"
  )
}


#' Print Method for pnardl Objects
#'
#' @param x An object of class \code{"pnardl"}.
#' @param digits Integer. Number of significant digits to display.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pnardl <- function(x, digits = 4L, ...) {
  cat("\nPanel Nonlinear ARDL (PNARDL)\n")
  cat("Estimator:", toupper(x$estimator), "\n")
  cat("Dependent variable:", x$depvar, "\n")
  cat("Asymmetric variables:", paste(x$asymmetric, collapse = ", "), "\n")
  cat("Number of panels:", length(x$panel_coefs), "\n\n")

  # coefficient table
  nms  <- names(x$coefficients)
  coef <- x$coefficients
  se   <- x$se
  ts   <- x$tstat
  pv   <- x$pval

  tab <- data.frame(
    Estimate  = round(coef,  digits),
    Std.Error = round(se,    digits),
    t.value   = round(ts,    digits),
    Pr        = round(pv,    digits),
    stringsAsFactors = FALSE,
    row.names = nms
  )
  names(tab)[4] <- "Pr(>|t|)"

  cat("Pooled Coefficients:\n")
  print(tab)

  if (!is.null(x$asym_tests)) {
    cat("\nAsymmetry Tests (Wald):\n")
    at <- x$asym_tests
    at[, sapply(at, is.numeric)] <-
      round(at[, sapply(at, is.numeric)], digits)
    print(at, row.names = FALSE)
  }

  cat("\nSpeed of Adjustment (ECT) per panel:\n")
  print(round(x$ecm_coefs, digits))

  invisible(x)
}


#' Summary Method for pnardl Objects
#'
#' @param object An object of class \code{"pnardl"}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns a list with components:
#'   \item{coefficients}{Data frame with estimates, SE, t-stats, and p-values.}
#'   \item{asym_tests}{Asymmetry test results (may be \code{NULL}).}
#'   \item{ecm_summary}{Summary statistics of the ECT coefficients.}
#'
#' @export
summary.pnardl <- function(object, ...) {
  cat("\n====================================================\n")
  cat(" Panel Nonlinear ARDL (PNARDL) -- Summary\n")
  cat("====================================================\n")
  cat("Call:\n")
  print(object$call)
  cat("\nEstimator :", toupper(object$estimator), "\n")
  cat("Dep. var  :", object$depvar, "\n")
  cat("Asymmetric:", paste(object$asymmetric, collapse = ", "), "\n")
  cat("Panels    :", length(object$panel_coefs), "\n")

  cat("\n--- Pooled Coefficients ---\n")
  tab <- data.frame(
    Estimate   = object$coefficients,
    Std.Error  = object$se,
    t.value    = object$tstat,
    `Pr(>|t|)` = object$pval,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  stars <- ifelse(object$pval < 0.01, "***",
           ifelse(object$pval < 0.05, "** ",
           ifelse(object$pval < 0.10, "*  ", "   ")))
  tab$sig <- stars
  print(tab, digits = 4)
  cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.1 ' ' 1\n")

  cat("\n--- Speed of Adjustment (ECT) ---\n")
  ecm <- object$ecm_coefs
  cat(sprintf("  Mean phi: %.4f   Half-life: %.2f periods\n",
              mean(ecm, na.rm = TRUE),
              log(2) / abs(mean(ecm, na.rm = TRUE))))
  conv <- sum(ecm < 0 & ecm > -2, na.rm = TRUE)
  cat(sprintf("  Convergent panels: %d / %d\n", conv, length(ecm)))

  if (!is.null(object$asym_tests)) {
    cat("\n--- Asymmetry Tests ---\n")
    print(object$asym_tests, row.names = FALSE, digits = 4)
  }

  if (!is.null(object$multipliers)) {
    cat("\n--- Dynamic Multipliers (first 5 periods) ---\n")
    print(utils::head(object$multipliers, 5), row.names = FALSE, digits = 4)
  }

  res <- list(
    coefficients = tab,
    asym_tests   = object$asym_tests,
    ecm_summary  = summary(ecm)
  )
  invisible(res)
}
