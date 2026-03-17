# pnardl

**Panel Nonlinear ARDL Estimation**

The `pnardl` R package implements the Panel Nonlinear ARDL (PNARDL) model following Shin, Yu and Greenwood-Nimmo (2014). It decomposes regressors into positive and negative partial sums to capture asymmetric long-run and short-run effects in panel data settings.

## Features

- Positive and negative partial-sum decomposition
- PMG, MG, and DFE pooled estimators
- Wald tests for long-run and short-run asymmetry
- Cumulative dynamic multipliers
- Impulse response functions

## Installation

```r
# Install from CRAN (once available):
install.packages("pnardl")
```

## Usage

```r
library(pnardl)

# Simulate panel data
set.seed(42)
N <- 5; TT <- 30
df <- data.frame(
  id   = rep(1:N, each = TT),
  time = rep(1:TT, N),
  y    = rnorm(N * TT),
  x    = rnorm(N * TT)
)

res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
              sr = character(0), asymmetric = "x",
              id = "id", time = "time")
summary(res)
```

## Reference

Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric cointegration and dynamic multipliers in a nonlinear ARDL framework. In R. C. Sickles & W. C. Horrace (Eds.), *Festschrift in Honor of Peter Schmidt* (pp. 281–314). Springer. https://doi.org/10.1007/978-1-4899-8008-3_9

## License

GPL-3
