
# scmBayesPost

<!-- badges: start -->
[![R-CMD-check](https://github.com/andrebonfrer/scmBayesPost/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andrebonfrer/scmBayesPost/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/andrebonfrer/scmBayesPost/graph/badge.svg)](https://app.codecov.io/gh/andrebonfrer/scmBayesPost)
<!-- badges: end -->

The goal of scmBayesPost is to separate weight 
estimation from post-treatment inference in synthetic control designs.

The package provides a flexible post-estimation framework that:

Accepts synthetic control weights from any source (SCM, ASCM, multi-outcome 
synth, or custom estimators),

Constructs treated-specific pseudo-panels,

Estimates unit-level treatment effects,

Models heterogeneity via moderators,

Supports instrumental variables and control-function approaches for endogenous moderators,

Allows both binary and continuous treatment intensities.

In contrast to implementations tied to a specific weight estimator, 
scmBayesPost treats the weight matrix as a first-class object and focuses on 
rigorous second-stage inference.

## Installation

You can install the development version of scmBayesPost from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("andrebonfrer/scmBayesPost")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(scmBayesPost)
## basic example code
```

