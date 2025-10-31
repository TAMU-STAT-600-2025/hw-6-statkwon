# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

testthat::test_that("Test for LassoInC.cpp", {
  # Generate mocking data
  n1 <- 20
  p1 <- 5
  X1 <- matrix(rnorm(n1 * p1), n1, p1)
  beta1 <- rnorm(p1)
  epsilon1 <- rnorm(n1)
  Y1 <- X1 %*% beta1 + epsilon1
  out1 <- standardizeXY(X1, Y1)
  Xtilde1 <- out1$Xtilde
  Ytilde1 <- out1$Ytilde
  lambda1 <- 0.1
  
  n2 <- 100
  p2 <- 10
  X2 <- matrix(rnorm(n2 * p2), n2, p2)
  beta2 <- rnorm(p2)
  epsilon2 <- rnorm(n2)
  Y2 <- X2 %*% beta2 + epsilon2
  out2 <- standardizeXY(X2, Y2)
  Xtilde2 <- out2$Xtilde
  Ytilde2 <- out2$Ytilde
  lambda2 <- 0.3
  
  beta_start <- rep(0, p1)
  
  # Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
  #################################################
  testthat::expect_equal(soft(3.5, 1), soft_c(3.5, 1))
  testthat::expect_equal(soft(-5, 1), soft_c(-5, 1))
  
  # Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
  #################################################
  testthat::expect_equal(
    lasso(Xtilde1, Ytilde1, beta1, lambda1),
    lasso_c(Xtilde1, Ytilde1, beta1, lambda1),
    check.attributes = FALSE
  )
  testthat::expect_equal(
    lasso(Xtilde2, Ytilde2, beta2, lambda2),
    lasso_c(Xtilde2, Ytilde2, beta2, lambda2),
    check.attributes = FALSE
  )
  
  # Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
  #################################################
  testthat::expect_equal(
    fitLASSOstandardized(Xtilde1, Ytilde1, lambda1, beta_start)$beta,
    fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda1, beta_start),
    check.attributes = FALSE
  )
  
  # Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
  #################################################
  out <- fitLASSOstandardized_seq(Xtilde1, Ytilde1)
  lambda_seq <- out$lambda_seq
  testthat::expect_equal(
    out$beta,
    fitLASSOstandardized_seq_c(Xtilde1, Ytilde1, lambda_seq),
    check.attributes = FALSE
  )
  
  # Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
  ######################################################################
  res <- microbenchmark::microbenchmark(
    fitLASSOstandardized(Xtilde1, Ytilde1, lambda1, beta_start),
    fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda1, beta_start),
    times = 10
  )
  print(res)
  
  # Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
  ######################################################################
  res <- microbenchmark::microbenchmark(
    out <- fitLASSOstandardized_seq(Xtilde1, Ytilde1, lambda_seq),
    fitLASSOstandardized_seq_c(Xtilde1, Ytilde1, out$lambda_seq),
    times = 10
  )
  print(res)
})


# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark::microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
