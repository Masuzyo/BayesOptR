library(testthat)

test_that("Simple optimization works end-to-end", {
  skip_on_cran()
  
  # Simple quadratic function
  fn <- function(x) -(x[1]^2 + x[2]^2)  # Maximum at (0, 0)
  bounds <- matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2)
  
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "gaussian",
    n_iter = 10,
    verbose = FALSE
  )
  
  expect_true(inherits(result, "bayesopt"))
  expect_equal(result$n_iter, 10)
  expect_true(result$best_y <= 0)  # Should find near (0,0)
  expect_true(abs(result$best_y) < 0.5)  # Reasonable convergence
})


test_that("Batch optimization works", {
  skip_on_cran()
  
  fn <- function(x) -(x[1]^2 + x[2]^2)
  bounds <- matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2)
  
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "student_t",
    batch_size = 4,
    n_iter = 5,
    verbose = FALSE
  )
  
  expect_equal(nrow(result$X), 5 * 4 + 10)  # 5 iters × 4 batch + 10 init
})


test_that("Blocking strategies work in optimization", {
  skip_on_cran()
  
  fn <- function(x) -(x[1]^2 + x[2]^2)
  bounds <- matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2)
  
  # Test each blocking strategy
  strategies <- c("dimwise_hard", "dimwise_soft", "distance_hard")
  
  for (strategy in strategies) {
    result <- bayes_optimize(
      fn = fn,
      bounds = bounds,
      batch_size = 2,
      blocking_strategy = strategy,
      n_iter = 5,
      verbose = FALSE
    )
    
    expect_true(inherits(result, "bayesopt"))
    expect_equal(result$blocking_strategy, strategy)
  }
})


test_that("Adaptive hybrid converges", {
  skip_on_cran()
  
  fn <- function(x) -(x[1]^2 + x[2]^2)
  bounds <- matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2)
  
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "adaptive_hybrid",
    nu_schedule = "linear",
    n_iter = 20,
    verbose = FALSE
  )
  
  # Check that convergence is improving
  expect_true(tail(result$convergence, 1) > result$convergence[1])
  
  # Check nu adaptation
  diagnostics <- result$model$get_diagnostics()
  expect_true(diagnostics$current_nu > 3)
})


test_that("Minimization works", {
  skip_on_cran()
  
  # Minimize quadratic (minimum at origin)
  fn <- function(x) x[1]^2 + x[2]^2
  bounds <- matrix(c(-1, -1, 1, 1), nrow = 2, ncol = 2)
  
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    minimize = TRUE,
    n_iter = 15,
    verbose = FALSE
  )
  
  expect_true(result$best_y >= 0)
  expect_true(result$best_y < 0.2)  # Should find near minimum
})
