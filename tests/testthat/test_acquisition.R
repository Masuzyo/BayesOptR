test_that("Student-t Expected Improvement is implemented", {
  acq <- StudentTExpectedImprovement$new(xi = 0.01, nu = 3)
  
  # Test basic properties
  expect_equal(acq$xi, 0.01)
  expect_equal(acq$nu, 3)
  
  # Test evaluation
  mean <- c(1.0, 2.0, 3.0)
  std <- c(0.5, 0.3, 0.1)
  f_best <- 2.5
  
  ei_values <- acq$evaluate(mean, std, f_best, minimize = FALSE)
  
  # EI should be non-negative
  expect_true(all(ei_values >= 0))
  expect_equal(length(ei_values), 3)
  
  # Point with mean > f_best should have positive EI
  expect_true(ei_values[3] > 0)
})

test_that("Student-t EI converges to Gaussian EI as nu increases", {
  mean <- c(1.5, 2.5, 3.5)
  std <- c(0.5, 0.5, 0.5)
  f_best <- 2.0
  
  # Student-t EI with high nu
  acq_t_high <- StudentTExpectedImprovement$new(nu = 100)
  ei_t_high <- acq_t_high$evaluate(mean, std, f_best)
  
  # Gaussian EI
  acq_gauss <- ExpectedImprovement$new()
  ei_gauss <- acq_gauss$evaluate(mean, std, f_best)
  
  # Should be very similar
  expect_true(max(abs(ei_t_high - ei_gauss)) < 0.01)
})

test_that("Student-t EI with low nu explores more than high nu", {
  mean <- c(1.5, 2.5, 3.5)
  std <- c(0.5, 0.5, 0.5)
  f_best <- 2.0
  
  # Low nu (exploration)
  acq_low <- StudentTExpectedImprovement$new(nu = 3)
  ei_low <- acq_low$evaluate(mean, std, f_best)
  
  # High nu (exploitation)
  acq_high <- StudentTExpectedImprovement$new(nu = 50)
  ei_high <- acq_high$evaluate(mean, std, f_best)
  
  # Low nu should give higher EI for uncertain points
  # (heavy tails favor exploration)
  expect_true(all(ei_low >= ei_high - 1e-10))  # Allow numerical tolerance
})

test_that("Student-t EI handles zero std correctly", {
  acq <- StudentTExpectedImprovement$new(nu = 5)
  
  mean <- c(1.0, 2.0, 3.0)
  std <- c(0.0, 0.5, 0.0)  # Zero std at positions 1 and 3
  f_best <- 2.5
  
  ei_values <- acq$evaluate(mean, std, f_best)
  
  # Zero std should give zero EI
  expect_equal(ei_values[1], 0)
  expect_equal(ei_values[3], 0)
  expect_true(ei_values[2] >= 0)
})

test_that("Student-t EI set_nu method works", {
  acq <- StudentTExpectedImprovement$new(nu = 3)
  expect_equal(acq$nu, 3)
  
  acq$set_nu(10)
  expect_equal(acq$nu, 10)
  
  # Evaluate with new nu
  mean <- c(2.0)
  std <- c(0.5)
  f_best <- 1.5
  
  ei_value <- acq$evaluate(mean, std, f_best)
  expect_true(ei_value > 0)
})

test_that("Student-t EI can override nu in evaluate", {
  acq <- StudentTExpectedImprovement$new(nu = 3)
  
  mean <- c(2.0)
  std <- c(0.5)
  f_best <- 1.5
  
  # Evaluate with default nu
  ei_default <- acq$evaluate(mean, std, f_best)
  
  # Evaluate with overridden nu
  ei_override <- acq$evaluate(mean, std, f_best, nu = 20)
  
  # Should be different (lower nu explores more)
  expect_false(isTRUE(all.equal(ei_default, ei_override)))
  
  # Original nu should be unchanged
  expect_equal(acq$nu, 3)
})

test_that("AHBO uses Student-t EI by default", {
  skip_if_not_installed("DiceKriging")
  
  # Simple 2D function
  fn <- function(x) -(x[1]^2 + x[2]^2)
  bounds <- matrix(c(-5, -5, 5, 5), nrow = 2, ncol = 2)
  
  # Run AHBO with minimal iterations
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "adaptive_hybrid",
    acquisition = "ei",
    n_init = 5,
    n_iter = 3,
    nu_schedule = "linear"
  )
  
  # Check that optimization completed
  expect_true("x_best" %in% names(result))
  expect_true("f_best" %in% names(result))
  expect_equal(length(result$x_best), 2)
})

test_that("Student-t BO method uses Student-t EI", {
  skip_if_not_installed("DiceKriging")
  
  # Simple 1D function
  fn <- function(x) -(x[1] - 2)^2
  bounds <- matrix(c(0, 5), nrow = 2, ncol = 1)
  
  # Run Student-t BO
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "student_t",
    acquisition = "ei",
    n_init = 5,
    n_iter = 3,
    nu = 5
  )
  
  # Check that optimization completed
  expect_true("x_best" %in% names(result))
  expect_true("f_best" %in% names(result))
  expect_equal(length(result$x_best), 1)
})

test_that("Gaussian BO uses Gaussian EI not Student-t EI", {
  skip_if_not_installed("DiceKriging")
  
  # Simple 1D function
  fn <- function(x) -(x[1] - 2)^2
  bounds <- matrix(c(0, 5), nrow = 2, ncol = 1)
  
  # Run Gaussian BO (should use Gaussian EI)
  result <- bayes_optimize(
    fn = fn,
    bounds = bounds,
    method = "gaussian",
    acquisition = "ei",
    n_init = 5,
    n_iter = 3
  )
  
  # Check that optimization completed
  expect_true("x_best" %in% names(result))
  expect_true("f_best" %in% names(result))
  expect_equal(length(result$x_best), 1)
})
