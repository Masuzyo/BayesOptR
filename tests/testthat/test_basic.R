library(testthat)

test_that("GaussianBO initialization works", {
  bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  model <- GaussianBO$new(bounds)
  
  expect_equal(nrow(model$bounds), 2)
  expect_equal(ncol(model$bounds), 2)
  expect_equal(model$kernel, "matern52")
  expect_equal(model$n_init, 10)  # 5 * n_dim
})

test_that("GaussianBO add_observation works", {
  bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  model <- GaussianBO$new(bounds)
  
  x <- c(0.5, 0.5)
  y <- 1.0
  
  model$add_observation(x, y)
  
  expect_equal(nrow(model$X), 1)
  expect_equal(model$y, 1.0)
  expect_equal(as.vector(model$X[1, ]), x)
})

test_that("GaussianBO prediction works", {
  bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  model <- GaussianBO$new(bounds)
  
  # Add some observations
  model$add_observation(c(0.2, 0.3), 1.0)
  model$add_observation(c(0.8, 0.7), 0.5)
  model$add_observation(c(0.5, 0.5), 0.8)
  
  # Predict
  X_new <- matrix(c(0.4, 0.4, 0.6, 0.6), nrow = 2, ncol = 2, byrow = TRUE)
  pred <- model$predict(X_new)
  
  expect_true("mean" %in% names(pred))
  expect_true("variance" %in% names(pred))
  expect_equal(length(pred$mean), 2)
  expect_equal(length(pred$variance), 2)
})


test_that("StudentTBO inherits from GaussianBO", {
  bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  model <- StudentTBO$new(bounds, nu = 3)
  
  expect_equal(model$nu, 3)
  expect_true(inherits(model, "StudentTBO"))
  expect_true(inherits(model, "GaussianBO"))
})


test_that("AdaptiveHybridBO nu adaptation works", {
  bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  model <- AdaptiveHybridBO$new(bounds, nu_min = 3, nu_max = 50, 
                                nu_schedule = "linear", max_iterations = 100)
  
  expect_equal(model$nu_min, 3)
  expect_equal(model$nu_max, 50)
  expect_equal(model$current_nu, 3)
  
  # Test nu progression
  model$update_nu(50)
  expect_true(model$current_nu > 3)
  expect_true(model$current_nu < 50)
  
  model$update_nu(100)
  expect_equal(model$current_nu, 50)
})


test_that("Acquisition functions work", {
  mean <- c(1.0, 2.0, 3.0)
  std <- c(0.1, 0.2, 0.3)
  f_best <- 2.5
  
  # Expected Improvement
  ei <- ExpectedImprovement$new()
  ei_vals <- ei$evaluate(mean, std, f_best)
  expect_equal(length(ei_vals), 3)
  expect_true(all(ei_vals >= 0))
  
  # UCB
  ucb <- UpperConfidenceBound$new(kappa = 2.0)
  ucb_vals <- ucb$evaluate(mean, std)
  expect_equal(length(ucb_vals), 3)
  expect_true(all(ucb_vals > mean))  # UCB should be higher than mean
  
  # PI
  pi <- ProbabilityOfImprovement$new()
  pi_vals <- pi$evaluate(mean, std, f_best)
  expect_equal(length(pi_vals), 3)
  expect_true(all(pi_vals >= 0 & pi_vals <= 1))
})


test_that("Blocking strategies initialize correctly", {
  # Dimension-wise hard
  blocking <- DimensionWiseHardBlocking$new(block_size = 5, n_recent = 3)
  expect_equal(blocking$block_size, 5)
  expect_equal(blocking$n_recent, 3)
  
  # Distance-based soft
  blocking <- DistanceSoftBlocking$new(radius = 0.1)
  expect_equal(blocking$radius, 0.1)
  
  # Local penalization
  blocking <- LocalPenalization$new(lipschitz_constant = 10)
  expect_equal(blocking$lipschitz_constant, 10)
  
  # Constant liar
  blocking <- ConstantLiar$new(strategy = "mean")
  expect_equal(blocking$strategy, "mean")
})
