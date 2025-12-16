#!/usr/bin/env Rscript

# BayesOptR Package Testing Script
# Comprehensive tests for all package components

cat("================================================================================\n")
cat("BayesOptR PACKAGE TESTING\n")
cat("================================================================================\n\n")

# Load required libraries
library(R6)

# Source all package files
cat("Loading package files...\n")
source("R/acquisition_functions.R")
source("R/gaussian_bo.R")
source("R/student_t_bo.R")
source("R/adaptive_hybrid_bo.R")
source("R/dimension_wise_blocking.R")
source("R/distance_blocking.R")
source("R/other_blocking_strategies.R")
source("R/batch_acquisition.R")
source("R/main_interface.R")
source("R/visualization.R")
cat("[OK] All files loaded\n\n")

# Define test functions
branin <- function(x) {
  x1 <- x[1] * 15 - 5
  x2 <- x[2] * 15
  -(1 * (x2 - 5.1/(4*pi^2) * x1^2 + 5/pi * x1 - 6)^2 + 
    10 * (1 - 1/(8*pi)) * cos(x1) + 10)
}

sphere <- function(x) {
  -(sum((x * 10 - 5)^2))
}

bounds_2d <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)

# Test counter
test_num <- 0
tests_passed <- 0
tests_failed <- 0

run_test <- function(test_name, test_func) {
  test_num <<- test_num + 1
  cat(sprintf("\n[Test %d] %s\n", test_num, test_name))
  cat(rep("-", 80), "\n", sep="")
  
  result <- tryCatch({
    test_func()
    cat("[PASS] Test passed\n")
    tests_passed <<- tests_passed + 1
    TRUE
  }, error = function(e) {
    cat("[FAIL] Test failed:", e$message, "\n")
    tests_failed <<- tests_failed + 1
    FALSE
  })
  
  return(result)
}

# ============================================================================
# TEST 1: Acquisition Functions
# ============================================================================

run_test("Gaussian Expected Improvement", function() {
  acq <- ExpectedImprovement$new(xi = 0.01)
  mean_vals <- c(1.5, 2.5, 3.5)
  std_vals <- c(0.5, 0.5, 0.5)
  f_best <- 2.0
  
  ei_values <- acq$evaluate(mean_vals, std_vals, f_best)
  
  cat("  EI values:", sprintf("%.4f ", ei_values), "\n")
  stopifnot(all(ei_values >= 0))
  stopifnot(length(ei_values) == 3)
})

run_test("Student-t Expected Improvement", function() {
  acq <- StudentTExpectedImprovement$new(xi = 0.01, nu = 5)
  mean_vals <- c(1.5, 2.5, 3.5)
  std_vals <- c(0.5, 0.5, 0.5)
  f_best <- 2.0
  
  ei_values <- acq$evaluate(mean_vals, std_vals, f_best)
  
  cat("  EI values (ν=5):", sprintf("%.4f ", ei_values), "\n")
  stopifnot(all(ei_values >= 0))
  
  # Test nu override
  ei_nu10 <- acq$evaluate(mean_vals, std_vals, f_best, nu = 10)
  cat("  EI values (ν=10):", sprintf("%.4f ", ei_nu10), "\n")
})

run_test("Student-t EI converges to Gaussian EI", function() {
  acq_t <- StudentTExpectedImprovement$new(nu = 100)
  acq_g <- ExpectedImprovement$new()
  
  mean_vals <- c(1.5, 2.5, 3.5)
  std_vals <- c(0.5, 0.5, 0.5)
  f_best <- 2.0
  
  ei_t <- acq_t$evaluate(mean_vals, std_vals, f_best)
  ei_g <- acq_g$evaluate(mean_vals, std_vals, f_best)
  
  diff <- max(abs(ei_t - ei_g))
  cat(sprintf("  Max difference: %.6f\n", diff))
  stopifnot(diff < 0.01)
})

run_test("Upper Confidence Bound", function() {
  acq <- UpperConfidenceBound$new(kappa = 2.0)
  mean_vals <- c(1.5, 2.5, 3.5)
  std_vals <- c(0.5, 0.5, 0.5)
  
  ucb_values <- acq$evaluate(mean_vals, std_vals, n_dim = 2)
  cat("  UCB values:", sprintf("%.4f ", ucb_values), "\n")
  stopifnot(all(!is.na(ucb_values)))
})

# ============================================================================
# TEST 2: Bayesian Optimization Methods
# ============================================================================

run_test("Gaussian BO", function() {
  set.seed(123)
  model <- GaussianBO$new(bounds_2d, n_init = 5)
  X_init <- model$generate_initial_samples()
  y_init <- apply(X_init, 1, branin)
  
  for (i in 1:nrow(X_init)) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  test_point <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
  pred <- model$predict(test_point)
  
  cat(sprintf("  Prediction: mean=%.4f, std=%.4f\n", pred$mean, pred$sd))
  stopifnot(!is.null(pred$mean))
  stopifnot(pred$sd > 0)
})

run_test("Student-t BO", function() {
  set.seed(123)
  model <- StudentTBO$new(bounds_2d, nu = 5, n_init = 5)
  X_init <- model$generate_initial_samples()
  y_init <- apply(X_init, 1, branin)
  
  for (i in 1:nrow(X_init)) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  test_point <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
  pred <- model$predict(test_point)
  
  cat(sprintf("  Prediction: mean=%.4f, std=%.4f, nu=%.1f\n", 
              pred$mean, pred$sd, pred$nu))
  stopifnot(!is.null(pred$nu))
})

run_test("Adaptive Hybrid BO - Nu Adaptation", function() {
  set.seed(123)
  model <- AdaptiveHybridBO$new(
    bounds_2d, 
    nu_min = 3, 
    nu_max = 50,
    nu_schedule = "linear",
    max_iterations = 20,
    n_init = 5
  )
  
  X_init <- model$generate_initial_samples()
  y_init <- apply(X_init, 1, branin)
  
  for (i in 1:nrow(X_init)) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  # Test nu adaptation
  nu_vals <- numeric(10)
  for (i in 1:10) {
    model$iteration <- i
    model$update_nu()
    nu_vals[i] <- model$nu
  }
  
  cat("  Nu progression:", sprintf("%.1f ", nu_vals), "\n")
  stopifnot(nu_vals[1] < nu_vals[10])
  stopifnot(nu_vals[1] >= 3 && nu_vals[10] <= 50)
})

# ============================================================================
# TEST 3: Blocking Strategies
# ============================================================================

run_test("Dimension-Wise Hard Blocking", function() {
  blocking <- DimensionWiseHardBlocking$new(block_size = 5, n_recent = 3)
  
  ref_points <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
  blocking$set_pending(ref_points)
  
  candidates <- matrix(runif(20), nrow = 10, ncol = 2)
  acq_values <- runif(10, 0.5, 1.5)
  
  blocked <- blocking$apply_blocking(acq_values, candidates, bounds_2d)
  
  n_blocked <- sum(is.infinite(blocked))
  cat(sprintf("  Blocked points: %d/10\n", n_blocked))
  stopifnot(n_blocked > 0)
})

run_test("Distance-Based Soft Blocking", function() {
  blocking <- DistanceSoftBlocking$new(
    radius = 0.2, 
    penalty_strength = 3.0,
    kernel = "exponential"
  )
  
  ref_points <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
  blocking$set_pending(ref_points)
  
  candidates <- matrix(runif(20), nrow = 10, ncol = 2)
  acq_values <- runif(10, 0.5, 1.5)
  
  blocked <- blocking$apply_blocking(acq_values, candidates, bounds_2d)
  
  penalty_factor <- mean(blocked / acq_values)
  cat(sprintf("  Mean penalty factor: %.3f\n", penalty_factor))
  stopifnot(all(blocked <= acq_values))
})

run_test("Targeted Blocking", function() {
  blocking <- TargetedBlocking$new()
  
  target_point <- c(0.7, 0.7)
  blocking$set_target(target_point)
  
  all_points <- matrix(c(0.2, 0.2, 0.3, 0.4, 0.5, 0.5), 
                       nrow = 3, ncol = 2, byrow = TRUE)
  candidates <- matrix(runif(20), nrow = 10, ncol = 2)
  acq_values <- runif(10, 0.5, 1.5)
  
  blocked <- blocking$apply_blocking(acq_values, candidates, all_points, bounds_2d)
  
  penalties <- blocking$calculate_penalization(candidates, all_points, bounds_2d)
  cat(sprintf("  Penalty range: [%.3f, %.3f]\n", min(penalties), max(penalties)))
  stopifnot(all(penalties >= 0))
})

run_test("Local Penalization", function() {
  blocking <- LocalPenalization$new(lipschitz_constant = 10.0, kappa = 2.0)
  
  ref_points <- matrix(c(0.5, 0.5), nrow = 1, ncol = 2)
  blocking$set_pending(ref_points)
  
  candidates <- matrix(runif(20), nrow = 10, ncol = 2)
  acq_values <- runif(10, 0.5, 1.5)
  mean_pred <- runif(10)
  std_pred <- runif(10, 0.1, 0.5)
  
  blocked <- blocking$apply_blocking(acq_values, candidates, mean_pred, 
                                     std_pred, bounds_2d)
  
  cat(sprintf("  Mean penalty: %.3f\n", mean(blocked / acq_values)))
  stopifnot(all(blocked <= acq_values))
})

# ============================================================================
# TEST 4: Batch Acquisition
# ============================================================================

run_test("Batch Acquisition without Blocking", function() {
  set.seed(456)
  
  model <- GaussianBO$new(bounds_2d, n_init = 8)
  X_init <- model$generate_initial_samples()
  y_init <- apply(X_init, 1, sphere)
  
  for (i in 1:nrow(X_init)) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  acq_fn <- ExpectedImprovement$new()
  batch_acq <- BatchAcquisition$new(model, acq_fn, NULL, batch_size = 4)
  
  batch <- batch_acq$select_batch()
  
  cat("  Batch points:\n")
  for (i in 1:nrow(batch)) {
    cat(sprintf("    Point %d: (%.3f, %.3f)\n", i, batch[i, 1], batch[i, 2]))
  }
  
  stopifnot(nrow(batch) == 4)
})

run_test("Batch Acquisition with Blocking", function() {
  set.seed(456)
  
  model <- GaussianBO$new(bounds_2d, n_init = 8)
  X_init <- model$generate_initial_samples()
  y_init <- apply(X_init, 1, sphere)
  
  for (i in 1:nrow(X_init)) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  acq_fn <- ExpectedImprovement$new()
  blocking <- DistanceSoftBlocking$new(radius = 0.15, penalty_strength = 5.0)
  batch_acq <- BatchAcquisition$new(model, acq_fn, blocking, batch_size = 4)
  
  batch <- batch_acq$select_batch()
  
  min_dist <- min(dist(batch))
  cat(sprintf("  Minimum inter-point distance: %.3f\n", min_dist))
  
  stopifnot(nrow(batch) == 4)
  stopifnot(min_dist > 0)
})

# ============================================================================
# TEST 5: Full Optimization Runs
# ============================================================================

run_test("Sequential Optimization (Gaussian BO)", function() {
  set.seed(789)
  
  result <- bayes_optimize(
    fn = branin,
    bounds = bounds_2d,
    method = "gaussian",
    acquisition = "ei",
    n_iter = 15,
    n_init = 8,
    batch_size = 1,
    blocking_strategy = "none",
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f (optimal: -0.398)\n", result$best_y))
  cat(sprintf("  Best point: (%.3f, %.3f)\n", result$best_x[1], result$best_x[2]))
  
  stopifnot(!is.null(result$best_y))
  stopifnot(result$best_y > -0.5)
})

run_test("Sequential Optimization (Student-t BO)", function() {
  set.seed(790)
  
  result <- bayes_optimize(
    fn = branin,
    bounds = bounds_2d,
    method = "student_t",
    acquisition = "ei",
    n_iter = 15,
    n_init = 8,
    nu = 5,
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  stopifnot(!is.null(result$best_y))
})

run_test("Sequential Optimization (AHBO)", function() {
  set.seed(791)
  
  result <- bayes_optimize(
    fn = sphere,
    bounds = bounds_2d,
    method = "adaptive_hybrid",
    acquisition = "ei",
    n_iter = 15,
    n_init = 8,
    nu_schedule = "linear",
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  cat(sprintf("  Convergence: %.4f → %.4f\n", 
              result$convergence[1], result$convergence[15]))
  
  stopifnot(result$convergence[15] >= result$convergence[1])
})

run_test("Batch Optimization with Distance Blocking", function() {
  set.seed(792)
  
  result <- bayes_optimize(
    fn = sphere,
    bounds = bounds_2d,
    method = "adaptive_hybrid",
    acquisition = "ei",
    n_iter = 10,
    n_init = 8,
    batch_size = 3,
    blocking_strategy = "distance_soft",
    radius = 0.15,
    penalty_strength = 5.0,
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  cat(sprintf("  Total evaluations: %d\n", nrow(result$X)))
  
  stopifnot(nrow(result$X) == 8 + 10 * 3)
})

run_test("Batch Optimization with Targeted Blocking", function() {
  set.seed(793)
  
  result <- bayes_optimize(
    fn = branin,
    bounds = bounds_2d,
    method = "adaptive_hybrid",
    acquisition = "ei",
    n_iter = 10,
    n_init = 8,
    batch_size = 3,
    blocking_strategy = "targeted",
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  stopifnot(!is.null(result$best_y))
})

# ============================================================================
# TEST 6: Different Acquisition Functions
# ============================================================================

run_test("Optimization with UCB", function() {
  set.seed(800)
  
  result <- bayes_optimize(
    fn = sphere,
    bounds = bounds_2d,
    method = "gaussian",
    acquisition = "ucb",
    n_iter = 12,
    n_init = 8,
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  stopifnot(!is.null(result$best_y))
})

run_test("Optimization with PI", function() {
  set.seed(801)
  
  result <- bayes_optimize(
    fn = sphere,
    bounds = bounds_2d,
    method = "gaussian",
    acquisition = "pi",
    n_iter = 12,
    n_init = 8,
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  stopifnot(!is.null(result$best_y))
})

run_test("Optimization with Thompson Sampling", function() {
  set.seed(802)
  
  result <- bayes_optimize(
    fn = sphere,
    bounds = bounds_2d,
    method = "gaussian",
    acquisition = "thompson",
    n_iter = 12,
    n_init = 8,
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  stopifnot(!is.null(result$best_y))
})

# ============================================================================
# TEST 7: Edge Cases
# ============================================================================

run_test("1D Optimization", function() {
  set.seed(900)
  
  bounds_1d <- matrix(c(0, 1), nrow = 1, ncol = 2)
  fn_1d <- function(x) -(x[1] - 0.5)^2
  
  result <- bayes_optimize(
    fn = fn_1d,
    bounds = bounds_1d,
    method = "gaussian",
    n_iter = 10,
    n_init = 3,
    verbose = FALSE
  )
  
  cat(sprintf("  Best point: %.3f (optimal: 0.5)\n", result$best_x))
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  
  stopifnot(abs(result$best_x - 0.5) < 0.2)
})

run_test("High-Dimensional (5D) Optimization", function() {
  set.seed(901)
  
  bounds_5d <- matrix(c(rep(0, 5), rep(1, 5)), nrow = 5, ncol = 2)
  fn_5d <- function(x) -(sum((x - 0.5)^2))
  
  result <- bayes_optimize(
    fn = fn_5d,
    bounds = bounds_5d,
    method = "adaptive_hybrid",
    n_iter = 15,
    n_init = 15,
    batch_size = 2,
    blocking_strategy = "distance_soft",
    verbose = FALSE
  )
  
  cat(sprintf("  Best value: %.4f\n", result$best_y))
  cat(sprintf("  Dimensions: %d\n", length(result$best_x)))
  
  stopifnot(length(result$best_x) == 5)
})

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat(rep("=", 80), "\n", sep="")
cat("TEST SUMMARY\n")
cat(rep("=", 80), "\n", sep="")
cat(sprintf("Total tests run:    %d\n", test_num))
cat(sprintf("Tests passed:       %d (%.1f%%)\n", 
            tests_passed, 100 * tests_passed / test_num))
cat(sprintf("Tests failed:       %d (%.1f%%)\n", 
            tests_failed, 100 * tests_failed / test_num))
cat(rep("=", 80), "\n", sep="")

if (tests_failed == 0) {
  cat("\n[SUCCESS] ALL TESTS PASSED!\n\n")
  quit(status = 0)
} else {
  cat("\n[ERROR] SOME TESTS FAILED\n\n")
  quit(status = 1)
}
