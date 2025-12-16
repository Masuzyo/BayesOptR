# Example 3: Parallel Batch Optimization
# Demonstrates parallel evaluation with different backends

library(BayesOptR)

# Define expensive Rosenbrock function with artificial delay
rosenbrock <- function(x) {
  Sys.sleep(0.1)  # Simulate expensive evaluation
  n <- length(x)
  sum(100 * (x[2:n] - x[1:(n-1)]^2)^2 + (1 - x[1:(n-1)])^2)
}

# Bounds: [-2, 2]^5
bounds <- matrix(c(rep(-2, 5), rep(2, 5)), nrow = 5, ncol = 2)

cat("========================================\n")
cat("Example 3: Parallel Batch Optimization\n")
cat("========================================\n\n")

# Test sequential vs parallel

cat("1. Sequential evaluation (batch_size=1)\n")
cat(rep("-", 40), "\n", sep = "")
t1 <- system.time({
  result_seq <- bayes_optimize(
    fn = rosenbrock,
    bounds = bounds,
    method = "adaptive_hybrid",
    batch_size = 1,
    n_iter = 10,
    minimize = TRUE,
    verbose = FALSE
  )
})

cat(sprintf("Time: %.2f seconds\n", t1[3]))
cat(sprintf("Best value: %.6f\n", result_seq$best_y))
cat(sprintf("Evaluations: %d\n\n", nrow(result_seq$X)))

cat("2. Parallel batch evaluation (batch_size=4, future backend)\n")
cat(rep("-", 40), "\n", sep = "")
t2 <- system.time({
  result_parallel_future <- bayes_optimize(
    fn = rosenbrock,
    bounds = bounds,
    method = "adaptive_hybrid",
    batch_size = 4,
    blocking_strategy = "dimwise_soft",
    parallel = TRUE,
    parallel_backend = "future",
    n_cores = 4,
    n_iter = 10,
    minimize = TRUE,
    verbose = FALSE
  )
})

cat(sprintf("Time: %.2f seconds\n", t2[3]))
cat(sprintf("Best value: %.6f\n", result_parallel_future$best_y))
cat(sprintf("Evaluations: %d\n", nrow(result_parallel_future$X)))
cat(sprintf("Speedup: %.2fx\n\n", t1[3] / t2[3]))

# Test different backends
if (requireNamespace("foreach", quietly = TRUE) && 
    requireNamespace("doParallel", quietly = TRUE)) {
  
  cat("3. Parallel batch evaluation (foreach backend)\n")
  cat(rep("-", 40), "\n", sep = "")
  t3 <- system.time({
    result_parallel_foreach <- bayes_optimize(
      fn = rosenbrock,
      bounds = bounds,
      method = "adaptive_hybrid",
      batch_size = 4,
      blocking_strategy = "dimwise_soft",
      parallel = TRUE,
      parallel_backend = "foreach",
      n_cores = 4,
      n_iter = 10,
      minimize = TRUE,
      verbose = FALSE
    )
  })
  
  cat(sprintf("Time: %.2f seconds\n", t3[3]))
  cat(sprintf("Best value: %.6f\n", result_parallel_foreach$best_y))
  cat(sprintf("Speedup: %.2fx\n\n", t1[3] / t3[3]))
}

# Performance summary
cat("\n========================================\n")
cat("Performance Summary\n")
cat("========================================\n")
cat(sprintf("Sequential (1 core):     %.2f seconds\n", t1[3]))
cat(sprintf("Parallel future (4 cores): %.2f seconds (%.2fx speedup)\n", 
            t2[3], t1[3]/t2[3]))

if (exists("t3")) {
  cat(sprintf("Parallel foreach (4 cores): %.2f seconds (%.2fx speedup)\n", 
              t3[3], t1[3]/t3[3]))
}

cat(sprintf("\nTheoretical speedup with 4 cores: ~4x\n"))
cat(sprintf("Achieved speedup: %.2fx\n", t1[3]/t2[3]))
cat(sprintf("Efficiency: %.1f%%\n", (t1[3]/t2[3])/4 * 100))

# Plot convergence comparison
par(mfrow = c(1, 1))
plot(1:10, result_seq$convergence, type = "l", lwd = 2, col = "blue",
     xlab = "Iteration", ylab = "Best Value (minimize)",
     main = "Sequential vs Parallel Convergence",
     ylim = range(c(result_seq$convergence, result_parallel_future$convergence)))
lines(1:10, result_parallel_future$convergence, lwd = 2, col = "red")

legend("topright",
       legend = c("Sequential", "Parallel (batch=4)"),
       col = c("blue", "red"),
       lwd = 2)
grid()

cat("\n\nNote: Parallel speedup depends on:\n")
cat("  - Function evaluation time (overhead vs computation)\n")
cat("  - Number of cores available\n")
cat("  - Batch size (larger = better parallelization)\n")
cat("  - Backend efficiency\n")

cat("\nExample completed!\n")
