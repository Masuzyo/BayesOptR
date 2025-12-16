# Example 2: Batch Optimization with Blocking Strategies
# Demonstrates different blocking strategies for batch acquisition

library(BayesOptR)

# Define Ackley function (3D, multimodal)
ackley3d <- function(x) {
  n <- length(x)
  sum1 <- sum(x^2)
  sum2 <- sum(cos(2 * pi * x))
  -(-20 * exp(-0.2 * sqrt(sum1/n)) - exp(sum2/n) + 20 + exp(1))
}

# Bounds: [-5, 5]^3
bounds <- matrix(c(-5, -5, -5, 5, 5, 5), nrow = 3, ncol = 2)

cat("========================================\n")
cat("Example 2: Batch Optimization with Blocking\n")
cat("========================================\n\n")

# Test different blocking strategies
strategies <- c("none", "dimwise_hard", "dimwise_soft", 
                "distance_hard", "distance_soft", "local_pen")

results <- list()

for (strategy in strategies) {
  cat(sprintf("\nTesting: %s\n", strategy))
  cat(rep("=", 40), "\n", sep = "")
  
  results[[strategy]] <- bayes_optimize(
    fn = ackley3d,
    bounds = bounds,
    method = "adaptive_hybrid",
    acquisition = "ei",
    batch_size = 4,
    blocking_strategy = strategy,
    n_iter = 20,  # 20 iterations × 4 batch = 80 evaluations
    minimize = TRUE,
    verbose = FALSE
  )
  
  cat(sprintf("Final best: %.6f\n", results[[strategy]]$best_y))
}

# Summary comparison
cat("\n\n========================================\n")
cat("Blocking Strategy Comparison\n")
cat("========================================\n")
cat(sprintf("%-20s %15s %15s\n", "Strategy", "Best Value", "Evaluations"))
cat(rep("-", 50), "\n", sep = "")

for (strategy in strategies) {
  cat(sprintf("%-20s %15.6f %15d\n", 
              strategy,
              results[[strategy]]$best_y,
              nrow(results[[strategy]]$X)))
}

cat(sprintf("\nKnown optimum: %.6f\n", 0.0))

# Plot convergence for all strategies
par(mfrow = c(2, 3))
for (strategy in strategies) {
  plot(1:length(results[[strategy]]$convergence),
       results[[strategy]]$convergence,
       type = "l", lwd = 2, col = "blue",
       xlab = "Iteration", ylab = "Best Value",
       main = sprintf("Strategy: %s", strategy))
  grid()
}

# Detailed comparison plot
par(mfrow = c(1, 1))
colors <- rainbow(length(strategies))
plot(1:20, results[[strategies[1]]]$convergence, type = "l", lwd = 2,
     col = colors[1], xlab = "Iteration", ylab = "Best Value (minimization)",
     main = "Convergence Comparison: Blocking Strategies",
     ylim = range(sapply(results, function(r) range(r$convergence))))

for (i in 2:length(strategies)) {
  lines(1:20, results[[strategies[i]]]$convergence, 
        lwd = 2, col = colors[i])
}

legend("topright", legend = strategies, col = colors, lwd = 2, cex = 0.8)
grid()

# Analyze batch diversity
cat("\n\nBatch Diversity Analysis\n")
cat("========================================\n")

for (strategy in strategies) {
  X <- results[[strategy]]$X
  
  # Calculate pairwise distances in each batch
  n_batches <- nrow(X) / 4
  diversities <- numeric(n_batches)
  
  for (i in 1:n_batches) {
    batch_idx <- ((i-1)*4 + 1):(i*4)
    batch_points <- X[batch_idx, ]
    
    # Mean pairwise distance
    dists <- dist(batch_points)
    diversities[i] <- mean(dists)
  }
  
  cat(sprintf("%-20s Mean diversity: %.4f (SD: %.4f)\n",
              strategy, mean(diversities), sd(diversities)))
}

cat("\nExample completed!\n")
