# Example: Diverse Batch Selection
# Demonstrates the batch selection utilities

library(BayesOptR)

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# Example 1: Basic diverse batch selection
# ============================================================================

cat("Example 1: Basic Diverse Batch Selection\n")
cat("==========================================\n\n")

# Generate sample data (3D observations)
n_samples <- 100
X <- matrix(rnorm(n_samples * 3), ncol = 3)

# Create scores (e.g., negative Rosenbrock function)
scores <- apply(X, 1, function(x) {
  -sum((1 - x[1:(length(x)-1)])^2 + 100 * (x[2:length(x)] - x[1:(length(x)-1)]^2)^2)
})

cat("Total samples:", n_samples, "\n")
cat("Score range:", round(min(scores), 2), "to", round(max(scores), 2), "\n\n")

# Select 10 diverse high-scoring points (maximize scores)
indices_max <- select_diverse_batch(
  X, scores, 
  batch_size = 10, 
  minimize = FALSE,
  diversity_weight = 0.5,
  method = "greedy"
)

cat("Selected indices (maximize):", indices_max, "\n")
cat("Selected scores:", round(scores[indices_max], 2), "\n")
cat("Mean selected score:", round(mean(scores[indices_max]), 2), "\n")
cat("Mean all scores:", round(mean(scores), 2), "\n\n")

# ============================================================================
# Example 2: Different diversity weights
# ============================================================================

cat("Example 2: Effect of Diversity Weight\n")
cat("======================================\n\n")

# More emphasis on diversity (0.8)
indices_div <- select_diverse_batch(
  X, scores, 
  batch_size = 10, 
  minimize = FALSE,
  diversity_weight = 0.8,
  method = "greedy"
)

# More emphasis on score quality (0.2)
indices_score <- select_diverse_batch(
  X, scores, 
  batch_size = 10, 
  minimize = FALSE,
  diversity_weight = 0.2,
  method = "greedy"
)

# Compute average pairwise distance for each selection
compute_avg_distance <- function(X, indices) {
  X_selected <- X[indices, , drop = FALSE]
  n <- nrow(X_selected)
  if (n <= 1) return(0)
  
  dists <- as.matrix(dist(X_selected))
  sum(dists) / (n * (n - 1))
}

dist_balanced <- compute_avg_distance(X, indices_max)
dist_diverse <- compute_avg_distance(X, indices_div)
dist_score <- compute_avg_distance(X, indices_score)

cat("Diversity weight = 0.5 (balanced):\n")
cat("  Mean score:", round(mean(scores[indices_max]), 2), "\n")
cat("  Avg pairwise distance:", round(dist_balanced, 3), "\n\n")

cat("Diversity weight = 0.8 (diversity focus):\n")
cat("  Mean score:", round(mean(scores[indices_div]), 2), "\n")
cat("  Avg pairwise distance:", round(dist_diverse, 3), "\n\n")

cat("Diversity weight = 0.2 (score focus):\n")
cat("  Mean score:", round(mean(scores[indices_score]), 2), "\n")
cat("  Avg pairwise distance:", round(dist_score, 3), "\n\n")

# ============================================================================
# Example 3: Different methods
# ============================================================================

cat("Example 3: Comparison of Selection Methods\n")
cat("==========================================\n\n")

# Greedy method
indices_greedy <- select_diverse_batch(
  X, scores, batch_size = 10, minimize = FALSE,
  diversity_weight = 0.5, method = "greedy"
)

# K-means method
indices_kmeans <- select_diverse_batch(
  X, scores, batch_size = 10, minimize = FALSE,
  diversity_weight = 0.5, method = "kmeans"
)

# Farthest point method
indices_farthest <- select_diverse_batch(
  X, scores, batch_size = 10, minimize = FALSE,
  diversity_weight = 0.5, method = "farthest"
)

cat("Greedy method:\n")
cat("  Mean score:", round(mean(scores[indices_greedy]), 2), "\n")
cat("  Avg distance:", round(compute_avg_distance(X, indices_greedy), 3), "\n\n")

cat("K-means method:\n")
cat("  Mean score:", round(mean(scores[indices_kmeans]), 2), "\n")
cat("  Avg distance:", round(compute_avg_distance(X, indices_kmeans), 3), "\n\n")

cat("Farthest point method:\n")
cat("  Mean score:", round(mean(scores[indices_farthest]), 2), "\n")
cat("  Avg distance:", round(compute_avg_distance(X, indices_farthest), 3), "\n\n")

# ============================================================================
# Example 4: Minimization vs Maximization
# ============================================================================

cat("Example 4: Minimization vs Maximization\n")
cat("=======================================\n\n")

# Create optimization problem scores
scores_opt <- -rowSums(X^2)  # Minimize sum of squares

# Select for minimization
indices_min <- select_diverse_batch(
  X, scores_opt, 
  batch_size = 10, 
  minimize = TRUE,
  diversity_weight = 0.5
)

# Select for maximization (same scores, different interpretation)
indices_max <- select_diverse_batch(
  X, scores_opt, 
  batch_size = 10, 
  minimize = FALSE,
  diversity_weight = 0.5
)

cat("Minimize (seeking low scores):\n")
cat("  Selected scores:", round(scores_opt[indices_min], 2), "\n")
cat("  Mean:", round(mean(scores_opt[indices_min]), 2), "\n\n")

cat("Maximize (seeking high scores):\n")
cat("  Selected scores:", round(scores_opt[indices_max], 2), "\n")
cat("  Mean:", round(mean(scores_opt[indices_max]), 2), "\n\n")

# ============================================================================
# Example 5: Pure diversity (max-min distance)
# ============================================================================

cat("Example 5: Maximum Diversity Selection\n")
cat("======================================\n\n")

indices_pure_div <- select_max_diversity_batch(X, batch_size = 10)

cat("Pure diversity selection (ignoring scores):\n")
cat("  Selected indices:", indices_pure_div, "\n")
cat("  Avg pairwise distance:", round(compute_avg_distance(X, indices_pure_div), 3), "\n")
cat("  Mean score:", round(mean(scores[indices_pure_div]), 2), "\n\n")

# Compare with score-focused selection
top_indices <- order(scores, decreasing = TRUE)[1:10]
cat("Top 10 by score only:\n")
cat("  Avg pairwise distance:", round(compute_avg_distance(X, top_indices), 3), "\n")
cat("  Mean score:", round(mean(scores[top_indices]), 2), "\n\n")

# ============================================================================
# Example 6: Score-weighted sampling
# ============================================================================

cat("Example 6: Score-Weighted Sampling\n")
cat("==================================\n\n")

# Sample with different temperatures
indices_cold <- select_score_weighted_batch(
  X, scores, 
  batch_size = 10, 
  minimize = FALSE,
  temperature = 0.5  # Lower temp = more deterministic (favor high scores)
)

indices_hot <- select_score_weighted_batch(
  X, scores, 
  batch_size = 10, 
  minimize = FALSE,
  temperature = 2.0  # Higher temp = more random
)

cat("Score-weighted sampling (temperature = 0.5, greedy):\n")
cat("  Mean score:", round(mean(scores[indices_cold]), 2), "\n\n")

cat("Score-weighted sampling (temperature = 2.0, exploratory):\n")
cat("  Mean score:", round(mean(scores[indices_hot]), 2), "\n\n")

# ============================================================================
# Example 7: Return actual points instead of indices
# ============================================================================

cat("Example 7: Returning Actual Points\n")
cat("===================================\n\n")

result <- select_diverse_batch(
  X, scores, 
  batch_size = 5, 
  minimize = FALSE,
  diversity_weight = 0.6,
  return_indices = FALSE
)

cat("Selected batch (first 3 points):\n")
print(round(result$X[1:3, ], 3))
cat("\nSelected scores:\n")
print(round(result$scores, 2))

cat("\n[SUCCESS] All examples completed!\n")
