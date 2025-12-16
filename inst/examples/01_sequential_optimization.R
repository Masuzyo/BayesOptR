# Example 1: Simple Sequential Optimization
# Demonstrates basic usage of BayesOptR

library(BayesOptR)

# Define Branin function (standard benchmark)
branin <- function(x) {
  x1 <- x[1] * 15 - 5
  x2 <- x[2] * 15
  a <- 1
  b <- 5.1 / (4 * pi^2)
  c <- 5 / pi
  r <- 6
  s <- 10
  t <- 1 / (8 * pi)
  -(a * (x2 - b * x1^2 + c * x1 - r)^2 + s * (1 - t) * cos(x1) + s)
}

# Define search bounds [0, 1]^2
bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)

cat("========================================\n")
cat("Example 1: Sequential Bayesian Optimization\n")
cat("========================================\n\n")

# Run with Gaussian BO
cat("Running Gaussian BO...\n")
result_gaussian <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "gaussian",
  acquisition = "ei",
  n_iter = 30,
  verbose = TRUE
)

cat("\nGaussian BO Results:\n")
print(result_gaussian)

# Run with Student-t BO
cat("\n\nRunning Student-t BO (nu=3)...\n")
result_student <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "student_t",
  nu = 3,
  acquisition = "ei",
  n_iter = 30,
  verbose = TRUE
)

cat("\nStudent-t BO Results:\n")
print(result_student)

# Run with Adaptive Hybrid BO
cat("\n\nRunning Adaptive Hybrid BO...\n")
result_ahbo <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "adaptive_hybrid",
  nu_schedule = "linear",
  acquisition = "ei",
  n_iter = 30,
  verbose = TRUE
)

cat("\nAHBO Results:\n")
print(result_ahbo)

# Compare results
cat("\n\n========================================\n")
cat("Comparison Summary\n")
cat("========================================\n")
cat(sprintf("Gaussian BO:  Best = %.6f\n", result_gaussian$best_y))
cat(sprintf("Student-t BO: Best = %.6f\n", result_student$best_y))
cat(sprintf("AHBO:         Best = %.6f\n", result_ahbo$best_y))
cat(sprintf("Known optimum: %.6f\n", -0.397887))

# Plot convergence comparison
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  df <- data.frame(
    iteration = rep(1:30, 3),
    best_value = c(result_gaussian$convergence, 
                   result_student$convergence,
                   result_ahbo$convergence),
    method = rep(c("Gaussian BO", "Student-t BO", "AHBO"), each = 30)
  )
  
  p <- ggplot(df, aes(x = iteration, y = best_value, color = method)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(title = "Convergence Comparison on Branin Function",
         x = "Iteration",
         y = "Best Value Found",
         color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  ggsave("branin_convergence_comparison.pdf", width = 8, height = 6)
} else {
  # Fallback to base R plotting
  par(mfrow = c(1, 1))
  plot(1:30, result_gaussian$convergence, type = "l", col = "blue", lwd = 2,
       xlab = "Iteration", ylab = "Best Value",
       main = "Convergence Comparison",
       ylim = range(c(result_gaussian$convergence, 
                     result_student$convergence,
                     result_ahbo$convergence)))
  lines(1:30, result_student$convergence, col = "red", lwd = 2)
  lines(1:30, result_ahbo$convergence, col = "green", lwd = 2)
  legend("bottomright", 
         legend = c("Gaussian BO", "Student-t BO", "AHBO"),
         col = c("blue", "red", "green"),
         lwd = 2)
  grid()
}

cat("\nExample completed!\n")
