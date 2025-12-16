#' Plot Convergence Curve
#'
#' @description
#' Plot convergence history of Bayesian Optimization.
#'
#' @param result bayesopt result object
#' @param log_scale Use log scale for y-axis (default: FALSE)
#' @param ... Additional arguments passed to plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- bayes_optimize(fn, bounds, n_iter = 50)
#' plot_convergence(result)
#' }
plot_convergence <- function(result, log_scale = FALSE, ...) {
  if (!inherits(result, "bayesopt")) {
    stop("result must be a bayesopt object")
  }
  
  iterations <- 1:length(result$convergence)
  
  if (log_scale) {
    plot(iterations, result$convergence, type = "l", lwd = 2, col = "blue",
         xlab = "Iteration", ylab = "Best Value (log scale)",
         main = "Bayesian Optimization Convergence",
         log = "y", ...)
  } else {
    plot(iterations, result$convergence, type = "l", lwd = 2, col = "blue",
         xlab = "Iteration", ylab = "Best Value",
         main = "Bayesian Optimization Convergence", ...)
  }
  
  grid()
  
  # Add markers for improvements
  improvements <- which(diff(c(-Inf, result$convergence)) > 0)
  points(improvements, result$convergence[improvements], 
         pch = 19, col = "red", cex = 1.2)
  
  # Add legend
  legend("bottomright", 
         legend = c("Best value", "Improvement"),
         col = c("blue", "red"),
         lty = c(1, NA),
         pch = c(NA, 19),
         lwd = c(2, NA))
}


#' Plot Acquisition Function
#'
#' @description
#' Plot acquisition function surface for 2D problems.
#'
#' @param model BO model object
#' @param acquisition_fn Acquisition function object
#' @param resolution Number of points per dimension (default: 50)
#' @param ... Additional arguments passed to contour
#'
#' @export
plot_acquisition <- function(model, acquisition_fn, resolution = 50, ...) {
  if (nrow(model$bounds) != 2) {
    stop("Can only plot acquisition for 2D problems")
  }
  
  # Generate grid
  x1_seq <- seq(model$bounds[1, 1], model$bounds[1, 2], length.out = resolution)
  x2_seq <- seq(model$bounds[2, 1], model$bounds[2, 2], length.out = resolution)
  
  grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)
  
  # Get predictions
  pred <- model$predict(as.matrix(grid))
  
  # Evaluate acquisition
  f_best <- model$get_best()$best_y
  
  if (inherits(acquisition_fn, "ThompsonSampling")) {
    acq_values <- acquisition_fn$evaluate(pred$mean, pred$sd)
  } else if (inherits(acquisition_fn, "UpperConfidenceBound")) {
    acq_values <- acquisition_fn$evaluate(pred$mean, pred$sd, 2)
  } else {
    acq_values <- acquisition_fn$evaluate(pred$mean, pred$sd, f_best)
  }
  
  # Reshape for contour plot
  acq_matrix <- matrix(acq_values, nrow = resolution, ncol = resolution)
  
  # Plot
  contour(x1_seq, x2_seq, acq_matrix, 
          xlab = "x1", ylab = "x2",
          main = "Acquisition Function",
          nlevels = 20,
          ...)
  
  # Add observed points
  points(model$X[, 1], model$X[, 2], pch = 19, col = "red", cex = 1.2)
  
  # Highlight best point
  best <- model$get_best()
  points(best$best_x[1], best$best_x[2], pch = 8, col = "green", cex = 2, lwd = 2)
  
  legend("topright",
         legend = c("Observed", "Best"),
         col = c("red", "green"),
         pch = c(19, 8))
}


#' Plot Model Predictions
#'
#' @description
#' Plot predictive mean and uncertainty for 2D problems.
#'
#' @param model BO model object
#' @param type Type of plot: "mean", "uncertainty", "both" (default: "both")
#' @param resolution Number of points per dimension (default: 50)
#'
#' @export
plot_predictions <- function(model, type = "both", resolution = 50) {
  if (nrow(model$bounds) != 2) {
    stop("Can only plot for 2D problems")
  }
  
  # Generate grid
  x1_seq <- seq(model$bounds[1, 1], model$bounds[1, 2], length.out = resolution)
  x2_seq <- seq(model$bounds[2, 1], model$bounds[2, 2], length.out = resolution)
  grid <- expand.grid(x1 = x1_seq, x2 = x2_seq)
  
  # Get predictions
  pred <- model$predict(as.matrix(grid))
  
  mean_matrix <- matrix(pred$mean, nrow = resolution, ncol = resolution)
  sd_matrix <- matrix(pred$sd, nrow = resolution, ncol = resolution)
  
  if (type %in% c("both", "mean")) {
    # Plot mean
    if (type == "both") par(mfrow = c(1, 2))
    
    contour(x1_seq, x2_seq, mean_matrix,
            xlab = "x1", ylab = "x2",
            main = "Predictive Mean",
            nlevels = 20)
    points(model$X[, 1], model$X[, 2], pch = 19, col = "red")
  }
  
  if (type %in% c("both", "uncertainty")) {
    # Plot uncertainty
    contour(x1_seq, x2_seq, sd_matrix,
            xlab = "x1", ylab = "x2",
            main = "Predictive Std Dev",
            nlevels = 20)
    points(model$X[, 1], model$X[, 2], pch = 19, col = "red")
  }
  
  if (type == "both") par(mfrow = c(1, 1))
}


#' Export Results to CSV
#'
#' @description
#' Export optimization results to CSV file.
#'
#' @param result bayesopt result object
#' @param file Output file path
#'
#' @export
export_results <- function(result, file) {
  df <- data.frame(
    iteration = rep(1:result$n_iter, each = result$batch_size),
    batch = rep(1:result$batch_size, result$n_iter),
    result$X,
    y = result$y,
    best_so_far = rep(result$convergence, each = result$batch_size)
  )
  
  colnames(df)[3:(2 + ncol(result$X))] <- paste0("x", 1:ncol(result$X))
  
  write.csv(df, file, row.names = FALSE)
  message("Results exported to: ", file)
}
