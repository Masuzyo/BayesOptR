#' Main Bayesian Optimization Interface
#'
#' @description
#' High-level interface for running Bayesian Optimization with various configurations.
#'
#' @param fn Objective function to maximize (or minimize if minimize=TRUE)
#' @param bounds Matrix with 2 columns (lower, upper) and n_dim rows
#' @param method BO method: "gaussian", "student_t", "adaptive_hybrid" (default: "adaptive_hybrid")
#' @param acquisition Acquisition function: "ei", "ucb", "pi", "thompson" (default: "ei")
#' @param n_iter Number of iterations (default: 50)
#' @param n_init Number of initial random samples (default: 5 * n_dim)
#' @param batch_size Batch size for parallel evaluation (default: 1, sequential)
#' @param blocking_strategy Blocking strategy: "none", "dimwise_hard", "dimwise_soft",
#'   "distance_hard", "distance_soft", "targeted", "local_pen", "constant_liar" (default: "none")
#' @param parallel Enable parallel evaluation (default: FALSE)
#' @param parallel_backend Backend for parallelization (default: "future")
#' @param n_cores Number of cores (default: detectCores() - 1)
#' @param minimize Minimize instead of maximize (default: FALSE)
#' @param verbose Print progress (default: TRUE)
#' @param nu Degrees of freedom for Student-t (default: 3)
#' @param nu_schedule Schedule for adaptive hybrid: "data_driven", "linear", "exponential", "sigmoid", "performance" (default: "data_driven")
#' @param ... Additional arguments passed to blocking strategies
#'
#' @return List with optimization results including:
#'   - best_x: Best point found
#'   - best_y: Best value found
#'   - X: All evaluated points
#'   - y: All function values
#'   - convergence: Convergence history
#'   - model: Final BO model
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Branin function
#' branin <- function(x) {
#'   x1 <- x[1] * 15 - 5
#'   x2 <- x[2] * 15
#'   a <- 1
#'   b <- 5.1 / (4 * pi^2)
#'   c <- 5 / pi
#'   r <- 6
#'   s <- 10
#'   t <- 1 / (8 * pi)
#'   -(a * (x2 - b * x1^2 + c * x1 - r)^2 + s * (1 - t) * cos(x1) + s)
#' }
#'
#' bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
#'
#' # Standard BO
#' result <- bayes_optimize(branin, bounds, method = "adaptive_hybrid", n_iter = 50)
#'
#' # Batch BO with dimension-wise blocking
#' result <- bayes_optimize(branin, bounds, method = "adaptive_hybrid",
#'                         batch_size = 4, blocking_strategy = "dimwise_soft",
#'                         parallel = TRUE, n_iter = 25)
#' }
bayes_optimize <- function(fn, bounds, 
                          method = "adaptive_hybrid",
                          acquisition = "ei",
                          n_iter = 50,
                          n_init = NULL,
                          batch_size = 1,
                          blocking_strategy = "none",
                          parallel = FALSE,
                          parallel_backend = "future",
                          n_cores = NULL,
                          minimize = FALSE,
                          verbose = TRUE,
                          nu = 3,
                          nu_schedule = "data_driven",
                          ...) {
  
  # Validate inputs
  if (!is.matrix(bounds) || ncol(bounds) != 2) {
    stop("bounds must be a matrix with 2 columns (lower, upper)")
  }
  
  n_dim <- nrow(bounds)
  if (is.null(n_init)) {
    n_init <- 5 * n_dim
  }
  
  # Wrap function for minimization
  fn_wrapped <- if (minimize) {
    function(x) -fn(x)
  } else {
    fn
  }
  
  # Initialize model
  model <- switch(
    method,
    "gaussian" = GaussianBO$new(bounds, n_init = n_init),
    "student_t" = StudentTBO$new(bounds, nu = nu, n_init = n_init),
    "adaptive_hybrid" = AdaptiveHybridBO$new(bounds, nu_min = 3, nu_max = 50,
                                            nu_schedule = nu_schedule,
                                            max_iterations = n_iter,
                                            n_init = n_init),
    stop("Unknown method: ", method)
  )
  
  # Initialize acquisition function
  # For AHBO, use Student-t EI by default; for others use Gaussian EI
  acq_fn <- switch(
    acquisition,
    "ei" = {
      if (method == "adaptive_hybrid" || method == "student_t") {
        StudentTExpectedImprovement$new(nu = nu)
      } else {
        ExpectedImprovement$new()
      }
    },
    "ucb" = UpperConfidenceBound$new(kappa = 2.0, adaptive = TRUE),
    "pi" = ProbabilityOfImprovement$new(),
    "thompson" = ThompsonSampling$new(n_samples = 10),
    stop("Unknown acquisition function: ", acquisition)
  )
  
  # Initialize blocking strategy
  blocking <- if (blocking_strategy != "none") {
    switch(
      blocking_strategy,
      "dimwise_hard" = DimensionWiseHardBlocking$new(...),
      "dimwise_soft" = DimensionWiseSoftBlocking$new(...),
      "distance_hard" = DistanceHardBlocking$new(...),
      "distance_soft" = DistanceSoftBlocking$new(...),
      "targeted" = TargetedBlocking$new(),
      "local_pen" = LocalPenalization$new(...),
      "constant_liar" = ConstantLiar$new(...),
      NULL
    )
  } else {
    NULL
  }
  
  # Initialize batch acquisition
  if (parallel && batch_size > 1) {
    batch_acq <- ParallelBatchAcquisition$new(
      model, acq_fn, blocking, batch_size, 
      parallel_backend = parallel_backend,
      n_cores = n_cores
    )
  } else {
    batch_acq <- BatchAcquisition$new(model, acq_fn, blocking, batch_size)
  }
  
  # Initial sampling
  if (verbose) cat(sprintf("Initializing with %d random samples...\n", n_init))
  X_init <- model$generate_initial_samples()
  
  if (parallel && batch_size > 1) {
    y_init <- batch_acq$evaluate_batch_parallel(X_init, fn_wrapped)
  } else {
    y_init <- apply(X_init, 1, fn_wrapped)
  }
  
  for (i in 1:n_init) {
    model$add_observation(X_init[i, ], y_init[i])
  }
  
  # Main optimization loop
  convergence <- numeric(n_iter)
  
  for (iter in 1:n_iter) {
    # Select batch
    batch <- batch_acq$select_batch()
    
    # Evaluate batch
    if (parallel && batch_size > 1) {
      y_batch <- batch_acq$evaluate_batch_parallel(batch, fn_wrapped)
    } else {
      y_batch <- apply(batch, 1, fn_wrapped)
    }
    
    # Add observations
    for (i in 1:batch_size) {
      model$add_observation(batch[i, ], y_batch[i])
    }
    
    # Track convergence
    best <- model$get_best()
    convergence[iter] <- if (minimize) -best$best_y else best$best_y
    
    if (verbose) {
      cat(sprintf("Iteration %d/%d: Best = %.6f\n", 
                  iter, n_iter, convergence[iter]))
    }
  }
  
  # Cleanup parallel resources
  if (parallel && batch_size > 1) {
    batch_acq$cleanup()
  }
  
  # Return results
  best <- model$get_best()
  
  result <- list(
    best_x = best$best_x,
    best_y = if (minimize) -best$best_y else best$best_y,
    X = model$X,
    y = if (minimize) -model$y else model$y,
    convergence = convergence,
    model = model,
    method = method,
    acquisition = acquisition,
    blocking_strategy = blocking_strategy,
    n_iter = n_iter,
    batch_size = batch_size
  )
  
  class(result) <- c("bayesopt", "list")
  result
}


#' Batch Bayesian Optimization (Convenience Wrapper)
#'
#' @description
#' Convenience wrapper for batch Bayesian Optimization with parallelization.
#'
#' @inheritParams bayes_optimize
#' @export
batch_bayes_optimize <- function(fn, bounds, batch_size = 4, 
                                 blocking_strategy = "dimwise_soft",
                                 parallel = TRUE, ...) {
  bayes_optimize(fn, bounds, batch_size = batch_size,
                blocking_strategy = blocking_strategy,
                parallel = parallel, ...)
}


#' Print method for bayesopt objects
#' @param x bayesopt object
#' @param ... Additional arguments
#' @export
print.bayesopt <- function(x, ...) {
  cat("Bayesian Optimization Results\n")
  cat("=============================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Acquisition: %s\n", x$acquisition))
  cat(sprintf("Blocking: %s\n", x$blocking_strategy))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Batch size: %d\n", x$batch_size))
  cat(sprintf("\nBest value: %.6f\n", x$best_y))
  cat("Best point:\n")
  print(x$best_x)
}


#' Summary method for bayesopt objects
#' @param object bayesopt object
#' @param ... Additional arguments
#' @export
summary.bayesopt <- function(object, ...) {
  cat("Bayesian Optimization Summary\n")
  cat("============================\n\n")
  
  print(object)
  
  cat("\n\nConvergence Statistics:\n")
  cat(sprintf("  Initial best: %.6f\n", object$convergence[1]))
  cat(sprintf("  Final best: %.6f\n", tail(object$convergence, 1)))
  cat(sprintf("  Total improvement: %.6f\n", 
              tail(object$convergence, 1) - object$convergence[1]))
  cat(sprintf("  Mean improvement per iteration: %.6f\n",
              mean(diff(object$convergence))))
  
  cat("\n\nEvaluations:\n")
  cat(sprintf("  Total: %d\n", nrow(object$X)))
  cat(sprintf("  Mean value: %.6f\n", mean(object$y)))
  cat(sprintf("  Std dev: %.6f\n", sd(object$y)))
}
