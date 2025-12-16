#' Gaussian Process Bayesian Optimization
#'
#' @description
#' Standard Bayesian Optimization using Gaussian Process surrogate model.
#' Assumes Gaussian likelihood and implements GP regression with various kernel options.
#'
#' @details
#' This class implements classic GP-based Bayesian Optimization with:
#' - Gaussian likelihood (normality assumption)
#' - Various kernel functions (Matern, RBF, etc.)
#' - Automatic hyperparameter optimization via MLE
#' - Numerical stability through jitter addition
#'
#' @field bounds Matrix of bounds (lower, upper) for each dimension
#' @field kernel Kernel function type ("matern52", "matern32", "rbf")
#' @field nugget Nugget parameter for numerical stability (default: 1e-8)
#' @field X Matrix of observed input points
#' @field y Vector of observed function values
#' @field gp_model Fitted GP model object
#'
#' @export
GaussianBO <- R6::R6Class(
  "GaussianBO",
  
  public = list(
    bounds = NULL,
    kernel = NULL,
    nugget = NULL,
    X = NULL,
    y = NULL,
    gp_model = NULL,
    n_init = NULL,
    
    #' @description
    #' Initialize Gaussian BO
    #' @param bounds Matrix with 2 columns (lower, upper) and n_dim rows
    #' @param kernel Kernel type: "matern52" (default), "matern32", "rbf"
    #' @param nugget Nugget for numerical stability (default: 1e-8)
    #' @param n_init Number of initial random samples (default: 5*n_dim)
    initialize = function(bounds, kernel = "matern52", nugget = 1e-8, n_init = NULL) {
      if (!is.matrix(bounds) || ncol(bounds) != 2) {
        stop("bounds must be a matrix with 2 columns (lower, upper)")
      }
      
      self$bounds <- bounds
      self$kernel <- kernel
      self$nugget <- nugget
      self$n_init <- if (is.null(n_init)) 5 * nrow(bounds) else n_init
      self$X <- matrix(nrow = 0, ncol = nrow(bounds))
      self$y <- numeric(0)
    },
    
    #' @description
    #' Add observation to the dataset
    #' @param x Input point (vector)
    #' @param y_obs Observed function value
    add_observation = function(x, y_obs) {
      if (length(x) != nrow(self$bounds)) {
        stop("x dimension mismatch with bounds")
      }
      
      self$X <- rbind(self$X, matrix(x, nrow = 1))
      self$y <- c(self$y, y_obs)
      
      # Refit GP model
      self$fit_gp()
    },
    
    #' @description
    #' Fit Gaussian Process model to current data
    fit_gp = function() {
      if (nrow(self$X) < 2) {
        return(NULL)
      }
      
      tryCatch({
        # Use DiceKriging for GP fitting
        self$gp_model <- DiceKriging::km(
          formula = ~1,
          design = as.data.frame(self$X),
          response = self$y,
          covtype = self$kernel,
          nugget = self$nugget,
          control = list(trace = FALSE)
        )
      }, error = function(e) {
        warning(paste("GP fitting failed:", e$message))
        self$gp_model <- NULL
      })
    },
    
    #' @description
    #' Predict mean and variance at new points
    #' @param X_new Matrix of new points to predict
    #' @return List with mean and variance vectors
    predict = function(X_new) {
      if (is.null(self$gp_model)) {
        # Return uninformative predictions
        n <- nrow(X_new)
        return(list(
          mean = rep(mean(self$y), n),
          sd = rep(sd(self$y), n),
          variance = rep(var(self$y), n)
        ))
      }
      
      pred <- DiceKriging::predict.km(
        self$gp_model,
        newdata = as.data.frame(X_new),
        type = "UK",
        checkNames = FALSE
      )
      
      list(
        mean = pred$mean,
        sd = pred$sd,
        variance = pred$sd^2
      )
    },
    
    #' @description
    #' Get current best observation
    #' @return List with best_x and best_y
    get_best = function() {
      if (length(self$y) == 0) {
        return(list(best_x = NULL, best_y = -Inf))
      }
      
      best_idx <- which.max(self$y)
      list(
        best_x = self$X[best_idx, ],
        best_y = self$y[best_idx]
      )
    },
    
    #' @description
    #' Generate initial random samples using Latin Hypercube Sampling
    #' @return Matrix of initial samples
    generate_initial_samples = function() {
      n_dim <- nrow(self$bounds)
      
      # Latin Hypercube Sampling
      lhs_samples <- lhs::randomLHS(self$n_init, n_dim)
      
      # Scale to bounds
      X_init <- matrix(nrow = self$n_init, ncol = n_dim)
      for (i in 1:n_dim) {
        X_init[, i] <- self$bounds[i, 1] + 
          lhs_samples[, i] * (self$bounds[i, 2] - self$bounds[i, 1])
      }
      
      X_init
    }
  )
)
