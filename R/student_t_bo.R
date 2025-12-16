#' Student-t Process Bayesian Optimization
#'
#' @description
#' Bayesian Optimization using Student-t Process as surrogate model.
#' More robust to outliers compared to Gaussian Process due to heavy tails.
#'
#' @details
#' Student-t Process provides:
#' - Heavy-tailed predictions (robust to outliers)
#' - Controlled by degrees of freedom parameter (nu)
#' - nu=3: Very heavy tails, strong exploration
#' - nu=5: Moderate tails, balanced
#' - nu->infinity: Converges to Gaussian Process
#'
#' @field nu Degrees of freedom parameter (controls tail heaviness)
#' @inheritParams GaussianBO
#'
#' @export
StudentTBO <- R6::R6Class(
  "StudentTBO",
  inherit = GaussianBO,
  
  public = list(
    nu = NULL,
    
    #' @description
    #' Initialize Student-t BO
    #' @param nu Degrees of freedom (default: 3, heavier tails)
    #' @inheritParams GaussianBO
    initialize = function(bounds, nu = 3, kernel = "matern52", 
                         nugget = 1e-8, n_init = NULL) {
      super$initialize(bounds, kernel, nugget, n_init)
      self$nu <- nu
    },
    
    #' @description
    #' Predict with Student-t Process
    #' @param X_new Matrix of new points
    #' @return List with mean, variance, and degrees of freedom
    predict = function(X_new) {
      # First get GP predictions
      gp_pred <- super$predict(X_new)
      
      if (is.null(self$gp_model)) {
        return(c(gp_pred, list(nu = self$nu)))
      }
      
      n <- nrow(self$X)
      nu_posterior <- self$nu + n
      
      # Student-t variance adjustment
      variance_scale <- (self$nu + mean((self$y - mean(self$y))^2)) / 
                       (self$nu + n - 2)
      
      list(
        mean = gp_pred$mean,
        variance = gp_pred$variance * variance_scale,
        sd = sqrt(gp_pred$variance * variance_scale),
        nu = nu_posterior
      )
    },
    
    #' @description
    #' Sample from Student-t predictive distribution
    #' @param X_new Points to sample at
    #' @param n_samples Number of samples
    #' @return Matrix of samples (n_samples × nrow(X_new))
    sample = function(X_new, n_samples = 1) {
      pred <- self$predict(X_new)
      
      n_points <- nrow(X_new)
      samples <- matrix(nrow = n_samples, ncol = n_points)
      
      for (i in 1:n_points) {
        # Student-t sampling: use scaled chi-square
        z <- rnorm(n_samples)
        chi_sq <- rchisq(n_samples, df = pred$nu)
        
        samples[, i] <- pred$mean[i] + 
          pred$sd[i] * z * sqrt(pred$nu / chi_sq)
      }
      
      samples
    }
  )
)
