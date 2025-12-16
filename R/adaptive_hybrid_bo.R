#' Adaptive Hybrid Bayesian Optimization (AHBO)
#'
#' @description
#' Advanced Bayesian Optimization that adapts the degrees of freedom (nu) parameter
#' dynamically during optimization, transitioning from heavy-tailed exploration to
#' Gaussian-like exploitation.
#'
#' @details
#' AHBO implements an adaptive schedule for the nu parameter:
#' - Early iterations: Low nu (heavy tails) → Strong exploration
#' - Late iterations: High nu → Converges to Gaussian behavior
#' 
#' Adaptation strategies:
#' - "linear": Linear increase from nu_min to nu_max
#' - "exponential": Exponential growth
#' - "sigmoid": Smooth sigmoid transition
#' - "performance": Adapts based on improvement rate
#'
#' @field nu_min Minimum nu value (initial exploration)
#' @field nu_max Maximum nu value (final exploitation)
#' @field nu_schedule Adaptation schedule type
#' @field current_nu Current nu value
#'
#' @export
AdaptiveHybridBO <- R6::R6Class(
  "AdaptiveHybridBO",
  inherit = StudentTBO,
  
  public = list(
    nu_min = NULL,
    nu_max = NULL,
    nu_schedule = NULL,
    current_nu = NULL,
    iteration = 0,
    max_iterations = NULL,
    improvement_history = NULL,
    
    #' @description
    #' Initialize Adaptive Hybrid BO
    #' @param nu_min Minimum nu (default: 3, heavy tails)
    #' @param nu_max Maximum nu (default: 50, near-Gaussian)
    #' @param nu_schedule Schedule type: "linear", "exponential", "sigmoid", "performance"
    #' @param max_iterations Total iterations for adaptation (required for non-performance schedules)
    #' @inheritParams GaussianBO
    initialize = function(bounds, nu_min = 3, nu_max = 50, 
                         nu_schedule = "linear", max_iterations = 100,
                         kernel = "matern52", nugget = 1e-8, n_init = NULL) {
      super$initialize(bounds, nu_min, kernel, nugget, n_init)
      
      self$nu_min <- nu_min
      self$nu_max <- nu_max
      self$nu_schedule <- nu_schedule
      self$current_nu <- nu_min
      self$max_iterations <- max_iterations
      self$improvement_history <- numeric(0)
    },
    
    #' @description
    #' Update nu parameter based on schedule
    #' @param iteration Current iteration number
    update_nu = function(iteration = NULL) {
      if (!is.null(iteration)) {
        self$iteration <- iteration
      } else {
        self$iteration <- self$iteration + 1
      }
      
      progress <- min(self$iteration / self$max_iterations, 1.0)
      
      self$current_nu <- switch(
        self$nu_schedule,
        
        "linear" = {
          self$nu_min + progress * (self$nu_max - self$nu_min)
        },
        
        "exponential" = {
          self$nu_min * exp(progress * log(self$nu_max / self$nu_min))
        },
        
        "sigmoid" = {
          # Sigmoid with steepness = 10, centered at 0.5
          x <- 10 * (progress - 0.5)
          sigmoid <- 1 / (1 + exp(-x))
          self$nu_min + sigmoid * (self$nu_max - self$nu_min)
        },
        
        "performance" = {
          # Adapt based on improvement rate
          if (length(self$improvement_history) < 5) {
            self$nu_min
          } else {
            recent_improvement <- mean(tail(self$improvement_history, 5))
            
            if (recent_improvement > 0.01) {
              # Good improvement: maintain exploration
              max(self$nu_min, self$current_nu * 0.95)
            } else {
              # Poor improvement: increase exploitation
              min(self$nu_max, self$current_nu * 1.1)
            }
          }
        },
        
        # Default: linear
        self$nu_min + progress * (self$nu_max - self$nu_min)
      )
      
      # Update parent nu
      self$nu <- self$current_nu
      
      invisible(self$current_nu)
    },
    
    #' @description
    #' Add observation and update nu schedule
    #' @inheritParams GaussianBO$add_observation
    add_observation = function(x, y_obs) {
      # Track improvement
      if (length(self$y) > 0) {
        improvement <- max(0, y_obs - max(self$y))
        self$improvement_history <- c(self$improvement_history, improvement)
      }
      
      super$add_observation(x, y_obs)
      self$update_nu()
    },
    
    #' @description
    #' Get adaptation diagnostics
    #' @return List with nu history and improvement metrics
    get_diagnostics = function() {
      list(
        current_nu = self$current_nu,
        iteration = self$iteration,
        nu_schedule = self$nu_schedule,
        improvement_history = self$improvement_history,
        mean_improvement = if (length(self$improvement_history) > 0) {
          mean(self$improvement_history)
        } else {
          NA
        }
      )
    }
  )
)
