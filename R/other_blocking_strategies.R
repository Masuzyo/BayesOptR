#' Local Penalization Strategy
#'
#' @description
#' Local penalization that modifies acquisition function based on distance
#' to pending evaluations. Implements the Gonzalez et al. (2016) approach.
#'
#' @details
#' Key features:
#' - Lipschitz-based penalization
#' - Considers predictive variance and distance
#' - Balances exploration and avoiding pending regions
#' - Theoretically grounded in GP-UCB framework
#'
#' Reference: González, J., et al. (2016). Batch Bayesian Optimization via 
#' Local Penalization. AISTATS.
#'
#' @export
LocalPenalization <- R6::R6Class(
  "LocalPenalization",
  
  public = list(
    lipschitz_constant = NULL,
    pending_points = NULL,
    kappa = NULL,
    
    #' @description
    #' Initialize local penalization
    #' @param lipschitz_constant Lipschitz constant estimate (default: 10.0)
    #' @param kappa Exploration parameter (default: 2.0)
    initialize = function(lipschitz_constant = 10.0, kappa = 2.0) {
      self$lipschitz_constant <- lipschitz_constant
      self$kappa <- kappa
      self$pending_points <- matrix(nrow = 0, ncol = 0)
    },
    
    #' @description
    #' Set pending evaluation points
    #' @param points Matrix of pending points
    set_pending = function(points) {
      if (is.null(points) || nrow(points) == 0) {
        self$pending_points <- matrix(nrow = 0, ncol = if(length(points) > 0) ncol(points) else 0)
      } else {
        self$pending_points <- points
      }
    },
    
    #' @description
    #' Calculate penalization factor
    #' @param candidates Matrix of candidate points
    #' @param mean_pred Vector of predictive means at candidates
    #' @param std_pred Vector of predictive standard deviations
    #' @param bounds Matrix of bounds for normalization
    #' @return Vector of penalization factors (0 to 1)
    calculate_penalization = function(candidates, mean_pred, std_pred, bounds) {
      if (nrow(self$pending_points) == 0) {
        return(rep(1.0, nrow(candidates)))
      }
      
      n_candidates <- nrow(candidates)
      n_pending <- nrow(self$pending_points)
      n_dim <- ncol(candidates)
      
      # Normalize points
      cand_norm <- candidates
      pend_norm <- self$pending_points
      
      for (d in 1:n_dim) {
        range_d <- bounds[d, 2] - bounds[d, 1]
        cand_norm[, d] <- (candidates[, d] - bounds[d, 1]) / range_d
        pend_norm[, d] <- (self$pending_points[, d] - bounds[d, 1]) / range_d
      }
      
      penalization <- rep(1.0, n_candidates)
      
      for (i in 1:n_candidates) {
        # Calculate distances to all pending points
        distances <- sqrt(rowSums((sweep(pend_norm, 2, cand_norm[i, ]))^2))
        
        # Local penalization formula
        # r(x) = min_{x' in pending} [μ(x') + κ·σ(x') - L·||x - x'||]
        for (j in 1:n_pending) {
          # Upper confidence bound at pending point (approximate)
          ucb_pending <- mean_pred[i] + self$kappa * std_pred[i]
          
          # Penalization based on distance
          penalty <- ucb_pending - self$lipschitz_constant * distances[j]
          
          # Apply penalty if positive
          if (penalty > 0) {
            penalization[i] <- min(penalization[i], 
                                  1.0 - penalty / (ucb_pending + 1e-8))
          }
        }
      }
      
      pmax(penalization, 0)  # Ensure non-negative
    },
    
    #' @description
    #' Apply local penalization to acquisition values
    #' @param acquisition_values Vector of acquisition values
    #' @param candidates Matrix of candidate points
    #' @param mean_pred Vector of predictive means
    #' @param std_pred Vector of predictive std deviations
    #' @param bounds Matrix of bounds
    #' @return Penalized acquisition values
    apply_blocking = function(acquisition_values, candidates, mean_pred, 
                             std_pred, bounds) {
      penalization <- self$calculate_penalization(candidates, mean_pred, 
                                                  std_pred, bounds)
      acquisition_values * penalization
    }
  )
)


#' Constant Liar Strategy
#'
#' @description
#' Heuristic batch acquisition using "constant liar" fantasy points.
#' Temporarily assigns a constant value to pending evaluations.
#'
#' @details
#' Strategy variants:
#' - "min": Pessimistic (assume pending → worst observed)
#' - "mean": Neutral (assume pending → mean observed)
#' - "max": Optimistic (assume pending → best observed)
#' - "kriging": Use Kriging believer (predictive mean)
#'
#' @export
ConstantLiar <- R6::R6Class(
  "ConstantLiar",
  
  public = list(
    strategy = NULL,
    pending_points = NULL,
    fantasy_values = NULL,
    
    #' @description
    #' Initialize Constant Liar
    #' @param strategy Lying strategy: "min", "mean", "max", "kriging" (default: "mean")
    initialize = function(strategy = "mean") {
      if (!strategy %in% c("min", "mean", "max", "kriging")) {
        stop("strategy must be one of: 'min', 'mean', 'max', 'kriging'")
      }
      
      self$strategy <- strategy
      self$pending_points <- matrix(nrow = 0, ncol = 0)
      self$fantasy_values <- numeric(0)
    },
    
    #' @description
    #' Generate fantasy values for pending points
    #' @param pending Matrix of pending points
    #' @param observed_y Vector of observed function values
    #' @param model BO model object (for kriging strategy)
    #' @return Vector of fantasy values
    generate_fantasy_values = function(pending, observed_y, model = NULL) {
      if (nrow(pending) == 0) {
        return(numeric(0))
      }
      
      n_pending <- nrow(pending)
      
      fantasy <- switch(
        self$strategy,
        
        "min" = rep(min(observed_y), n_pending),
        
        "mean" = rep(mean(observed_y), n_pending),
        
        "max" = rep(max(observed_y), n_pending),
        
        "kriging" = {
          if (is.null(model)) {
            # Fallback to mean if no model
            rep(mean(observed_y), n_pending)
          } else {
            # Use predictive mean
            pred <- model$predict(pending)
            pred$mean
          }
        },
        
        # Default: mean
        rep(mean(observed_y), n_pending)
      )
      
      fantasy
    },
    
    #' @description
    #' Set pending points and compute fantasy values
    #' @param pending Matrix of pending points
    #' @param observed_y Vector of observed values
    #' @param model BO model (optional, for kriging)
    set_pending = function(pending, observed_y, model = NULL) {
      self$pending_points <- pending
      self$fantasy_values <- self$generate_fantasy_values(pending, observed_y, model)
    },
    
    #' @description
    #' Get augmented dataset with fantasy points
    #' @param X Matrix of observed inputs
    #' @param y Vector of observed outputs
    #' @return List with X_augmented and y_augmented
    get_augmented_data = function(X, y) {
      if (nrow(self$pending_points) == 0) {
        return(list(X = X, y = y))
      }
      
      list(
        X = rbind(X, self$pending_points),
        y = c(y, self$fantasy_values)
      )
    }
  )
)
