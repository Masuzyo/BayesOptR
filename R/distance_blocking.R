#' Distance-Based Hard Blocking Strategy
#'
#' @description
#' Hard blocking based on Euclidean distance from pending/recent evaluations.
#' Creates exclusion spheres around reference points.
#'
#' @details
#' - Calculates Euclidean distance to nearest pending/recent point
#' - Hard exclusion within radius threshold
#' - Effective for preventing tight clustering
#' - Respects full-dimensional geometry
#'
#' @export
DistanceHardBlocking <- R6::R6Class(
  "DistanceHardBlocking",
  
  public = list(
    radius = NULL,
    pending_points = NULL,
    n_recent = NULL,
    recent_points = NULL,
    use_targeted = NULL,
    all_evaluated_points = NULL,
    
    #' @description
    #' Initialize distance-based hard blocking
    #' @param radius Exclusion radius (as fraction of search space diagonal, default: 0.1)
    #' @param n_recent Number of recent points to consider (default: 5)
    #' @param use_targeted Use targeted blocking (default: TRUE)
    initialize = function(radius = 0.1, n_recent = 5, use_targeted = TRUE) {
      self$radius <- radius
      self$n_recent <- n_recent
      self$use_targeted <- use_targeted
      self$pending_points <- matrix(nrow = 0, ncol = 0)
      self$recent_points <- matrix(nrow = 0, ncol = 0)
      self$all_evaluated_points <- NULL
    },
    
    #' @description
    #' Set pending evaluation points
    #' @param points Matrix of pending points
    set_pending = function(points) {
      if (is.null(points) || nrow(points) == 0) {
        self$pending_points <- matrix(nrow = 0, ncol = if(ncol(points) > 0) ncol(points) else 0)
      } else {
        self$pending_points <- points
      }
    },
    
    #' @description
    #' Update recent points history
    #' @param all_points Matrix of all evaluated points
    update_recent = function(all_points) {
      if (is.null(all_points) || nrow(all_points) == 0) {
        return()
      }
      
      n <- nrow(all_points)
      start_idx <- max(1, n - self$n_recent + 1)
      self$recent_points <- all_points[start_idx:n, , drop = FALSE]
      self$all_evaluated_points <- all_points
    },
    
    #' @description
    #' Calculate minimum distance to reference points
    #' @param candidates Matrix of candidate points
    #' @param bounds Matrix of bounds for normalization
    #' @return Vector of minimum distances (normalized)
    calculate_distances = function(candidates, bounds) {
      reference_points <- rbind(self$pending_points, self$recent_points)
      
      if (nrow(reference_points) == 0) {
        return(rep(Inf, nrow(candidates)))
      }
      
      # Normalize to [0, 1] for distance calculation
      n_dim <- ncol(candidates)
      cand_norm <- candidates
      ref_norm <- reference_points
      
      for (d in 1:n_dim) {
        range_d <- bounds[d, 2] - bounds[d, 1]
        cand_norm[, d] <- (candidates[, d] - bounds[d, 1]) / range_d
        ref_norm[, d] <- (reference_points[, d] - bounds[d, 1]) / range_d
      }
      
      # Calculate minimum distance for each candidate
      min_distances <- apply(cand_norm, 1, function(x) {
        distances <- sqrt(rowSums((sweep(ref_norm, 2, x))^2))
        min(distances)
      })
      
      min_distances
    },
    
    #' @description
    #' Apply hard distance-based blocking with targeted penalization
    #' @param acquisition_values Vector of acquisition values
    #' @param candidates Matrix of candidate points
    #' @param bounds Matrix of bounds
    #' @param current_best Optional current best point (for targeted blocking)
    #' @return Modified acquisition values (blocked points get -Inf)
    apply_blocking = function(acquisition_values, candidates, bounds, current_best = NULL) {
      # If using targeted blocking and have evaluated points
      if (self$use_targeted && !is.null(self$all_evaluated_points) && 
          nrow(self$all_evaluated_points) > 0 && !is.null(current_best)) {
        
        n_candidates <- nrow(candidates)
        n_dim <- ncol(candidates)
        
        # f^+ is the current best point
        f_plus <- current_best
        
        # Find closest known point to f^+
        f_plus_norm <- f_plus
        points_norm <- self$all_evaluated_points
        
        for (d in 1:n_dim) {
          range_d <- bounds[d, 2] - bounds[d, 1]
          f_plus_norm[d] <- (f_plus[d] - bounds[d, 1]) / range_d
          points_norm[, d] <- (self$all_evaluated_points[, d] - bounds[d, 1]) / range_d
        }
        
        distances_to_target <- sqrt(rowSums((sweep(points_norm, 2, f_plus_norm))^2))
        g_min <- self$all_evaluated_points[which.min(distances_to_target), ]
        
        # Normalize all points
        g_min_norm <- g_min
        candidates_norm <- candidates
        
        for (d in 1:n_dim) {
          range_d <- bounds[d, 2] - bounds[d, 1]
          g_min_norm[d] <- (g_min[d] - bounds[d, 1]) / range_d
          candidates_norm[, d] <- (candidates[, d] - bounds[d, 1]) / range_d
        }
        
        dist_target_to_closest <- sqrt(sum((f_plus_norm - g_min_norm)^2))
        
        for (i in 1:n_candidates) {
          f <- candidates_norm[i, ]
          dist_target_to_f <- sqrt(sum((f_plus_norm - f)^2))
          dist_closest_to_f <- sqrt(sum((g_min_norm - f)^2))
          
          # Check if in blocking region
          if (dist_target_to_f < dist_target_to_closest && 
              dist_closest_to_f < dist_target_to_closest) {
            # In blocking region - apply hard block
            acquisition_values[i] <- -Inf
          }
        }
        
        return(acquisition_values)
      }
      
      # Fallback to standard distance blocking
      distances <- self$calculate_distances(candidates, bounds)
      blocked <- distances < self$radius
      acquisition_values[blocked] <- -Inf
      
      acquisition_values
    }
  )
)


#' Distance-Based Soft Blocking Strategy
#'
#' @description
#' Soft blocking with smooth penalty based on distance from reference points.
#'
#' @details
#' - Penalty decreases smoothly with distance
#' - Uses exponential decay or smooth kernel
#' - More gradual than hard blocking
#' - Allows controlled exploration near pending points
#'
#' @export
DistanceSoftBlocking <- R6::R6Class(
  "DistanceSoftBlocking",
  inherit = DistanceHardBlocking,
  
  public = list(
    penalty_strength = NULL,
    kernel = NULL,
    
    #' @description
    #' Initialize distance-based soft blocking
    #' @param penalty_strength Strength of distance penalty (default: 3.0)
    #' @param kernel Kernel type: "exponential" or "gaussian" (default: "exponential")
    #' @inheritParams DistanceHardBlocking
    initialize = function(radius = 0.2, n_recent = 5, 
                         penalty_strength = 3.0, kernel = "exponential") {
      super$initialize(radius, n_recent)
      self$penalty_strength <- penalty_strength
      self$kernel <- kernel
    },
    
    #' @description
    #' Apply soft distance-based penalty with targeted penalization
    #' @inheritParams DistanceHardBlocking$apply_blocking
    apply_blocking = function(acquisition_values, candidates, bounds, current_best = NULL) {
      # If using targeted blocking and have evaluated points
      if (self$use_targeted && !is.null(self$all_evaluated_points) && 
          nrow(self$all_evaluated_points) > 0 && !is.null(current_best)) {
        
        n_candidates <- nrow(candidates)
        n_dim <- ncol(candidates)
        
        # f^+ is the current best point
        f_plus <- current_best
        
        # Find closest known point to f^+
        f_plus_norm <- f_plus
        points_norm <- self$all_evaluated_points
        
        for (d in 1:n_dim) {
          range_d <- bounds[d, 2] - bounds[d, 1]
          f_plus_norm[d] <- (f_plus[d] - bounds[d, 1]) / range_d
          points_norm[, d] <- (self$all_evaluated_points[, d] - bounds[d, 1]) / range_d
        }
        
        distances_to_target <- sqrt(rowSums((sweep(points_norm, 2, f_plus_norm))^2))
        g_min <- self$all_evaluated_points[which.min(distances_to_target), ]
        
        # Normalize all points
        g_min_norm <- g_min
        candidates_norm <- candidates
        
        for (d in 1:n_dim) {
          range_d <- bounds[d, 2] - bounds[d, 1]
          g_min_norm[d] <- (g_min[d] - bounds[d, 1]) / range_d
          candidates_norm[, d] <- (candidates[, d] - bounds[d, 1]) / range_d
        }
        
        dist_target_to_closest <- sqrt(sum((f_plus_norm - g_min_norm)^2))
        
        for (i in 1:n_candidates) {
          f <- candidates_norm[i, ]
          dist_target_to_f <- sqrt(sum((f_plus_norm - f)^2))
          dist_closest_to_f <- sqrt(sum((g_min_norm - f)^2))
          
          # Check if in blocking region
          if (dist_target_to_f < dist_target_to_closest && 
              dist_closest_to_f < dist_target_to_closest) {
            # Penalty = |f - g_min| / |f^+ - g_min|
            penalty <- dist_closest_to_f / (dist_target_to_closest + 1e-10)
            acquisition_values[i] <- acquisition_values[i] * penalty
          }
        }
        
        return(acquisition_values)
      }
      
      # Fallback to standard distance blocking
      distances <- self$calculate_distances(candidates, bounds)
      
      penalties <- switch(
        self$kernel,
        "exponential" = {
          exp(-self$penalty_strength * distances / self$radius)
        },
        "gaussian" = {
          exp(-self$penalty_strength * (distances / self$radius)^2)
        },
        exp(-self$penalty_strength * distances / self$radius)
      )
      
      acquisition_values * (1 - penalties)
    }
  )
)


#' Targeted Blocking Strategy
#'
#' @description
#' Blocking based on the region between selected target (f^+) and closest known point (g_min).
#' Penalization is applied only in the region between these two points.
#'
#' @details
#' Penalization formula:
#' - penalty = |f - g_min| / |f^+ - g_min| if |f^+ - f| < |f^+ - g_min|
#' - penalty = 1 otherwise
#' 
#' This creates a focused blocking region between the target and nearest point,
#' rather than blocking the entire space.
#'
#' @export
TargetedBlocking <- R6::R6Class(
  "TargetedBlocking",
  
  public = list(
    pending_points = NULL,
    target_point = NULL,
    
    #' @description
    #' Initialize targeted blocking
    initialize = function() {
      self$pending_points <- matrix(nrow = 0, ncol = 0)
      self$target_point <- NULL
    },
    
    #' @description
    #' Set pending evaluation points
    #' @param points Matrix of pending points
    set_pending = function(points) {
      if (is.null(points) || nrow(points) == 0) {
        self$pending_points <- matrix(nrow = 0, ncol = if(!is.null(points)) ncol(points) else 0)
      } else {
        self$pending_points <- points
      }
    },
    
    #' @description
    #' Set the target point (f^+) for current batch
    #' @param point Vector representing the target point
    set_target = function(point) {
      self$target_point <- point
    },
    
    #' @description
    #' Find closest known point (g_min) to target
    #' @param all_points Matrix of all evaluated points
    #' @param bounds Matrix of bounds for normalization
    #' @return Index of closest point
    find_closest_to_target = function(all_points, bounds) {
      if (is.null(self$target_point) || nrow(all_points) == 0) {
        return(NULL)
      }
      
      n_dim <- length(self$target_point)
      
      # Normalize points
      target_norm <- self$target_point
      points_norm <- all_points
      
      for (d in 1:n_dim) {
        range_d <- bounds[d, 2] - bounds[d, 1]
        target_norm[d] <- (self$target_point[d] - bounds[d, 1]) / range_d
        points_norm[, d] <- (all_points[, d] - bounds[d, 1]) / range_d
      }
      
      # Calculate distances to target
      distances <- sqrt(rowSums((sweep(points_norm, 2, target_norm))^2))
      which.min(distances)
    },
    
    #' @description
    #' Calculate penalization based on position relative to f^+ and g_min
    #' @param candidates Matrix of candidate points
    #' @param all_points Matrix of all evaluated points
    #' @param bounds Matrix of bounds for normalization
    #' @return Vector of penalization factors (0 to 1)
    calculate_penalization = function(candidates, all_points, bounds) {
      if (is.null(self$target_point) || nrow(all_points) == 0) {
        return(rep(1.0, nrow(candidates)))
      }
      
      n_candidates <- nrow(candidates)
      n_dim <- ncol(candidates)
      
      # Find closest known point (g_min)
      g_min_idx <- self$find_closest_to_target(all_points, bounds)
      if (is.null(g_min_idx)) {
        return(rep(1.0, n_candidates))
      }
      
      g_min <- all_points[g_min_idx, ]
      f_plus <- self$target_point
      
      # Normalize all points
      candidates_norm <- candidates
      g_min_norm <- g_min
      f_plus_norm <- f_plus
      
      for (d in 1:n_dim) {
        range_d <- bounds[d, 2] - bounds[d, 1]
        candidates_norm[, d] <- (candidates[, d] - bounds[d, 1]) / range_d
        g_min_norm[d] <- (g_min[d] - bounds[d, 1]) / range_d
        f_plus_norm[d] <- (f_plus[d] - bounds[d, 1]) / range_d
      }
      
      # Distance from f^+ to g_min
      dist_target_to_closest <- sqrt(sum((f_plus_norm - g_min_norm)^2))
      
      # Calculate penalization for each candidate
      penalties <- numeric(n_candidates)
      
      for (i in 1:n_candidates) {
        f <- candidates_norm[i, ]
        
        # Distance from f^+ to candidate f
        dist_target_to_f <- sqrt(sum((f_plus_norm - f)^2))
        
        # Distance from g_min to candidate f
        dist_closest_to_f <- sqrt(sum((g_min_norm - f)^2))
        
        # Check if f is in the blocking region
        # Condition 1: |f^+ - f| < |f^+ - g_min| (f is within sphere around f^+)
        # Condition 2: |f - g_min| < |f^+ - g_min| (f is between g_min and f^+)
        if (dist_target_to_f < dist_target_to_closest && 
            dist_closest_to_f < dist_target_to_closest) {
          # Penalization = |f - g_min| / |f^+ - g_min|
          penalties[i] <- dist_closest_to_f / (dist_target_to_closest + 1e-10)
        } else {
          # Outside blocking region
          penalties[i] <- 1.0
        }
      }
      
      penalties
    },
    
    #' @description
    #' Apply targeted blocking to acquisition values
    #' @param acquisition_values Vector of acquisition values
    #' @param candidates Matrix of candidate points
    #' @param all_points Matrix of all evaluated points
    #' @param bounds Matrix of bounds
    #' @return Penalized acquisition values
    apply_blocking = function(acquisition_values, candidates, all_points, bounds) {
      penalties <- self$calculate_penalization(candidates, all_points, bounds)
      acquisition_values * penalties
    }
  )
)
