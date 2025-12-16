#' Dimension-Wise Hard Blocking Strategy
#'
#' @description
#' Hard blocking strategy that divides search space into dimension-wise blocks
#' and completely excludes regions around pending/recent evaluations.
#'
#' @details
#' For each dimension independently:
#' - Identifies blocks containing pending/recent points
#' - Creates exclusion zones with hard boundaries
#' - Forces exploration in unblocked regions
#' - Useful for preventing dimension-wise clustering
#'
#' @export
DimensionWiseHardBlocking <- R6::R6Class(
  "DimensionWiseHardBlocking",
  
  public = list(
    block_size = NULL,
    pending_points = NULL,
    n_recent = NULL,
    recent_points = NULL,
    use_targeted = NULL,
    all_evaluated_points = NULL,
    
    #' @description
    #' Initialize dimension-wise hard blocking
    #' @param block_size Number of blocks per dimension (default: 5)
    #' @param n_recent Number of recent points to consider (default: 5)
    #' @param use_targeted Use targeted blocking (default: TRUE)
    initialize = function(block_size = 5, n_recent = 5, use_targeted = TRUE) {
      self$block_size <- block_size
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
        self$pending_points <- matrix(nrow = 0, ncol = ncol(points))
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
    #' Get blocked dimensions for each candidate point
    #' @param candidates Matrix of candidate points
    #' @param bounds Matrix of bounds (n_dim × 2)
    #' @return Matrix indicating blocked dimensions (1=blocked, 0=free)
    get_blocked_dimensions = function(candidates, bounds) {
      n_candidates <- nrow(candidates)
      n_dim <- ncol(candidates)
      
      blocked <- matrix(0, nrow = n_candidates, ncol = n_dim)
      
      # Combine pending and recent points
      reference_points <- rbind(self$pending_points, self$recent_points)
      
      if (nrow(reference_points) == 0) {
        return(blocked)
      }
      
      # For each dimension
      for (d in 1:n_dim) {
        # Get dimension range
        dim_min <- bounds[d, 1]
        dim_max <- bounds[d, 2]
        dim_range <- dim_max - dim_min
        
        # Create blocks
        block_width <- dim_range / self$block_size
        
        # Find which blocks contain reference points
        ref_blocks <- floor((reference_points[, d] - dim_min) / block_width)
        ref_blocks <- unique(pmin(ref_blocks, self$block_size - 1))
        
        # Check each candidate
        cand_blocks <- floor((candidates[, d] - dim_min) / block_width)
        cand_blocks <- pmin(cand_blocks, self$block_size - 1)
        
        # Mark as blocked if in same block as reference
        for (i in 1:n_candidates) {
          if (cand_blocks[i] %in% ref_blocks) {
            blocked[i, d] <- 1
          }
        }
      }
      
      blocked
    },
    
    #' @description
    #' Apply blocking penalty to acquisition values with targeted penalization
    #' @param acquisition_values Vector of acquisition function values
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
      
      # Fallback to standard dimension-wise blocking
      blocked <- self$get_blocked_dimensions(candidates, bounds)
      
      # Hard blocking: any blocked dimension → -Inf
      any_blocked <- rowSums(blocked) > 0
      acquisition_values[any_blocked] <- -Inf
      
      acquisition_values
    }
  )
)


#' Dimension-Wise Soft Blocking Strategy
#'
#' @description
#' Soft blocking that penalizes (but doesn't completely exclude) regions
#' in dimension-wise blocks around pending/recent evaluations.
#'
#' @details
#' Similar to hard blocking but uses penalty weights instead of hard exclusion:
#' - Penalty decreases with distance in each dimension
#' - Allows some exploration in penalized regions
#' - More flexible than hard blocking
#'
#' @export
DimensionWiseSoftBlocking <- R6::R6Class(
  "DimensionWiseSoftBlocking",
  inherit = DimensionWiseHardBlocking,
  
  public = list(
    penalty_strength = NULL,
    
    #' @description
    #' Initialize dimension-wise soft blocking
    #' @param penalty_strength Penalty multiplier (default: 0.5)
    #' @inheritParams DimensionWiseHardBlocking
    initialize = function(block_size = 5, n_recent = 5, penalty_strength = 0.5, use_targeted = TRUE) {
      super$initialize(block_size, n_recent, use_targeted)
      self$penalty_strength <- penalty_strength
    },
    
    #' @description
    #' Apply soft blocking penalty with targeted penalization
    #' @inheritParams DimensionWiseHardBlocking$apply_blocking
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
      
      # Fallback to standard dimension-wise blocking
      n_candidates <- nrow(candidates)
      n_dim <- ncol(candidates)
      
      reference_points <- rbind(self$pending_points, self$recent_points)
      
      if (nrow(reference_points) == 0) {
        return(acquisition_values)
      }
      
      penalties <- rep(1.0, n_candidates)
      
      # For each dimension
      for (d in 1:n_dim) {
        dim_min <- bounds[d, 1]
        dim_max <- bounds[d, 2]
        dim_range <- dim_max - dim_min
        block_width <- dim_range / self$block_size
        
        # Find blocks with reference points
        ref_blocks <- floor((reference_points[, d] - dim_min) / block_width)
        ref_blocks <- unique(pmin(ref_blocks, self$block_size - 1))
        
        # Calculate candidate blocks
        cand_blocks <- floor((candidates[, d] - dim_min) / block_width)
        cand_blocks <- pmin(cand_blocks, self$block_size - 1)
        
        # Apply penalty based on block distance
        for (i in 1:n_candidates) {
          min_block_dist <- min(abs(cand_blocks[i] - ref_blocks))
          
          # Exponential decay penalty
          dim_penalty <- exp(-self$penalty_strength * min_block_dist)
          penalties[i] <- penalties[i] * (1 - dim_penalty)
        }
      }
      
      acquisition_values * penalties
    }
  )
)
