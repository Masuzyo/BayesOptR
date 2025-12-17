#' Diverse Batch Selection
#'
#' @description
#' Select a diverse batch of points that balances diversity in observations 
#' with score optimization. Useful for selecting representative subsets from 
#' large candidate pools or creating diverse training batches.
#'
#' @param X Matrix of observations (n × d), where n is number of points and d is dimensions
#' @param scores Vector of scores/values (length n)
#' @param batch_size Number of points to select
#' @param minimize Logical; if TRUE, minimize scores, if FALSE, maximize scores (default: FALSE)
#' @param diversity_weight Weight for diversity vs score quality (0-1, default: 0.5)
#' @param method Selection method: "greedy" (default), "kmeans", or "farthest"
#' @param return_indices If TRUE, return indices; if FALSE, return actual points (default: TRUE)
#'
#' @return If return_indices=TRUE, vector of selected indices; 
#'         if FALSE, list with selected points (X) and scores (y)
#'
#' @details
#' Selection Methods:
#' \itemize{
#'   \item \strong{greedy}: Sequentially select points maximizing diversity + score quality
#'   \item \strong{kmeans}: K-means clustering then select best from each cluster
#'   \item \strong{farthest}: Farthest point sampling biased by scores
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(100 * 3), ncol = 3)
#' scores <- apply(X, 1, function(x) -sum(x^2))
#'
#' # Select 10 diverse high-scoring points
#' indices <- select_diverse_batch(X, scores, batch_size = 10, minimize = FALSE)
#' selected_X <- X[indices, ]
#'
#' # With more emphasis on diversity
#' indices_div <- select_diverse_batch(X, scores, batch_size = 10, 
#'                                      diversity_weight = 0.8)
#' }
#'
#' @export
select_diverse_batch <- function(X, scores, batch_size, 
                                 minimize = FALSE,
                                 diversity_weight = 0.5,
                                 method = "greedy",
                                 return_indices = TRUE) {
  
  # Input validation
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame")
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (length(scores) != n) {
    stop("Length of scores must equal number of rows in X")
  }
  
  if (batch_size < 1 || batch_size > n) {
    stop("batch_size must be between 1 and nrow(X)")
  }
  
  if (diversity_weight < 0 || diversity_weight > 1) {
    stop("diversity_weight must be between 0 and 1")
  }
  
  # Normalize X to [0, 1]^d for distance computation
  X_norm <- X
  for (i in 1:d) {
    x_min <- min(X[, i])
    x_max <- max(X[, i])
    if (x_max > x_min) {
      X_norm[, i] <- (X[, i] - x_min) / (x_max - x_min)
    }
  }
  
  # Normalize scores to [0, 1]
  scores_norm <- scores
  if (minimize) {
    scores_norm <- -scores_norm
  }
  score_min <- min(scores_norm)
  score_max <- max(scores_norm)
  if (score_max > score_min) {
    scores_norm <- (scores_norm - score_min) / (score_max - score_min)
  } else {
    scores_norm <- rep(0.5, n)
  }
  
  # Select method
  selected_indices <- switch(
    method,
    "greedy" = select_greedy(X_norm, scores_norm, batch_size, diversity_weight),
    "kmeans" = select_kmeans(X_norm, scores_norm, batch_size, minimize),
    "farthest" = select_farthest(X_norm, scores_norm, batch_size, diversity_weight),
    stop("Unknown method: ", method, ". Use 'greedy', 'kmeans', or 'farthest'")
  )
  
  # Return
  if (return_indices) {
    return(selected_indices)
  } else {
    return(list(
      X = X[selected_indices, , drop = FALSE],
      scores = scores[selected_indices]
    ))
  }
}


#' Greedy Diverse Batch Selection
#'
#' @description
#' Greedily select points that maximize a combination of diversity and score quality.
#'
#' @keywords internal
select_greedy <- function(X_norm, scores_norm, batch_size, diversity_weight) {
  n <- nrow(X_norm)
  selected <- integer(batch_size)
  available <- 1:n
  
  # Start with the best scoring point
  selected[1] <- which.max(scores_norm)
  available <- available[-selected[1]]
  
  # Greedily add points
  for (i in 2:batch_size) {
    if (length(available) == 0) break
    
    # Compute diversity scores (minimum distance to selected points)
    diversity_scores <- sapply(available, function(idx) {
      dists <- sqrt(rowSums((X_norm[selected[1:(i-1)], , drop = FALSE] - 
                             matrix(X_norm[idx, ], nrow = i-1, ncol = ncol(X_norm), byrow = TRUE))^2))
      min(dists)
    })
    
    # Normalize diversity scores
    if (max(diversity_scores) > min(diversity_scores)) {
      diversity_scores <- (diversity_scores - min(diversity_scores)) / 
                          (max(diversity_scores) - min(diversity_scores))
    } else {
      diversity_scores <- rep(0.5, length(diversity_scores))
    }
    
    # Combine diversity and score quality
    quality_scores <- scores_norm[available]
    combined_scores <- diversity_weight * diversity_scores + 
                       (1 - diversity_weight) * quality_scores
    
    # Select best
    best_idx <- which.max(combined_scores)
    selected[i] <- available[best_idx]
    available <- available[-best_idx]
  }
  
  return(selected)
}


#' K-means Based Batch Selection
#'
#' @description
#' Use k-means clustering to ensure spatial diversity, then select best from each cluster.
#'
#' @keywords internal
select_kmeans <- function(X_norm, scores_norm, batch_size, minimize) {
  n <- nrow(X_norm)
  
  if (batch_size >= n) {
    return(1:n)
  }
  
  # Handle case where batch_size = 1
  if (batch_size == 1) {
    return(which.max(scores_norm))
  }
  
  # K-means clustering
  set.seed(123)  # For reproducibility
  kmeans_result <- kmeans(X_norm, centers = batch_size, nstart = 10)
  clusters <- kmeans_result$cluster
  
  # Select best point from each cluster
  selected <- integer(batch_size)
  for (k in 1:batch_size) {
    cluster_indices <- which(clusters == k)
    if (length(cluster_indices) > 0) {
      # Select point with best score in this cluster
      selected[k] <- cluster_indices[which.max(scores_norm[cluster_indices])]
    } else {
      # Fallback if cluster is empty (shouldn't happen with good k-means)
      remaining <- setdiff(1:n, selected[1:(k-1)])
      selected[k] <- remaining[which.max(scores_norm[remaining])]
    }
  }
  
  return(selected)
}


#' Farthest Point Sampling with Score Bias
#'
#' @description
#' Iteratively select points that are farthest from already selected points,
#' with bias towards better scores.
#'
#' @keywords internal
select_farthest <- function(X_norm, scores_norm, batch_size, diversity_weight) {
  n <- nrow(X_norm)
  selected <- integer(batch_size)
  
  # Start with best scoring point
  selected[1] <- which.max(scores_norm)
  
  # Iteratively select farthest points
  for (i in 2:batch_size) {
    # Compute minimum distance to selected points for all remaining points
    min_dists <- sapply(1:n, function(idx) {
      if (idx %in% selected[1:(i-1)]) {
        return(0)  # Already selected
      }
      dists <- sqrt(rowSums((X_norm[selected[1:(i-1)], , drop = FALSE] - 
                            matrix(X_norm[idx, ], nrow = i-1, ncol = ncol(X_norm), byrow = TRUE))^2))
      min(dists)
    })
    
    # Normalize distances
    if (max(min_dists) > 0) {
      min_dists_norm <- min_dists / max(min_dists)
    } else {
      min_dists_norm <- min_dists
    }
    
    # Combine distance with score quality
    combined <- diversity_weight * min_dists_norm + (1 - diversity_weight) * scores_norm
    
    # Zero out already selected points
    combined[selected[1:(i-1)]] <- -Inf
    
    # Select point with maximum combined score
    selected[i] <- which.max(combined)
  }
  
  return(selected)
}


#' Maximum Diversity Batch Selection
#'
#' @description
#' Select batch that maximizes minimum pairwise distance (max-min diversity).
#' Pure diversity optimization without considering scores.
#'
#' @param X Matrix of observations (n × d)
#' @param batch_size Number of points to select
#' @param return_indices If TRUE, return indices; if FALSE, return actual points
#'
#' @return Vector of selected indices or matrix of selected points
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 2), ncol = 2)
#' indices <- select_max_diversity_batch(X, batch_size = 10)
#' }
#'
#' @export
select_max_diversity_batch <- function(X, batch_size, return_indices = TRUE) {
  
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame")
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (batch_size < 1 || batch_size > n) {
    stop("batch_size must be between 1 and nrow(X)")
  }
  
  # Normalize X
  X_norm <- X
  for (i in 1:d) {
    x_min <- min(X[, i])
    x_max <- max(X[, i])
    if (x_max > x_min) {
      X_norm[, i] <- (X[, i] - x_min) / (x_max - x_min)
    }
  }
  
  # Farthest point sampling
  selected <- integer(batch_size)
  
  # Start with a random point
  selected[1] <- sample(1:n, 1)
  
  # Iteratively select farthest points
  for (i in 2:batch_size) {
    # Compute minimum distance to selected points
    min_dists <- sapply(1:n, function(idx) {
      if (idx %in% selected[1:(i-1)]) {
        return(0)
      }
      dists <- sqrt(rowSums((X_norm[selected[1:(i-1)], , drop = FALSE] - 
                            matrix(X_norm[idx, ], nrow = i-1, ncol = d, byrow = TRUE))^2))
      min(dists)
    })
    
    # Select point with maximum minimum distance
    selected[i] <- which.max(min_dists)
  }
  
  # Return
  if (return_indices) {
    return(selected)
  } else {
    return(X[selected, , drop = FALSE])
  }
}


#' Score-Weighted Sampling
#'
#' @description
#' Sample a batch with probability proportional to scores (softmax sampling).
#' Useful for exploration with bias towards better scores.
#'
#' @param X Matrix of observations (n × d)
#' @param scores Vector of scores/values (length n)
#' @param batch_size Number of points to select
#' @param minimize Logical; if TRUE, minimize scores, if FALSE, maximize scores
#' @param temperature Temperature parameter for softmax (default: 1.0)
#' @param return_indices If TRUE, return indices; if FALSE, return actual points
#'
#' @return Vector of selected indices or list with X and scores
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 2), ncol = 2)
#' scores <- apply(X, 1, function(x) -sum(x^2))
#' indices <- select_score_weighted_batch(X, scores, batch_size = 10, 
#'                                         temperature = 0.5)
#' }
#'
#' @export
select_score_weighted_batch <- function(X, scores, batch_size, 
                                       minimize = FALSE,
                                       temperature = 1.0,
                                       return_indices = TRUE) {
  
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame")
  }
  
  n <- nrow(X)
  if (length(scores) != n) {
    stop("Length of scores must equal number of rows in X")
  }
  
  if (batch_size < 1 || batch_size > n) {
    stop("batch_size must be between 1 and nrow(X)")
  }
  
  # Adjust scores for minimization
  if (minimize) {
    scores <- -scores
  }
  
  # Softmax probabilities with temperature
  scores_scaled <- scores / temperature
  scores_exp <- exp(scores_scaled - max(scores_scaled))  # Numerical stability
  probs <- scores_exp / sum(scores_exp)
  
  # Sample without replacement
  selected_indices <- sample(1:n, size = batch_size, replace = FALSE, prob = probs)
  
  # Return
  if (return_indices) {
    return(selected_indices)
  } else {
    return(list(
      X = X[selected_indices, , drop = FALSE],
      scores = scores[selected_indices]
    ))
  }
}
