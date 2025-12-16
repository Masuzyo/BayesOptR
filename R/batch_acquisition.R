#' Batch Acquisition with Blocking Strategies
#'
#' @description
#' Sequential greedy batch acquisition with blocking strategies to ensure diversity.
#'
#' @export
BatchAcquisition <- R6::R6Class(
  "BatchAcquisition",
  
  public = list(
    model = NULL,
    acquisition_fn = NULL,
    blocking_strategy = NULL,
    batch_size = NULL,
    n_candidates = NULL,
    
    #' @description
    #' Initialize batch acquisition
    #' @param model BO model object (GaussianBO, StudentTBO, or AdaptiveHybridBO)
    #' @param acquisition_fn Acquisition function object
    #' @param blocking_strategy Blocking strategy object (optional)
    #' @param batch_size Number of points per batch (default: 4)
    #' @param n_candidates Number of candidate points to evaluate (default: 1000)
    initialize = function(model, acquisition_fn, blocking_strategy = NULL,
                         batch_size = 4, n_candidates = 1000) {
      self$model <- model
      self$acquisition_fn <- acquisition_fn
      self$blocking_strategy <- blocking_strategy
      self$batch_size <- batch_size
      self$n_candidates <- n_candidates
    },
    
    #' @description
    #' Generate candidate points using Sobol sequence
    #' @return Matrix of candidate points
    generate_candidates = function() {
      n_dim <- nrow(self$model$bounds)
      
      # Use random sampling (can be replaced with Sobol/Halton)
      candidates <- matrix(runif(self$n_candidates * n_dim), 
                          nrow = self$n_candidates, 
                          ncol = n_dim)
      
      # Scale to bounds
      for (d in 1:n_dim) {
        candidates[, d] <- self$model$bounds[d, 1] + 
          candidates[, d] * (self$model$bounds[d, 2] - self$model$bounds[d, 1])
      }
      
      candidates
    },
    
    #' @description
    #' Select next batch of points using greedy sequential selection
    #' @return Matrix of selected points (batch_size × n_dim)
    select_batch = function() {
      batch <- matrix(nrow = self$batch_size, ncol = nrow(self$model$bounds))
      pending <- matrix(nrow = 0, ncol = nrow(self$model$bounds))
      
      for (i in 1:self$batch_size) {
        # Generate candidates
        candidates <- self$generate_candidates()
        
        # Get predictions
        pred <- self$model$predict(candidates)
        
        # Evaluate acquisition function
        f_best <- self$model$get_best()$best_y
        
        # Get nu parameter if using Student-t models
        nu_param <- if (!is.null(pred$nu)) pred$nu else NULL
        
        if (inherits(self$acquisition_fn, "ThompsonSampling")) {
          acq_values <- self$acquisition_fn$evaluate(pred$mean, pred$sd)
        } else if (inherits(self$acquisition_fn, "UpperConfidenceBound")) {
          acq_values <- self$acquisition_fn$evaluate(pred$mean, pred$sd, 
                                                     nrow(self$model$bounds))
        } else if (inherits(self$acquisition_fn, "StudentTExpectedImprovement")) {
          # Pass nu parameter for Student-t EI
          acq_values <- self$acquisition_fn$evaluate(pred$mean, pred$sd, f_best, 
                                                     nu = nu_param)
        } else {
          acq_values <- self$acquisition_fn$evaluate(pred$mean, pred$sd, f_best)
        }
        
        # Apply blocking strategy
        if (!is.null(self$blocking_strategy)) {
          if (i > 1) {
            self$blocking_strategy$set_pending(pending)
          }
          
          if (inherits(self$blocking_strategy, "ConstantLiar")) {
            # Constant Liar: augment training data
            if (i > 1) {
              self$blocking_strategy$set_pending(pending, self$model$y, self$model)
              augmented <- self$blocking_strategy$get_augmented_data(
                self$model$X, self$model$y
              )
              
              # Temporarily update model (would need re-fitting)
              # For now, just apply to acquisition
            }
          } else if (inherits(self$blocking_strategy, "TargetedBlocking")) {
            # Targeted blocking: requires setting target point and all evaluated points
            # Set first point in batch as target for subsequent selections
            if (i == 1) {
              # For first point, use best current point as target
              best_current <- self$model$get_best()$best_x
              self$blocking_strategy$set_target(best_current)
            } else {
              # Use first batch point as target
              self$blocking_strategy$set_target(batch[1, ])
            }
            
            acq_values <- self$blocking_strategy$apply_blocking(
              acq_values, candidates, self$model$X, self$model$bounds
            )
          } else if (inherits(self$blocking_strategy, "LocalPenalization")) {
            acq_values <- self$blocking_strategy$apply_blocking(
              acq_values, candidates, pred$mean, pred$sd, self$model$bounds
            )
          } else {
            # Standard blocking with targeted behavior
            self$blocking_strategy$update_recent(self$model$X)
            
            # Determine current best point for targeted blocking
            current_best <- if (i == 1) {
              self$model$get_best()$best_x
            } else {
              batch[1, ]
            }
            
            acq_values <- self$blocking_strategy$apply_blocking(
              acq_values, candidates, self$model$bounds, current_best
            )
          }
        }
        
        # Select best point
        best_idx <- which.max(acq_values)
        
        if (length(best_idx) == 0 || is.infinite(acq_values[best_idx])) {
          # Fallback to random if all blocked
          best_idx <- sample(1:self$n_candidates, 1)
        }
        
        batch[i, ] <- candidates[best_idx, ]
        pending <- rbind(pending, batch[i, , drop = FALSE])
      }
      
      batch
    }
  )
)


#' Parallel Batch Acquisition
#'
#' @description
#' Parallel batch acquisition with support for multiple parallelization backends.
#'
#' @export
ParallelBatchAcquisition <- R6::R6Class(
  "ParallelBatchAcquisition",
  inherit = BatchAcquisition,
  
  public = list(
    parallel_backend = NULL,
    n_cores = NULL,
    cluster = NULL,
    
    #' @description
    #' Initialize parallel batch acquisition
    #' @param parallel_backend Backend: "future", "parallel", "foreach" (default: "future")
    #' @param n_cores Number of cores (default: all available - 1)
    #' @inheritParams BatchAcquisition
    initialize = function(model, acquisition_fn, blocking_strategy = NULL,
                         batch_size = 4, n_candidates = 1000,
                         parallel_backend = "future", n_cores = NULL) {
      super$initialize(model, acquisition_fn, blocking_strategy, 
                      batch_size, n_candidates)
      
      self$parallel_backend <- parallel_backend
      self$n_cores <- if (is.null(n_cores)) {
        parallel::detectCores() - 1
      } else {
        n_cores
      }
      
      self$setup_parallel()
    },
    
    #' @description
    #' Setup parallel backend
    setup_parallel = function() {
      if (self$parallel_backend == "future") {
        future::plan(future::multisession, workers = self$n_cores)
      } else if (self$parallel_backend == "parallel") {
        self$cluster <- parallel::makeCluster(self$n_cores)
      } else if (self$parallel_backend == "foreach") {
        cl <- parallel::makeCluster(self$n_cores)
        doParallel::registerDoParallel(cl)
        self$cluster <- cl
      }
    },
    
    #' @description
    #' Cleanup parallel resources
    cleanup = function() {
      if (!is.null(self$cluster)) {
        parallel::stopCluster(self$cluster)
        self$cluster <- NULL
      }
      
      if (self$parallel_backend == "future") {
        future::plan(future::sequential)
      }
    },
    
    #' @description
    #' Evaluate objective function in parallel
    #' @param batch Matrix of points to evaluate
    #' @param fn Objective function
    #' @return Vector of function values
    evaluate_batch_parallel = function(batch, fn) {
      n_batch <- nrow(batch)
      
      if (self$parallel_backend == "future") {
        # Future backend
        results <- future.apply::future_lapply(1:n_batch, function(i) {
          fn(batch[i, ])
        }, future.seed = TRUE)
        
      } else if (self$parallel_backend == "parallel") {
        # Parallel backend
        parallel::clusterExport(self$cluster, 
                               varlist = c("fn"), 
                               envir = environment())
        results <- parallel::parLapply(self$cluster, 1:n_batch, function(i) {
          fn(batch[i, ])
        })
        
      } else if (self$parallel_backend == "foreach") {
        # Foreach backend
        results <- foreach::foreach(i = 1:n_batch) %dopar% {
          fn(batch[i, ])
        }
        
      } else {
        # Fallback to sequential
        results <- lapply(1:n_batch, function(i) fn(batch[i, ]))
      }
      
      unlist(results)
    }
  )
)
