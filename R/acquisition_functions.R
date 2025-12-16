#' Expected Improvement Acquisition Function (Gaussian)
#'
#' @description
#' Classic Expected Improvement (EI) for Bayesian Optimization with Gaussian distribution.
#' Balances exploitation of good regions with exploration of uncertain regions.
#'
#' @export
ExpectedImprovement <- R6::R6Class(
  "ExpectedImprovement",
  
  public = list(
    xi = NULL,
    
    #' @description
    #' Initialize EI
    #' @param xi Exploration parameter (default: 0.01)
    initialize = function(xi = 0.01) {
      self$xi <- xi
    },
    
    #' @description
    #' Evaluate EI at candidate points
    #' @param mean Vector of predictive means
    #' @param std Vector of predictive standard deviations
    #' @param f_best Current best observed value
    #' @param minimize Logical, TRUE for minimization (default: FALSE)
    #' @param nu Degrees of freedom (ignored for Gaussian EI, kept for compatibility)
    #' @return Vector of EI values
    evaluate = function(mean, std, f_best, minimize = FALSE, nu = NULL) {
      if (minimize) {
        mean <- -mean
        f_best <- -f_best
      }
      
      # Avoid division by zero
      std <- pmax(std, 1e-10)
      
      # Calculate improvement
      improvement <- mean - f_best - self$xi
      Z <- improvement / std
      
      # Gaussian EI formula
      ei <- improvement * pnorm(Z) + std * dnorm(Z)
      ei[std == 0] <- 0  # No improvement where no uncertainty
      
      ei
    }
  )
)


#' Student-t Expected Improvement Acquisition Function
#'
#' @description
#' Expected Improvement using Student-t distribution instead of Gaussian.
#' More robust and allows for heavy-tailed exploration controlled by degrees of freedom (nu).
#'
#' @details
#' Student-t EI uses the Student-t CDF and PDF instead of normal distribution:
#' - Low nu (e.g., 3): Heavy tails, aggressive exploration
#' - High nu (e.g., 50): Approaches Gaussian EI
#' - Adaptive nu: Transitions from exploration to exploitation
#'
#' Formula:
#' EI_t(x) = (μ(x) - f* - ξ) * T_ν(Z) + σ(x) * t_ν(Z)
#' where T_ν is Student-t CDF, t_ν is Student-t PDF, Z = (μ(x) - f* - ξ)/σ(x)
#'
#' @export
StudentTExpectedImprovement <- R6::R6Class(
  "StudentTExpectedImprovement",
  
  public = list(
    xi = NULL,
    nu = NULL,
    
    #' @description
    #' Initialize Student-t EI
    #' @param xi Exploration parameter (default: 0.01)
    #' @param nu Degrees of freedom (default: 3 for heavy tails)
    initialize = function(xi = 0.01, nu = 3) {
      self$xi <- xi
      self$nu <- nu
    },
    
    #' @description
    #' Evaluate Student-t EI at candidate points
    #' @param mean Vector of predictive means
    #' @param std Vector of predictive standard deviations
    #' @param f_best Current best observed value
    #' @param minimize Logical, TRUE for minimization (default: FALSE)
    #' @param nu Degrees of freedom (overrides default if provided)
    #' @return Vector of Student-t EI values
    evaluate = function(mean, std, f_best, minimize = FALSE, nu = NULL) {
      if (minimize) {
        mean <- -mean
        f_best <- -f_best
      }
      
      # Use provided nu or default
      nu_use <- if (!is.null(nu)) nu else self$nu
      
      # Avoid division by zero
      std <- pmax(std, 1e-10)
      
      # Calculate improvement
      improvement <- mean - f_best - self$xi
      Z <- improvement / std
      
      # Student-t EI formula
      # T_nu(Z) = Student-t CDF
      # t_nu(Z) = Student-t PDF
      t_cdf <- pt(Z, df = nu_use)
      t_pdf <- dt(Z, df = nu_use)
      
      ei_t <- improvement * t_cdf + std * t_pdf
      ei_t[std == 0] <- 0  # No improvement where no uncertainty
      
      ei_t
    },
    
    #' @description
    #' Update nu parameter (for adaptive schedules)
    #' @param nu New degrees of freedom
    set_nu = function(nu) {
      self$nu <- nu
    }
  )
)


#' Upper Confidence Bound Acquisition Function
#'
#' @description
#' GP-UCB acquisition function. Optimistic bound on the objective.
#'
#' @details
#' UCB(x) = μ(x) + κ·σ(x)
#' - κ controls exploration-exploitation trade-off
#' - Larger κ → more exploration
#' - Theoretical guarantees for cumulative regret
#'
#' @export
UpperConfidenceBound <- R6::R6Class(
  "UpperConfidenceBound",
  
  public = list(
    kappa = NULL,
    adaptive = NULL,
    iteration = 0,
    
    #' @description
    #' Initialize UCB
    #' @param kappa Exploration parameter (default: 2.0)
    #' @param adaptive Use adaptive kappa based on iteration (default: FALSE)
    initialize = function(kappa = 2.0, adaptive = FALSE) {
      self$kappa <- kappa
      self$adaptive <- adaptive
    },
    
    #' @description
    #' Get current kappa value
    #' @param n_dim Number of dimensions
    #' @param n_iter Current iteration
    #' @return Kappa value
    get_kappa = function(n_dim = 1, n_iter = NULL) {
      if (!is.null(n_iter)) {
        self$iteration <- n_iter
      }
      
      if (self$adaptive && self$iteration > 0) {
        # GP-UCB theoretical bound: κ_t = sqrt(2 log(t^(d/2 + 2) π^2 / 3δ))
        delta <- 0.1
        kappa_t <- sqrt(2 * log(self$iteration^(n_dim/2 + 2) * pi^2 / (3 * delta)))
        return(kappa_t)
      }
      
      self$kappa
    },
    
    #' @description
    #' Evaluate UCB at candidate points
    #' @param mean Vector of predictive means
    #' @param std Vector of predictive standard deviations
    #' @param n_dim Number of dimensions (for adaptive kappa)
    #' @param minimize Logical, TRUE for minimization (default: FALSE)
    #' @return Vector of UCB values
    evaluate = function(mean, std, n_dim = 1, minimize = FALSE) {
      kappa <- self$get_kappa(n_dim)
      
      if (minimize) {
        # LCB for minimization
        return(mean - kappa * std)
      }
      
      # UCB for maximization
      mean + kappa * std
    }
  )
)


#' Probability of Improvement Acquisition Function
#'
#' @description
#' Probability that a point improves upon current best.
#'
#' @export
ProbabilityOfImprovement <- R6::R6Class(
  "ProbabilityOfImprovement",
  
  public = list(
    xi = NULL,
    
    #' @description
    #' Initialize PI
    #' @param xi Exploration parameter (default: 0.01)
    initialize = function(xi = 0.01) {
      self$xi <- xi
    },
    
    #' @description
    #' Evaluate PI at candidate points
    #' @param mean Vector of predictive means
    #' @param std Vector of predictive standard deviations
    #' @param f_best Current best observed value
    #' @param minimize Logical, TRUE for minimization (default: FALSE)
    #' @return Vector of PI values
    evaluate = function(mean, std, f_best, minimize = FALSE) {
      if (minimize) {
        mean <- -mean
        f_best <- -f_best
      }
      
      # Avoid division by zero
      std <- pmax(std, 1e-10)
      
      # Calculate Z-score
      Z <- (mean - f_best - self$xi) / std
      
      # Probability of improvement
      pi <- pnorm(Z)
      pi[std == 0] <- 0
      
      pi
    }
  )
)


#' Thompson Sampling Acquisition
#'
#' @description
#' Sample-based acquisition using Thompson Sampling.
#' Draws samples from predictive distribution.
#'
#' @export
ThompsonSampling <- R6::R6Class(
  "ThompsonSampling",
  
  public = list(
    n_samples = NULL,
    use_student_t = NULL,
    nu = NULL,
    
    #' @description
    #' Initialize Thompson Sampling
    #' @param n_samples Number of samples to draw (default: 1)
    #' @param use_student_t Use Student-t sampling (default: FALSE)
    #' @param nu Degrees of freedom for Student-t (default: 3)
    initialize = function(n_samples = 1, use_student_t = FALSE, nu = 3) {
      self$n_samples <- n_samples
      self$use_student_t <- use_student_t
      self$nu <- nu
    },
    
    #' @description
    #' Evaluate Thompson Sampling acquisition
    #' @param mean Vector of predictive means
    #' @param std Vector of predictive standard deviations
    #' @param minimize Logical, TRUE for minimization (default: FALSE)
    #' @return Vector of sampled values (averaged over samples)
    evaluate = function(mean, std, minimize = FALSE) {
      n <- length(mean)
      
      if (self$use_student_t) {
        # Student-t sampling
        samples <- matrix(nrow = self$n_samples, ncol = n)
        
        for (i in 1:n) {
          z <- rnorm(self$n_samples)
          chi_sq <- rchisq(self$n_samples, df = self$nu)
          samples[, i] <- mean[i] + std[i] * z * sqrt(self$nu / chi_sq)
        }
      } else {
        # Gaussian sampling
        samples <- matrix(nrow = self$n_samples, ncol = n)
        for (i in 1:n) {
          samples[, i] <- rnorm(self$n_samples, mean = mean[i], sd = std[i])
        }
      }
      
      # Average samples
      acquisition <- colMeans(samples)
      
      if (minimize) {
        acquisition <- -acquisition
      }
      
      acquisition
    }
  )
)
