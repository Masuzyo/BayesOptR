# BayesOptR: Advanced Bayesian Optimization in R

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)

Advanced Bayesian Optimization framework implementing **Gaussian Process**, **Student-t Process**, and **Adaptive Hybrid Bayesian Optimization** with **targeted blocking strategies** and **parallelization support** for efficient batch acquisition.

## Features

### Bayesian Optimization Methods

- **Gaussian BO**: Classic GP-based optimization with Gaussian likelihood
- **Student-t BO**: Robust optimization with heavy-tailed distributions (configurable ν)
- **Adaptive Hybrid BO (AHBO)**: GP surrogate + Adaptive Student-t EI acquisition

### Targeted Blocking Strategies

All distance-based and dimension-wise blocking strategies now support **targeted blocking** by default:

- **Targeted Blocking**: Penalizes only the region between current best (f⁺) and closest known point (g_min)
  - Prevents redundant sampling between known points
  - Maintains full exploration elsewhere
  - Formula: `penalty = |f - g_min| / |f⁺ - g_min|` for points between f⁺ and g_min

**Available Strategies:**
- **Dimension-Wise Hard Blocking**: Complete exclusion with targeted behavior (`use_targeted = TRUE`)
- **Dimension-Wise Soft Blocking**: Penalty-based with targeted behavior (`use_targeted = TRUE`)
- **Distance-Based Hard Blocking**: Exclusion spheres with targeted behavior (`use_targeted = TRUE`)
- **Distance-Based Soft Blocking**: Smooth penalties with targeted behavior (`use_targeted = TRUE`)
- **Local Penalization**: Lipschitz-based penalization
- **Constant Liar**: Fantasy point heuristics (min/mean/max/kriging)

### Acquisition Functions

- **Expected Improvement (EI)** - Gaussian and Student-t variants
- **Student-t Expected Improvement** - Robust heavy-tailed acquisition for AHBO
- **Upper Confidence Bound (UCB)** - With adaptive κ
- **Probability of Improvement (PI)**
- **Thompson Sampling** - Gaussian and Student-t variants

### Parallelization

- **Multiple backends**: `future`, `parallel`, `foreach`
- **Batch evaluation**: Efficient parallel function evaluation
- **Automatic resource management**: Cluster setup and cleanup

## Mathematical Foundation

### Student-t Expected Improvement

**AHBO Key Innovation**: GP surrogate model + Adaptive Student-t EI acquisition function

**Student-t Expected Improvement:**
```
EI_t(x) = (μ(x) - f* - ξ) × T_ν(Z) + σ(x) × t_ν(Z) × √[(ν + Z² - 1) / ν]
where Z = (μ(x) - f* - ξ) / σ(x)
```

- `T_ν(Z)`: Student-t CDF with ν degrees of freedom
- `t_ν(Z)`: Student-t PDF with ν degrees of freedom
- **ν → ∞**: Converges to Gaussian EI
- **Low ν (3-5)**: Heavy tails for aggressive exploration
- **High ν (30-50)**: Near-Gaussian for exploitation

**AHBO Adaptation**: ν transitions from exploration (low) to exploitation (high) over iterations

### Targeted Blocking

**Concept**: Focus penalization only between f⁺ (selected target) and g_min (closest known point)

**Mathematical Formulation:**
```
Blocking Region = {f : |f⁺ - f| < |f⁺ - g_min| AND |f - g_min| < |f⁺ - g_min|}

penalty(f) = |f - g_min| / |f⁺ - g_min|    if f in blocking region
penalty(f) = 1.0                            otherwise
```

**Advantages:**
- More focused than global blocking
- Prevents redundant sampling between known points
- Maintains full exploration elsewhere
- Automatic target and closest point detection

## Installation

```r
# Install from GitHub
devtools::install_github("Masuzyo/BayesOptR")

# Or install dependencies manually
install.packages(c("DiceKriging", "lhs", "parallel", "future", 
                   "future.apply", "foreach", "doParallel",
                   "MASS", "mvtnorm", "Matrix", "matrixStats", "R6"))
```

## Quick Start

### Simple Sequential Optimization

```r
library(BayesOptR)

# Define objective function (Branin)
branin <- function(x) {
  x1 <- x[1] * 15 - 5
  x2 <- x[2] * 15
  -(1 * (x2 - 5.1/(4*pi^2) * x1^2 + 5/pi * x1 - 6)^2 + 
    10 * (1 - 1/(8*pi)) * cos(x1) + 10)
}

# Define bounds [0, 1]^2
bounds <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)

# Run AHBO optimization
result <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "adaptive_hybrid",
  acquisition = "ei",
  n_iter = 50,
  n_init = 10,
  verbose = TRUE
)

# View results
print(result$best_x)
print(result$best_y)
```

### Batch Optimization with Targeted Blocking

```r
# Parallel batch optimization with targeted distance blocking
result <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "adaptive_hybrid",
  acquisition = "ei",
  n_iter = 25,
  batch_size = 4,
  blocking_strategy = "distance_soft",
  use_targeted = TRUE,  # Default: TRUE
  parallel = TRUE,
  n_cores = 4,
  verbose = TRUE
)
```

### Disable Targeted Blocking

```r
# Use traditional distance blocking without targeted behavior
result <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "gaussian",
  acquisition = "ei",
  n_iter = 25,
  batch_size = 4,
  blocking_strategy = "distance_hard",
  use_targeted = FALSE,  # Disable targeted blocking
  radius = 0.1,
  verbose = TRUE
)
```

### Low-Level API

```r
# Create model directly
model <- AdaptiveHybridBO$new(
  bounds = bounds,
  nu_min = 3,
  nu_max = 50,
  nu_schedule = "linear",
  max_iterations = 50,
  n_init = 10
)

# Create acquisition function
acq_fn <- StudentTExpectedImprovement$new(nu = 3)

# Create blocking strategy with targeted behavior
blocking <- DistanceSoftBlocking$new(
  radius = 0.1,
  penalty_strength = 2.0,
  use_targeted = TRUE  # Automatic target detection
)

# Create batch acquisition
batch_acq <- BatchAcquisition$new(
  model = model,
  acquisition_fn = acq_fn,
  blocking_strategy = blocking,
  batch_size = 4
)

# Optimize
for (i in 1:50) {
  # Get next batch
  X_next <- batch_acq$acquire_batch()
  
  # Evaluate function
  y_next <- apply(X_next, 1, branin)
  
  # Update model
  model$update(X_next, y_next)
  
  cat(sprintf("Iteration %d: Best = %.4f\n", i, max(model$y)))
}
```

## Advanced Examples

### Compare Blocking Strategies

```r
strategies <- c("none", "distance_hard", "distance_soft", 
                "dimwise_hard", "dimwise_soft")

results <- lapply(strategies, function(strat) {
  bayes_optimize(
    fn = branin,
    bounds = bounds,
    method = "adaptive_hybrid",
    n_iter = 30,
    batch_size = 4,
    blocking_strategy = strat,
    verbose = FALSE
  )
})

names(results) <- strategies

# Compare convergence
best_values <- sapply(results, function(r) r$best_y)
print(best_values)
```

### Custom ν Schedule for AHBO

```r
# Exponential decay schedule
result_exp <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "adaptive_hybrid",
  nu_schedule = "exponential",
  n_iter = 50
)

# Sigmoid schedule (smooth transition)
result_sigmoid <- bayes_optimize(
  fn = branin,
  bounds = bounds,
  method = "adaptive_hybrid",
  nu_schedule = "sigmoid",
  n_iter = 50
)
```

## Package Structure

```
BayesOptR/
├── R/
│   ├── acquisition_functions.R      # EI, UCB, PI, Thompson
│   ├── batch_acquisition.R          # Batch acquisition with blocking
│   ├── blocking_strategies.R        # Base blocking classes
│   ├── bo_methods.R                 # GaussianBO, StudentTBO, AHBO
│   ├── dimension_wise_blocking.R    # Dimension-wise strategies
│   ├── distance_blocking.R          # Distance-based strategies
│   ├── main_interface.R             # High-level bayes_optimize()
│   └── utils.R                      # Helper functions
├── tests/
│   └── testthat/
│       └── test_*.R                 # Unit tests
├── man/                             # Documentation
├── DESCRIPTION                      # Package metadata
├── NAMESPACE                        # Exports
├── README.md                        # This file
└── LICENSE                          # MIT License
```

## Testing

```r
# Run all tests
devtools::test()

# Run comprehensive test script
source("test_package.R")

# Render test report
rmarkdown::render("test_package_comprehensive.Rmd")
```

## Performance Tips

1. **Use AHBO for complex landscapes**: Better exploration-exploitation balance
2. **Enable targeted blocking**: More efficient batch acquisition (`use_targeted = TRUE` by default)
3. **Choose batch size wisely**: Balance parallelization overhead vs. speedup
4. **Tune ν parameter**: Lower ν (3-5) for noisy/multimodal, higher ν (30-50) for smooth
5. **Use soft blocking**: More flexible than hard blocking for diverse batches
6. **Adjust blocking radius**: Larger radius for smoother functions, smaller for complex landscapes

## References

- **Adaptive Hybrid BO**: Combining GP with Student-t acquisition
- **Targeted Blocking**: Focused penalization between best and closest points
- **Student-t Processes**: Heavy-tailed GPs for robust optimization
- **González et al. (2016)**: "Batch Bayesian Optimization via Local Penalization"
- **Ginsbourger et al. (2010)**: "Kriging is well-suited to parallelize optimization"

## Citation

```bibtex
@software{bayesoptr2025,
  author = {Masuzyo},
  title = {BayesOptR: Advanced Bayesian Optimization in R},
  year = {2025},
  url = {https://github.com/Masuzyo/BayesOptR}
}
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions welcome! Please open an issue or pull request.

## Contact

For questions or issues, please open a GitHub issue.
