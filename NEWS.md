# NEWS

## BayesOptR 0.1.0 (2025-12-16)

### Major Features

* Initial release of BayesOptR
* Three Bayesian Optimization methods:
  - Gaussian BO with Gaussian Process surrogate
  - Student-t BO with Student-t Process surrogate
  - Adaptive Hybrid BO (AHBO) with GP surrogate + Adaptive Student-t EI

### Targeted Blocking Innovation

* **Targeted Blocking**: New default behavior for distance-based and dimension-wise blocking
  - Automatically detects current best point (f⁺) and closest known point (g_min)
  - Penalizes only region between f⁺ and g_min
  - Maintains full exploration elsewhere
  - Prevents redundant sampling between known points
  - Controlled by `use_targeted` parameter (default: TRUE)

* **Blocking Strategies**:
  - Distance-based hard blocking (with targeted mode)
  - Distance-based soft blocking (with targeted mode)
  - Dimension-wise hard blocking (with targeted mode)
  - Dimension-wise soft blocking (with targeted mode)
  - Local penalization (Lipschitz-based)
  - Constant liar heuristics

### Student-t Expected Improvement

* Correct implementation of Student-t Expected Improvement acquisition function
* Formula: `EI_t(x) = (μ(x) - f*) × T_ν(Z) + σ(x) × t_ν(Z) × √[(ν + Z² - 1) / ν]`
* Convergence to Gaussian EI as ν → ∞
* Used by AHBO with adaptive ν parameter

### Acquisition Functions

* Expected Improvement (Gaussian and Student-t variants)
* Upper Confidence Bound (UCB) with adaptive κ
* Probability of Improvement (PI)
* Thompson Sampling (Gaussian and Student-t variants)

### Parallelization

* Batch acquisition support with multiple parallel backends
* Sequential greedy batch acquisition
* Parallel batch acquisition with future, parallel, and foreach backends
* Automatic cluster management

### High-Level Interface

* `bayes_optimize()` function for easy optimization
* Comprehensive parameter control
* Support for minimization and maximization
* Flexible acquisition and blocking configuration

### Documentation

* Complete README with examples
* Mathematical foundations explained
* Usage examples for all major features
* Package structure documented

### Testing

* Comprehensive test suite with 25+ tests
* RMarkdown test report for visual validation
* Integration tests for full optimization workflows
* Edge case testing

### Known Limitations

* Currently optimized for continuous optimization
* No built-in constraint handling (planned for v0.2.0)
* Single-objective only (multi-objective planned)
