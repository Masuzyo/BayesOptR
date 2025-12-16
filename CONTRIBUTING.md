# Contributing to BayesOptR

Thank you for considering contributing to BayesOptR! This document provides guidelines for contributing to the project.

## Code of Conduct

Please be respectful and constructive in all interactions.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue with:
- Clear description of the problem
- Minimal reproducible example
- Expected vs actual behavior
- R version and package versions
- Operating system

### Suggesting Enhancements

Enhancement suggestions are welcome! Please:
- Clearly describe the enhancement
- Explain why it would be useful
- Provide examples if possible

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass (`devtools::test()`)
6. Update documentation
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "testthat", "roxygen2", "pkgdown"))

# Clone repository
# git clone https://github.com/Masuzyo/BayesOptR.git
# cd BayesOptR

# Install package dependencies
devtools::install_deps()

# Load package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Coding Standards

### R Code Style

- Follow [tidyverse style guide](https://style.tidyverse.org/)
- Use meaningful variable and function names
- Add comments for complex logic
- Keep functions focused and modular

### Documentation

- Use roxygen2 for function documentation
- Include `@param`, `@return`, `@examples`
- Update README.md for user-facing changes
- Add vignettes for major features

### Testing

- Write tests for new functions using testthat
- Aim for high code coverage
- Test edge cases and error handling
- Include integration tests for complex features

Example test structure:
```r
test_that("function does what it should", {
  # Setup
  input <- create_test_input()
  
  # Execute
  result <- your_function(input)
  
  # Verify
  expect_equal(result, expected_value)
  expect_true(is.valid(result))
})
```

## Package Structure

```
BayesOptR/
├── R/                          # Source code
│   ├── bo_methods.R           # BO method classes
│   ├── acquisition_functions.R # Acquisition functions
│   ├── distance_blocking.R    # Distance-based blocking
│   ├── dimension_wise_blocking.R # Dimension-wise blocking
│   ├── batch_acquisition.R    # Batch acquisition
│   └── main_interface.R       # High-level API
├── tests/                      # Test files
│   └── testthat/
├── man/                        # Documentation (auto-generated)
├── inst/                       # Additional files
│   └── examples/              # Example scripts
├── DESCRIPTION                 # Package metadata
├── NAMESPACE                   # Package exports
└── README.md                   # Main documentation
```

## Areas for Contribution

### High Priority

- Additional blocking strategies
- More acquisition functions (e.g., Knowledge Gradient)
- Constraint handling for constrained optimization
- Multi-objective optimization support
- Better visualization tools

### Documentation

- More examples and use cases
- Vignettes explaining theory and practice
- Benchmarks against other BO packages
- Case studies with real-world applications

### Performance

- Optimize computational bottlenecks
- Better parallelization strategies
- Efficient memory usage for large datasets
- GPU support for expensive computations

### Testing

- Expand test coverage
- Add property-based tests
- Stress tests for edge cases
- Benchmarking suite

## Questions?

Feel free to open an issue for any questions about contributing!

## License

By contributing, you agree that your contributions will be licensed under the MIT License.
