# Contributing to REdiCAT

We welcome contributions to REdiCAT! This document provides guidelines for contributing to the project.

## Table of Contents

- [Contributing to REdiCAT](#contributing-to-redicat)
  - [Table of Contents](#table-of-contents)
  - [Getting Started](#getting-started)
  - [How to Contribute](#how-to-contribute)
    - [Reporting Bugs](#reporting-bugs)
    - [Suggesting Enhancements](#suggesting-enhancements)
    - [Code Contributions](#code-contributions)
  - [Development Setup](#development-setup)
    - [Prerequisites](#prerequisites)
    - [Building the Project](#building-the-project)
  - [Coding Standards](#coding-standards)
    - [Rust Style Guide](#rust-style-guide)
    - [Documentation Comments](#documentation-comments)
  - [Testing](#testing)
    - [Writing Tests](#writing-tests)
    - [Running Tests](#running-tests)
  - [Documentation](#documentation)
    - [Code Documentation](#code-documentation)
    - [User Documentation](#user-documentation)
  - [Pull Request Process](#pull-request-process)
    - [Pull Request Guidelines](#pull-request-guidelines)
    - [Code Review Process](#code-review-process)
  - [Community](#community)
    - [Communication](#communication)
    - [Code of Conduct](#code-of-conduct)
  - [License](#license)

## Getting Started

1. Fork the repository
2. Clone your fork: `git clone https://github.com/your-username/redicat.git`
3. Create a new branch: `git checkout -b feature/your-feature-name`
4. Make your changes
5. Commit your changes: `git commit -am 'Add some feature'`
6. Push to the branch: `git push origin feature/your-feature-name`
7. Create a new Pull Request

## How to Contribute

### Reporting Bugs

Before reporting a bug, please check if it has already been reported in the [issue tracker](https://github.com/aStudyingTurtle/perbase_redicat/issues).

When reporting a bug, please include:

- A clear and descriptive title
- A detailed description of the problem
- Steps to reproduce the issue
- Expected behavior
- Actual behavior
- System information (OS, Rust version, etc.)
- Any relevant error messages or logs

### Suggesting Enhancements

Feature requests are welcome! Please provide:

- A clear and descriptive title
- A detailed explanation of the proposed feature
- The problem it solves or the benefit it provides
- Any implementation suggestions (if applicable)

### Code Contributions

We welcome code contributions in the form of bug fixes, new features, and improvements to existing functionality.

## Development Setup

### Prerequisites

- Rust (latest stable version)
- Cargo
- Git
- samtools (for creating FASTA index files)

### Building the Project

```bash
# Clone the repository
git clone https://github.com/aStudyingTurtle/perbase_redicat.git
cd perbase_redicat

# Build the project
cargo build

# Run tests
cargo test
```

## Coding Standards

### Rust Style Guide

We follow the standard Rust coding conventions:

- Use `rustfmt` to format your code: `cargo fmt`
- Follow the [Rust API Guidelines](https://rust-lang.github.io/api-guidelines/)
- Write documentation for public APIs
- Use meaningful variable and function names
- Keep functions focused and small
- Handle errors appropriately

### Documentation Comments

All public functions, structs, and modules should have documentation comments:

```rust
/// Brief description of what this function does
///
/// More detailed explanation if needed.
///
/// # Arguments
///
/// * `arg1` - Description of arg1
/// * `arg2` - Description of arg2
///
/// # Returns
///
/// Description of return value
///
/// # Example
///
/// ```
/// let result = my_function(arg1, arg2);
/// ```
pub fn my_function(arg1: Type1, arg2: Type2) -> ReturnType {
    // Implementation
}
```

## Testing

### Writing Tests

All new functionality should include appropriate tests. Tests should be:

- Comprehensive and cover edge cases
- Fast and reliable
- Clear and readable
- Independent of each other

We use several testing frameworks:

1. **Unit Tests**: Standard Rust unit tests using `#[cfg(test)]`
2. **Integration Tests**: Tests in the `tests/` directory
3. **Property-based Tests**: Using `proptest` for property-based testing
4. **Benchmark Tests**: Using `criterion` for performance benchmarks

Example unit test:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_my_function() {
        let result = my_function(2, 3);
        assert_eq!(result, 5);
    }
}
```

Example property-based test:

```rust
#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_my_function_prop(a in 0..1000i32, b in 0..1000i32) {
            let result = my_function(a, b);
            assert_eq!(result, a + b);
        }
    }
}
```

### Running Tests

```bash
# Run all tests
cargo test

# Run specific tests
cargo test test_name

# Run tests with output
cargo test -- --nocapture

# Run only unit tests (exclude integration tests)
cargo test --lib

# Run only integration tests
cargo test --test '*'

# Run benchmarks
cargo bench
```

## Documentation

### Code Documentation

All public APIs should be documented with Rust documentation comments.

To generate and view documentation locally:

```bash
# Generate documentation
cargo doc --open
```

### User Documentation

User-facing documentation should be updated when:

- Adding new features
- Changing existing functionality
- Fixing bugs that affect user experience

Documentation files are located in the `docs/` directory.

To check documentation quality:

```bash
# Check documentation for lints
cargo doc --dry-run
```

## Pull Request Process

1. Ensure your code follows the coding standards
2. Add tests for new functionality
3. Update documentation as needed
4. Run all tests to ensure nothing is broken: `cargo test`
5. Format your code: `cargo fmt`
6. Create a pull request with a clear title and description
7. Address any feedback from reviewers

### Pull Request Guidelines

- Keep PRs focused on a single feature or bug fix
- Write a clear, descriptive title
- Include a detailed description of the changes
- Reference any related issues
- Ensure all CI checks pass
- Be responsive to feedback during review

### Code Review Process

All pull requests must be reviewed and approved by at least one maintainer before merging. The review process focuses on:

- Code quality and correctness
- Adherence to coding standards
- Test coverage
- Documentation completeness
- Performance considerations
- Security implications

Reviewers will use the following checklist:

- [ ] Code follows Rust style guidelines
- [ ] Public APIs are well documented
- [ ] Tests cover new functionality
- [ ] Existing tests still pass
- [ ] Documentation is updated
- [ ] No unnecessary dependencies added
- [ ] Error handling is appropriate
- [ ] Performance considerations addressed

## Community

### Communication

- For general discussion, use the project's issue tracker
- For security vulnerabilities, contact the maintainers directly

### Code of Conduct

Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms.

## License

By contributing to REdiCAT, you agree that your contributions will be licensed under the MIT License.