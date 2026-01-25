# Contributing to CasPINS

Thank you for your interest in contributing to CasPINS! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

1. Check existing [Issues](https://github.com/InnovationLine/CasPINS/issues) to avoid duplicates
2. Create a new issue with:
   - Clear, descriptive title
   - Steps to reproduce the problem
   - Expected vs. actual behavior
   - Python version and operating system
   - Relevant error messages or screenshots

### Suggesting Enhancements

1. Open an issue describing your proposed enhancement
2. Explain the use case and why it would benefit users
3. If possible, outline a potential implementation approach

### Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Make your changes following our coding standards
4. Add tests for new functionality
5. Ensure all tests pass: `pytest tests/ -v`
6. Update documentation if needed
7. Submit a pull request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/CasPINS.git
cd CasPINS

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
pip install pytest pytest-cov flake8

# Run tests
pytest tests/ -v
```

## Coding Standards

- Follow PEP 8 style guidelines
- Maximum line length: 120 characters
- Use descriptive variable and function names
- Add docstrings to functions and classes
- Write unit tests for new functionality

## Testing

All contributions must include appropriate tests:

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=src

# Run specific test file
pytest tests/test_grna_design.py -v
```

## Code Review Process

1. All pull requests require review before merging
2. Reviewers will check for:
   - Code quality and style
   - Test coverage
   - Documentation updates
   - Breaking changes

## Questions?

If you have questions about contributing, please open an issue or contact the maintainers.

## License

By contributing to CasPINS, you agree that your contributions will be licensed under the MIT License.
