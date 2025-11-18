# Contributing to airsspy

Thank you for your interest in contributing to airsspy! This document provides guidelines for contributors.

## Development Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/zhubonan/airsspy.git
   cd airsspy
   ```

2. **Set up the development environment:**
   ```bash
   # Using uv (recommended)
   uv venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   uv pip install -e .[dev,test]

   # Or using pip
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   pip install -e .[dev,test]
   ```

3. **Install pre-commit hooks:**
   ```bash
   pre-commit install
   ```

## Code Style and Quality

We use modern Python tooling to maintain code quality:

- **Formatting**: `ruff format` (auto-formats code)
- **Linting**: `ruff check` (finds and fixes issues)
- **Type checking**: `mypy` (static type analysis)

### Running Checks

```bash
# Format code
ruff format airsspy/

# Check for issues
ruff check airsspy/

# Type checking
mypy airsspy/ --ignore-missing-imports

# Run all checks at once
pre-commit run --all-files
```

## Testing

We use `pytest` for testing with coverage reporting:

```bash
# Run all tests
pytest airsspy/

# Run with coverage
pytest airsspy/ --cov=airsspy --cov-report=html

# Run specific test file
pytest airsspy/tests/test_seed.py
```

## Development Workflow

1. Create a new branch for your feature/bugfix
2. Make your changes
3. Run tests and code quality checks
4. Commit your changes (pre-commit hooks will run automatically)
5. Push to your fork and create a pull request

## Code Structure

- `airsspy/seed.py`: Core SeedAtoms class
- `airsspy/build.py`: Buildcell wrapper
- `airsspy/restools.py`: RES file utilities
- `airsspy/common.py`: Shared utilities

When adding new features:
- Follow existing code patterns
- Add type hints where appropriate
- Include tests for new functionality
- Update documentation if needed

## Submitting Changes

- Ensure all tests pass
- Follow the existing code style
- Write clear commit messages
- Include tests for new functionality
- Update documentation as needed

## Questions

Feel free to open an issue for questions or discussion before making major changes.
