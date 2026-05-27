# 1. Package, Math, and Model (Rust Part)

This part of the documentation covers the Rust infrastructure, model execution, and how to run tests.

## Architecture

The project consists of three main directories:
- `src/rust/`: A Rust crate utilizing PyO3 to wrap the IRI2020 mathematical models written in pure Rust and expose them to Python.
- `src/data/`: Centralized `.dat` and `.asc` files required to operate the thermosphere/ionosphere mathematical models. The Rust core loads the required files directly from here.
- `src/iri2020/`: The Python module and tests wrapper.

## Running Tests

Tests are maintained in `src/iri2020/tests` (Python integration tests) and `src/rust/tests` (Rust integration tests).

To run the Python test suite:
```bash
poetry run pytest
```

To run the Rust integration tests:
```bash
cd src/rust
cargo test
```

These tests verify bounds checking, extreme altitude inputs, and compare calculation outputs against the golden regression fixtures.

## Key Features & Benefits
- **Thread Safety**: The original Fortran process-global common blocks and mutable globals have been replaced with stateless or encapsulated Rust structs, allowing safe parallel execution.
- **Fast Build Times**: No Fortran compilers (like `gfortran`) or `cmake` wrappers are needed anymore. The package builds directly using standard Rust and Python tools.
