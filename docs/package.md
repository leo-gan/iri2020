# 1. Package, Math, and Model (Rust Part)

This part of the documentation covers the new Rust infrastructure, how the data is shared between Fortran and Rust, and how you can run tests for the new Rust code.

## Architecture

The project consists of three main directories:
- `src/fortran/`: The raw Fortran 77 source code that implements the complex matrix math, models (e.g., CIRA, IGRF).
- `src/rust/`: A Rust crate utilizing PyO3 and `cmake` to bind to the Fortran libraries statically via C FFI.
- `src/data/`: Centralized `.dat` and `.asc` files required to operate the thermosphere/ionosphere mathematical models. Both Fortran C Bindings and future Rust rewrites source files from here.

## Running Tests

Tests are maintained in `src/iri2020/tests`. You can easily run tests side-by-side on the old mathematical models that use PyO3 integration:

```bash
pytest src/iri2020/tests -v
```

These tests will span edge cases, altitude limits, bounds checking for deep bit-and-byte math matching exactly the original text output produced by `iri2020_driver`.

## Future Improvements

Right now, the architecture executes the identical Fortran code. However, the Rust integration natively enables the following high-potential optional improvements:
- **Parallelization via Rayon**: Instead of running a serial loop across `glat` or `time`, the `altkm` computations can be split into chunks across threads.
- **SIMD**: Rewriting the deep bit-and-byte routines (e.g., `irifun.for` integration methods) directly in Rust with AVX/SIMD instructions.
- **Togglable Strict Math**: The option to enforce Rust's strict IEEE 754 float guarantees versus Fortran's `-ffast-math` which often drops precision in `IRI_SUB`.
