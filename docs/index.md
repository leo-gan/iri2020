# IRI2020 Rust Port

Welcome to the IRI2020 (International Reference Ionosphere) Python and Rust package documentation!

The project is a 100% pure Rust port of the IRI2020 model, removing the Fortran compilation and compiler dependency entirely.
But the Fortran code is still available for reference. It is fully functional and can be used to verify the Rust implementation.

## Architecture

The project consists of the following main directories:

- `src/data/`: Centralized `.dat` and `.asc` files required to operate the thermosphere/ionosphere mathematical models. The Rust core loads the required files directly from here.
- `src/iri2020/`: The Python module and tests wrapper.
- `src/rust/`: A Rust crate utilizing PyO3 to wrap the IRI2020 mathematical models written in pure Rust and expose them to Python.
- `src/fortran/`: The original Fortran module.
- `docs/`: Documentation for the project.
