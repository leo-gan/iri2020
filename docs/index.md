# IRI2020 Rust Port

Welcome to the IRI2020 (International Reference Ionosphere) Python and Rust package documentation!

The project is a 100% pure Rust port of the IRI2020 model, removing the Fortran compilation and compiler dependency entirely.

## Installation

This package is built with **Maturin** (PyO3). It compiles the pure Rust `iri2020-rust` core library and binds it into a Python extension.

```bash
pip install .
```

To run the implementation, refer to the [Porting Documentation](porting.md) and [Package Documentation](package.md).
