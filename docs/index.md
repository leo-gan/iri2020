# IRI2020 Rust Port

Welcome to the IRI2020 (International Reference Ionosphere) Python, Fortran, and Rust wrapper documentation!

The primary goal of this refactoring was to reproduce all existing functionality "as is" and provide an avenue for Rust-Fortran side-by-side comparison, leading to a potential full port.

## Installation

This package is now built with **Maturin** (PyO3). It bundles the Fortran `IRI2020` codebase statically linked into a Rust shared extension via C FFI, providing a massive speed-up compared to the previous subprocess CLI execution.

```bash
pip install .
```

To run both implementations or strictly Rust, refer to the [Porting Documentation](porting.md) and [Package Documentation](package.md).
