# 2. Porting Documentation

This section provides the details of the porting process from Fortran 77 to safe, pure Rust.

## Porting Strategy

The codebase was ported incrementally from Fortran to Rust:
1. **Side-by-Side FFI Testing**: Initially, a Fortran compilation and C FFI bridge was established. Subroutines were translated to Rust one by one and checked for numerical equivalence against the original Fortran implementation under strict tolerances.
2. **Full Translation**: After validating all subroutines (including IGRF, CIRA, Rocdrift, CCIR/URSI data parsers, and all physical models like electron temperature and ion composition), we replaced all remaining FFI calls with native Rust implementations.
3. **Fortran Removal**: The Fortran source code, FFI bindings, and build-time Fortran compiler settings were entirely removed.

## Architecture and Numerical Precision
- **Multidimensional Arrays**: Fortran passes multidimensional arrays as flat column-major blocks. The Rust codebase uses standard layouts and flat array indexing to match the original behavior and data layout.
- **File Parsing**: Handled natively in Rust with custom data readers that parse the data files in `src/data/` exactly as the original Fortran code did.
- **Numerical Correctness**: The pure Rust implementation has been verified to match original Fortran calculations within a strict `1e-4` relative tolerance across multiple regression test scenarios.
