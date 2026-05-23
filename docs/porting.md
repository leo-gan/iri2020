# 2. Porting Documentation

This section provides the most rigorous and detailed descriptions of the decisions made during the architecture and the porting of Fortran to Rust via FFI.

## Options and Trade-Offs

When porting a 3MB Fortran 77 project, we have a few options:
1. **Full Line-by-Line Translation**: Takes a significant amount of time and is highly prone to subtle floating-point truncation regressions.
2. **Automated Transpilation (f2c -> c2rust)**: Fast, but results in thousands of lines of unmaintainable `unsafe` code with raw pointers.
3. **Incremental FFI (Chosen Approach)**: We implemented an architectural shell where Rust drives the Fortran code natively.

### Architectural Decisions

We used `iso_c_binding` in Fortran to expose subroutine signatures:

```fortran
subroutine c_iri_sub(...) bind(C, name="iri_sub_c")
```

This prevents name-mangling problems across different OS/CPU targets. Rust uses a standard `extern "C"` block with `cc/cmake` build scripts to statically link against `libiri.a`.

### FFI Details and Edge Cases

- **Multidimensional Arrays**: Fortran passes multidimensional arrays as 1D column-major blocks. `outf_c(20, 1000)` had to be dynamically flattened and expanded in Rust to maintain accurate column-indexing natively (e.g., `outf[col * 20 + row]`).
- **File Parsing**: We unified `src/data` out of the python specific packaging logic directly to the source level, using exact OS path mapping to fix deep `open(12, ...)` statements nested inside `irifun.for`.

By setting up this FFI, developers can now incrementally port individual subroutines (e.g., `igrf.for` or `cira.for`) into safe Rust.
