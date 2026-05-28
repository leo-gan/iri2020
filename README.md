# IRI2020 ionosphere model ported from Fortran to Rust with Python UI

[![GitHub Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue?style=flat-square)](https://leo-gan.github.io/iri2020/)
[![View iri2020 on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/81056-iri2020)

Python interface to the International Reference Ionosphere (IRI) 2020 model.
The original Fortran code was ported to Rust for better performance and reliability.
The Fortran code is still a part of the repository for reference. 
But the whole Fortran codebase was ported to Rust and the Rust code is now used in the Python interface.

## Install

Prerequisites

* Rust compiler--any modern Rust compiler will do. Here's how to get Rust:
  * Linux: `apt install rustc`
  * Mac: `brew install rust`
  *  Windows Subsystem for Linux

and then install latest release:

```bash
git clone https://github.com/leo-gan/iri2020
cd iri2020
pip install .
```

## Data files

`src/data/`
[data files](https://irimodel.org/indices/IRI-Format-indices-files.pdf)
are
[regularly updated](http://irimodel.org/indices/).
Currently we don't auto-update those.

## Rust Port

This project now includes a 100% pure Rust port of the IRI2020 model, removing the Fortran compilation and compiler dependency entirely. The Rust implementation provides:

- **Thread Safety**: Eliminates Fortran process-global common blocks and mutable globals, allowing safe parallel execution
- **Fast Build Times**: No Fortran compilers (like `gfortran`) or `cmake` wrappers needed
- **Numerical Precision**: Maintains equivalence with the original Fortran implementation through rigorous testing

### Documentation

Detailed documentation is available in the `docs/` directory:

- **[Algorithmic Reference](docs/algorithms.md)**: Scientific foundations and mathematical models implemented in the Rust codebase, including IGRF-13, CIRA-86, ROCSAT-1 drift, and electron density profile construction
- **[Porting Documentation](docs/porting.md)**: Software engineering process for the Fortran to Rust translation, including bug fixes and verification methodology
- **[Package Documentation](docs/package.md)**: Architecture overview, testing procedures, and key features of the Rust implementation

The Rust port is built with `Maturin` (`PyO3`) and can be installed with:

```bash
pip install .
```
