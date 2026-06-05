# Package, Math, and Model (Fortran Part)

This part of the documentation covers the legacy Fortran infrastructure, compilation, driver execution, and how the original reference model is run.

## Architecture

The Fortran codebase is located in:

- `src/fortran/`: Contains the original Fortran 77/90 source files, a CMake build configuration, and command-line driver/test scripts.

Key source files include:

- `irisub.for`: Core subroutine (`IRI_SUB`) that orchestrates all sub-model calculations.
- `irifun.for`: Auxiliary mathematical functions, interpolation routines, and profile calculation helpers.
- `iritec.for`: Subroutine (`IRITEC`) for numerical integration of Total Electron Content.
- `igrf.for` & `cira.for` & `rocdrift.for`: External physical and empirical models (IGRF magnetic field, CIRA/MSIS neutral atmosphere, ROCSAT vertical plasma drift).
- `iridreg.for`: Danilov D-region model and lookup table logic.
- `iri_driver.f90` & `test.f90`: Standalone command-line program drivers to execute the simulation and print outputs.
- `CMakeLists.txt`: CMake build configuration to compile the Fortran source code into a static library.

## Compiling and Running

To compile the Fortran library and drivers, you need a Fortran compiler (e.g., `gfortran`) and `CMake`.

### Using CMake

From the repository root:

```bash
mkdir -p src/fortran/build
cd src/fortran/build
cmake ..
make
```

### Direct Compilation with gfortran

Alternatively, you can compile the standalone basic test program directly using `gfortran`:

```bash
cd src/fortran
gfortran -o basictest test.f90 irisub.for irifun.for iritec.for iridreg.for iriflip.for igrf.for cira.for rocdrift.for
./basictest
```

Similarly, the command-line driver `iri_driver` can be compiled and executed:

```bash
gfortran -o iri_driver iri_driver.f90 irisub.for irifun.for iritec.for iridreg.for iriflip.for igrf.for cira.for rocdrift.for
./iri_driver 1980 03 21 12 00 00 0.0 0.0 100.0 500.0 20.0
```

### (Arguments are: year, month, day, hour, minute, second, latitude, longitude, min_alt_km, max_alt_km, step_alt_km)

## Running Tests

The legacy basic test is defined in `test.f90`. When compiled and executed as shown above, it prints key parameters like $N_e$ profile, $N_mF_2$, and $h_mF_2$ values, verifying that the outputs conform to expected ranges and checking for any output length discrepancies.

## Legacy Limitations

- **No Thread Safety**: Extensive use of global mutable variables, `/COMMON/` blocks, and `SAVE` statements prevents concurrent execution across threads.
- **Compiler Dependency**: Requires a Fortran compiler (`gfortran` or similar) and build systems like `CMake` to be installed on the host system.
- **Undefined/Implicit State**: Some variables are implicitly saved or left uninitialized, introducing potential state leakage bugs across sequential invocations.
