cmake_minimum_required(VERSION 3.14...3.27)

project(iri2020 LANGUAGES Fortran C)

add_library(iri STATIC
    irisub.for
    irifun.for
    iritec.for
    iridreg.for
    iriflip.for
    igrf.for
    cira.for
    rocdrift.for
    iri_c_bindings.f90
)

# Position independent code is required when linking static libraries into shared Rust libs (.so)
set_property(TARGET iri PROPERTY POSITION_INDEPENDENT_CODE ON)
