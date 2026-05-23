use std::env;
use std::path::PathBuf;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let fortran_dir = PathBuf::from(manifest_dir).join("../fortran");

    let dst = cmake::Config::new(fortran_dir).build_target("iri").build();

    println!("cargo:rustc-link-search=native={}/build", dst.display());
    println!("cargo:rustc-link-lib=static=iri");
    // Link gfortran
    println!("cargo:rustc-link-lib=dylib=gfortran");

    println!("cargo:rerun-if-changed=../fortran");
}
