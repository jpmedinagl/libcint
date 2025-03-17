#![feature(autodiff)]
#![allow(unused_variables,improper_ctypes_definitions,static_mut_refs,non_snake_case)]

use std::path::PathBuf;

// pub mod p2r; // pyo3 bindings

pub fn get_path() -> PathBuf {
    let p = std::env::current_dir().unwrap();
    // now apend to find sto3g.txt
    let path = p.join("molecules/h2o/sto3g.txt");
    assert!(path.exists());
    path
}

pub mod scf;

pub mod dscf;
pub mod p2c;

pub mod utils;

pub mod cint;

pub mod cart2sph;
pub mod cint1e;
pub mod intor1;
pub mod cint2e;
pub mod g1e;
pub mod g2e;
pub mod cint_bas;
pub mod eigh;
pub mod fblas;
pub mod find_roots;
pub mod fmt;
pub mod optimizer;
pub mod rys_roots;
pub mod rys_wheeler;