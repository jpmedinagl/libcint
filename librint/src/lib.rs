#![feature(autodiff)]
#![allow(unused_variables,improper_ctypes_definitions,static_mut_refs,non_snake_case)]

use std::path::PathBuf;

// pub mod gen;
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
pub mod linalg;

pub mod cint;

// pub mod c2f;
pub mod cint1e;
pub mod cint2e;
pub mod cint_bas;
pub mod cart2sph;
pub mod optimizer;
pub mod fblas;
pub mod g1e;
pub mod g2e;
pub mod rys_roots;
pub mod rys_wheeler;
pub mod find_roots;
pub mod fmt;
pub mod eigh;
// pub mod misc;
// pub mod g3c1e;
// pub mod g1e_grids;
// pub mod g2c2e;
// pub mod g3c2e;
pub mod intor1;
