#![allow(non_snake_case, non_upper_case_globals)]
// #![feature(autodiff)]

use std::hint::black_box;

use std::io;
use std::time::Instant;

use librint::utils::read_basis;

use librint::cint1e::cint1e_ovlp_cart;
use librint::cint1e::cint1e_nuc_cart;
use librint::intor1::cint1e_kin_cart;
use librint::cint2e::cint2e_cart;
use librint::cint_bas::CINTcgto_cart;

use librint::scf::{angl, integral1e, integral2e, nmol};
use librint::utils::print_arr;

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

// #[no_mangle]
// #[autodiff(dovlp, Reverse, Duplicated, Const, Const, Const, Const, Duplicated)]
// pub fn ovlp(
//     out: &mut Vec<f64>,
//     shls: [i32; 4],
//     atm: &mut Vec<i32>,
//     bas: &mut Vec<i32>,
//     env1: &mut Vec<f64>,
//     env2: &mut Vec<f64>,
// ) {
//     let (natm, nbas) = nmol(atm, bas);
//     let mut env: Vec<f64> = combine(&env1, &env2);
//     cint1e_ovlp_cart(out, shls, atm, natm as i32, bas, nbas as i32, &mut env);
// }

fn main() -> io::Result<()> {
    let mut atm: Vec<i32> = Vec::new();
    let mut bas: Vec<i32> = Vec::new();
    let mut env: Vec<f64> = Vec::new();

    let path = "/u/jpmedina/libcint/librint/molecules/h2o/sto3g.txt";
    read_basis(path, &mut atm, &mut bas, &mut env)?;

    println!("Single");

    let (natm, nbas) = nmol(&atm, &bas);

    let i = 0;
    let j = 0;
    let mut shls: [i32; 4] = [0, 0, 0, 0];
    shls[0] = i as i32;
    shls[1] = j as i32;
    let di = black_box(CINTcgto_cart(i, &bas) as usize);
    let dj = black_box(CINTcgto_cart(j, &bas) as usize);

    let mut buf = vec![0.0; di * dj];
    black_box(cint1e_ovlp_cart(&mut buf, shls, &mut atm, natm as i32, &mut bas, nbas as i32, &mut env));

    // let mut buf = vec![0.0; (di * dj) as usize];
    // let mut dbuf = vec![0.0; (di * dj) as usize];
    // let mut denv = vec![0.0; env2.len()];
    // dovlp(&mut buf, &mut dbuf, shls, &mut atm, &mut bas, &mut env1, &mut env2, &mut denv);

    // println!("Full");

    // let nshells_cart = angl(&mut bas, 0);
    // let nshells_sph = angl(&mut bas, 1);

    // println!("ovlp cart");
    // let mut S = integral1e(&mut atm, &mut bas, &mut env, 0, 0);
    // print_arr(nshells_cart, 2, &mut S);

    // println!("repulsion cart");
    // let mut R = integral2e(&mut atm, &mut bas, &mut env, 0);
    // print_arr(nshells_cart, 4, &mut R);

    // println!("repulsion sph");
    // R = integral2e(&mut atm, &mut bas, &mut env, 1);
    // print_arr(nshells_sph, 4, &mut R);

    Ok(())
}
