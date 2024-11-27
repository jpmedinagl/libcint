#![allow(non_snake_case, non_upper_case_globals)]
#![feature(autodiff)]

use std::autodiff::autodiff;

use librint::utils::read_basis;

use librint::cint1e::cint1e_ovlp_cart;
use librint::cint_bas::CINTcgto_cart;

use librint::scf::{angl, nmol};

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

#[no_mangle]
#[autodiff(
    dcart, Reverse, Duplicated, Const, Const, Const, Const, Const, Duplicated
)]
fn cart(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) {
    cint1e_ovlp_cart(out, shls, atm, natm, bas, nbas, env);
}

fn main() {
    let mut atm = Vec::new();
    let mut bas = Vec::new();
    let mut env = Vec::new();

    let path = "/u/jpmedina/libcint/librint/molecules/h2o/sto3g.txt";
    read_basis(path, &mut atm, &mut bas, &mut env);

    let (natm, nbas) = nmol(&atm, &bas);
    let nshells = angl(&bas, 0);

    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut buf;
    let mut dbuf;

    let i = 0;
    let j = 0;

    shls[0] = i as i32;
    let di = CINTcgto_cart(i, &bas);
    shls[1] = j as i32;
    let dj = CINTcgto_cart(j, &bas);

    println!("ijdidj {}{}{}{}", i, j, di, dj);

    buf = vec![0.0; (di * dj) as usize];
    dbuf = vec![0.0; (di * dj) as usize];

    println!("denv:");
    for k in 0..((di * dj) as usize) {
        dbuf[k] = 1.0;

        let mut denv: [f64; 56] = [0.0; 56];
        dcart(
            &mut buf,
            &mut dbuf,
            shls,
            &mut atm,
            natm as i32,
            &mut bas,
            nbas as i32,
            &mut env,
            &mut denv,
        );
        for i in 49..50 {
            print!("{:.6} ", denv[i]);
        }
        dbuf[k] = 0.0;
    }
    println!();

    // let h = 0.000001;

    // let mut b1 = vec![0.0; (di * dj) as usize];
    // let mut b2 = vec![0.0; (di * dj) as usize];

    // println!("finite d:");
    // for k in 49..50 {
    //     env[k] += h;
    //     cint1e_ovlp_cart(
    //         &mut buf,
    //         shls,
    //         &mut atm,
    //         natm as i32,
    //         &mut bas,
    //         nbas as i32,
    //         &mut env,
    //     );
    //     for l in 0..(di * dj) as usize {
    //         b1[l] = buf[l];
    //     }

    //     env[k] -= 2.0 * h;
    //     cint1e_ovlp_cart(
    //         &mut buf,
    //         shls,
    //         &mut atm,
    //         natm as i32,
    //         &mut bas,
    //         nbas as i32,
    //         &mut env,
    //     );
    //     for l in 0..(di * dj) as usize {
    //         b2[l] = buf[l];
    //     }

    //     env[k] += h;

    //     for l in 0..(di * dj) as usize {
    //         let grad = (b1[l] - b2[l]) / (2.0 * h);
    //         print!("{:.6} ", grad);
    //     }
    //     println!();
    // }
}
