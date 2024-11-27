#![allow(non_snake_case, non_upper_case_globals)]

use std::io;

use librint::utils::read_basis;

use librint::cint_bas::CINTcgto_cart;
use librint::cint1e::cint1e_nuc_cart;
use librint::cint1e::cint1e_ovlp_cart;
use librint::intor1::cint1e_kin_cart;
use librint::cint2e::cint2e_cart;

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

fn main() -> io::Result<()> {
    let mut atm: Vec<i32> = Vec::new();
    let mut bas: Vec<i32> = Vec::new();
    let mut env: Vec<f64> = Vec::new();

    let path = "/u/jpmedina/libcint/librint/molecules/h2/sto3g.txt";
    read_basis(path, &mut atm, &mut bas, &mut env)?;

    const natm: usize = 2;
    const nbas: usize = 2;

    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut buf;

    println!("ovlp");
    for i in 0..nbas {
        for j in 0..nbas {
            shls[0] = i as i32;
            shls[1] = j as i32;

            let di = CINTcgto_cart(i, &bas);
            let dj = CINTcgto_cart(j, &bas);

            buf = vec![0.0; (di * dj) as usize];
            cint1e_ovlp_cart(
                &mut buf,
                shls,
                &mut atm,
                natm as i32,
                &mut bas,
                nbas as i32,
                &mut env,
            );

            for i in 0..((di * dj) as usize) {
                print!("{} ", buf[i]);
            }
        }
        println!();
    }

    println!("kin");
    for i in 0..nbas {
        for j in 0..nbas {
            shls[0] = i as i32;
            shls[1] = j as i32;

            let di = CINTcgto_cart(i, &bas);
            let dj = CINTcgto_cart(j, &bas);

            buf = vec![0.0; (di * dj) as usize];
            cint1e_kin_cart(
                &mut buf,
                shls,
                &mut atm,
                natm as i32,
                &mut bas,
                nbas as i32,
                &mut env,
            );

            for i in 0..((di * dj) as usize) {
                print!("{} ", buf[i]);
            }
        }
        println!();
    }

    println!("nuc");
    for i in 0..nbas {
        for j in 0..nbas {
            shls[0] = i as i32;
            shls[1] = j as i32;

            let di = CINTcgto_cart(i, &bas);
            let dj = CINTcgto_cart(j, &bas);

            buf = vec![0.0; (di * dj) as usize];
            cint1e_nuc_cart(
                &mut buf,
                shls,
                &mut atm,
                natm as i32,
                &mut bas,
                nbas as i32,
                &mut env,
            );

            for i in 0..((di * dj) as usize) {
                print!("{} ", buf[i]);
            }
        }
        println!();
    }

    println!("rep");
    for i in 0..nbas {
        for j in 0..nbas {
            for k in 0..nbas {
                for l in 0..nbas {
                    shls[0] = i as i32;
                    shls[1] = j as i32;
                    shls[2] = k as i32;
                    shls[3] = l as i32;

                    let di = CINTcgto_cart(i, &bas);
                    let dj = CINTcgto_cart(j, &bas);
                    let dk = CINTcgto_cart(k, &bas);
                    let dl = CINTcgto_cart(l, &bas);

                    buf = vec![0.0; (di * dj * dk * dl) as usize];

                    cint2e_cart(
                        &mut buf,
                        shls,
                        &mut atm,
                        natm as i32,
                        &mut bas,
                        nbas as i32,
                        &mut env,
                    );

                    for m in 0..((di * dj * dk * dl) as usize) {
                        print!("{} ", buf[m]);
                    }
                }
                println!();
            }
            println!();
        }
        println!();
    }

    Ok(())
}
