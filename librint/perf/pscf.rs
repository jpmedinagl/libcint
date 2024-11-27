#![allow(non_snake_case, non_upper_case_globals)]

use std::io;

use librint::utils::{print_arr, read_basis, split};

use librint::scf::{nmol, angl, density, energy};

fn main() -> io::Result<()> {
    const nelec: usize = 2;

    let mut atm = Vec::new();
    let mut bas = Vec::new();
    let mut env = Vec::new();

    let path = "/u/jpmedina/libcint/librint/molecules/h2/sto3g.txt";
    read_basis(path, &mut atm, &mut bas, &mut env)?;

    let (natm, nbas) = nmol(&atm, &bas);
    let nshells = angl(&bas, 0);
    println!("{} {} {}", natm, nbas, nshells);

    const imax: i32 = 20;
    const conv: f64 = 0.000001;

    let mut P = density(&mut atm, &mut bas, &mut env, nelec, imax, conv);
    print_arr(nshells, 2, &mut P);

    let Etot: f64 = energy(&mut atm, &mut bas, &mut env, &mut P);
    println!("E: {}", Etot);

    Ok(())
}
