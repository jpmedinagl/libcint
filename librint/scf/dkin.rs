#![allow(non_snake_case, non_upper_case_globals)]
#![feature(autodiff)]
use std::autodiff::autodiff;
use std::io;
use std::time::Instant;

use librint::utils::read_basis;

use librint::cint_bas::CINTcgto_cart;
use librint::cint1e::cint1e_ovlp_cart;
use librint::cint1e::cint1e_nuc_cart;
use librint::intor1::cint1e_kin_cart;

use librint::cint_bas::CINTcgto_spheric;
use librint::cint1e::cint1e_ovlp_sph;
use librint::cint1e::cint1e_nuc_sph;
use librint::intor1::cint1e_kin_sph;

use librint::scf::{nmol, angl, density};
use librint::utils::split;

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

#[no_mangle]
#[autodiff(dkinw, Reverse, Duplicated, Const, Const, Const, Const, Duplicated)]
fn kinw(
    out: &mut Vec<f64>, 
    shls: &mut Vec<i32>, 
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>, 
    env1: &mut Vec<f64>,
    env2: &mut Vec<f64>,
) {
    let mut env: Vec<f64> = vec![0.0; env1.len() + env2.len()];
    let mut c = 0;
    for i in 0..env1.len() {
        env[c] = env1[i];
        c += 1;
    }
    for j in 0..env2.len() {
        env[c] = env2[j];
        c += 1;
    }

    let (natm, nbas) = nmol(atm, bas);

    cint1e_kin_cart(out, shls, atm, natm as i32, bas, nbas as i32, &mut env);
}

#[no_mangle]
#[autodiff(dnucw, Reverse, Duplicated, Const, Const, Const, Const, Duplicated)]
fn nucw(
    out: &mut Vec<f64>, 
    shls: &mut Vec<i32>, 
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>, 
    env1: &mut Vec<f64>,
    env2: &mut Vec<f64>,
) {
    let mut env: Vec<f64> = vec![0.0; env1.len() + env2.len()];
    let mut c = 0;
    for i in 0..env1.len() {
        env[c] = env1[i];
        c += 1;
    }
    for j in 0..env2.len() {
        env[c] = env2[j];
        c += 1;
    }

    let (natm, nbas) = nmol(atm, bas);

    cint1e_nuc_cart(out, shls, atm, natm as i32, bas, nbas as i32, &mut env);
}

#[no_mangle]
fn dTk(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env1: &mut Vec<f64>,
    env2: &mut Vec<f64>,
    P: &Vec<f64>,
) -> Vec<f64> {
    let (_, nbas) = nmol(atm, bas);
    let nshells = angl(bas, 0);

    let mut dS = vec![0.0; env2.len()];

    let mut buf;
    let mut dbuf;
    let mut denv;
    let mut shls = vec![0; 4];

    let mut mu;
    let mut nu;

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32; let di = CINTcgto_cart(i, &bas) as usize;
        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32; let dj = CINTcgto_cart(j, &bas) as usize;

            buf = vec![0.0; di * dj];
            dbuf = vec![0.0; di * dj];

            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    dbuf[c] = 1.0;

                    denv = vec![0.0; env2.len()];
                    dkinw(&mut buf, &mut dbuf, &mut shls, atm, bas, env1, env2, &mut denv);
                    for l in 0..env2.len() {
                        dS[l] += P[nuj * nshells + mui] * denv[l];
                    }
                    
                    dbuf[c] = 0.0;
                    c += 1;
                }
            }
            nu += dj;
        }
        mu += di;
    }
    
    return dS;
}

#[no_mangle]
pub fn dVk(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env1: &mut Vec<f64>,
    env2: &mut Vec<f64>,
    P: &Vec<f64>,
) -> Vec<f64> {
    let (_, nbas) = nmol(atm, bas);
    let nshells = angl(bas, 0);

    let mut dS = vec![0.0; env2.len()];

    let mut buf;
    let mut dbuf;
    let mut denv;
    let mut shls = vec![0; 4];

    let mut mu;
    let mut nu;

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32; let di = CINTcgto_cart(i, &bas) as usize;
        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32; let dj = CINTcgto_cart(j, &bas) as usize;

            buf = vec![0.0; di * dj];
            dbuf = vec![0.0; di * dj];

            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    dbuf[c] = 1.0;

                    denv = vec![0.0; env2.len()];
                    dnucw(&mut buf, &mut dbuf, &mut shls, atm, bas, env1, env2, &mut denv);
                    for l in 0..env2.len() {
                        dS[l] += P[nuj * nshells + mui] * denv[l];
                    }
                    
                    dbuf[c] = 0.0;
                    c += 1;
                }
            }
            nu += dj;
        }
        mu += di;
    }
    
    return dS;
}

pub fn dHk(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env1: &mut Vec<f64>,
    env2: &mut Vec<f64>,
    P: &Vec<f64>,
) -> Vec<f64> {
    let mut dT = dTk(atm, bas, env1, env2, P);
    let mut dV = dVk(atm, bas, env1, env2, P);

    let mut dH = vec![0.0; dT.len()];
    for i in 0..dT.len() {
        dH[i] = dT[i] + dV[i];
    }
    return dH;
}

fn main() -> io::Result<()> {
    let mut atm = Vec::new();
    let mut bas = Vec::new();
    let mut env = Vec::new();

    let path = librint::get_path();
    read_basis(&path, &mut atm, &mut bas, &mut env);

    // f vs df time
    let (s1, s2) = split(&mut bas);

    let mut env1: Vec<f64> = env[0..s1].to_vec();
    let mut env2: Vec<f64> = env[s1..s2].to_vec();

    let i = 0;
    let j = 0;

    let mut shls = vec![0; 4];
    shls[0] = i as i32; let di = CINTcgto_cart(i, &bas) as usize;
    shls[1] = j as i32; let dj = CINTcgto_cart(j, &bas) as usize;

    let mut buf = vec![0.0; di * dj];
    let mut dbuf = vec![0.0; di * dj];

    let mut denv = vec![0.0; s2-s1];

    let dnow = Instant::now();
    kinw(&mut buf, &mut shls, &mut atm, &mut bas, &mut env1, &mut env2);
    let delapsed_time = dnow.elapsed();

    println!("kin time: {}", delapsed_time.as_micros());

    dbuf[0] = 1.0;

    let dnow = Instant::now();
    dkinw(&mut buf, &mut dbuf, &mut shls, &mut atm, &mut bas, &mut env1, &mut env2, &mut denv);
    let delapsed_time = dnow.elapsed();

    println!("dkin time: {}", delapsed_time.as_micros());

    // f vs df time
    let (s1, s2) = split(&mut bas);

    let mut env1: Vec<f64> = env[0..s1].to_vec();
    let mut env2: Vec<f64> = env[s1..s2].to_vec();

    let i = 0;
    let j = 0;

    let mut shls = vec![0; 4];
    shls[0] = i as i32; let di = CINTcgto_cart(i, &bas) as usize;
    shls[1] = j as i32; let dj = CINTcgto_cart(j, &bas) as usize;

    let mut buf = vec![0.0; di * dj];
    let mut dbuf = vec![0.0; di * dj];

    let mut denv = vec![0.0; s2-s1];

    let dnow = Instant::now();
    nucw(&mut buf, &mut shls, &mut atm, &mut bas, &mut env1, &mut env2);
    let delapsed_time = dnow.elapsed();

    println!("nuc time: {}", delapsed_time.as_micros());

    dbuf[0] = 1.0;

    let dnow = Instant::now();
    dnucw(&mut buf, &mut dbuf, &mut shls, &mut atm, &mut bas, &mut env1, &mut env2, &mut denv);
    let delapsed_time = dnow.elapsed();

    println!("dnuc time: {}", delapsed_time.as_micros());

    // full time

    let nelec = 2;

    let P = density(&mut atm, &mut bas, &mut env, nelec, 20, 1e-6);

    let dnow = Instant::now();
    let mut dH = dHk(&mut atm, &mut bas, &mut env1, &mut env2, &P);
    let delapsed_time = dnow.elapsed();

    println!("tensor time contracted: {}", delapsed_time.as_micros());

    for i in 0..dH.len() {
        print!("{:.5} ", dH[i]);
    }
    println!();

    Ok(())
}
