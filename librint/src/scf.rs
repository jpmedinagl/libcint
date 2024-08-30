#![allow(non_snake_case, non_upper_case_globals, non_camel_case_types)]

use crate::cint1e::{cint1e_nuc_cart, cint1e_nuc_sph, cint1e_ovlp_cart, cint1e_ovlp_sph};
use crate::cint2e::{cint2e_cart, cint2e_sph};
use crate::cint_bas::{CINTcgto_cart, CINTcgto_spheric};
use crate::intor1::{cint1e_kin_cart, cint1e_kin_sph};

use crate::linalg::{dcopya, matmult, sort, transpose};

use faer::{complex_native::c64, linalg::solvers::Eigendecomposition, mat};

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

type inte = fn(
    buf: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32;

type cgto = fn(bas_id: usize, bas: &[i32]) -> i32;

#[no_mangle]
pub fn nmol(atm: &Vec<i32>, bas: &Vec<i32>) -> (usize, usize) {
    let natm: usize = atm.len() / ATM_SLOTS;
    let nbas: usize = bas.len() / BAS_SLOTS;
    return (natm, nbas);
}

#[no_mangle]
pub fn angl(bas: &Vec<i32>, coord: i32) -> usize {
    let mut nshells: usize = 0;
    for i in (0..bas.len()).step_by(BAS_SLOTS) {
        let l = bas[i + 1] as usize;
        if coord == 0 {
            nshells += (l + 1) * (l + 2) / 2;
        } else if coord == 1 {
            nshells += 2 * l + 1;
        }
    }
    return nshells;
}

#[no_mangle]
pub fn integral1e(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env: &mut Vec<f64>,
    coord: i32,
    typec: i32,
) -> Vec<f64> {
    let (natm, nbas) = nmol(atm, bas);
    let nshells = angl(bas, coord);

    let mut R = vec![0.0; nshells * nshells];

    let mut buf: Vec<f64>;
    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut mu;
    let mut nu;

    let mut di;
    let mut dj;

    let intcgto: cgto;
    let func: inte;

    if coord == 0 {
        intcgto = CINTcgto_cart;

        if typec == 0 {
            func = cint1e_ovlp_cart;
        } else if typec == 1 {
            func = cint1e_kin_cart;
        } else if typec == 2 {
            func = cint1e_nuc_cart;
        } else {
            std::process::exit(1);
        }
    } else if coord == 1 {
        intcgto = CINTcgto_spheric;

        if typec == 0 {
            func = cint1e_ovlp_sph;
        } else if typec == 1 {
            func = cint1e_kin_sph;
        } else if typec == 2 {
            func = cint1e_nuc_sph;
        } else {
            std::process::exit(1);
        }
    } else {
        std::process::exit(1);
    }

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32;
        di = intcgto(i, &bas) as usize;

        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32;
            dj = intcgto(j, &bas) as usize;

            buf = vec![0.0; di * dj];

            func(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    R[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            nu += dj;
        }
        mu += di;
    }

    return R;
}

#[no_mangle]
pub fn integral2e(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env: &mut Vec<f64>,
    coord: i32,
) -> Vec<f64> {
    let (natm, nbas) = nmol(atm, bas);
    let nshells = angl(bas, coord);

    let mut R = vec![0.0; nshells * nshells * nshells * nshells];

    let mut buf: Vec<f64>;
    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut mu;
    let mut nu;
    let mut sig;
    let mut lam;

    let mut di;
    let mut dj;
    let mut dk;
    let mut dl;

    let intcgto: cgto;
    let func: inte;

    if coord == 0 {
        intcgto = CINTcgto_cart;
        func = cint2e_cart;
    } else if coord == 1 {
        intcgto = CINTcgto_spheric;
        func = cint2e_sph;
    } else {
        std::process::exit(1);
    }

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32;
        di = intcgto(i, &bas) as usize;

        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32;
            dj = intcgto(j, &bas) as usize;

            sig = 0;
            for k in 0..nbas {
                shls[2] = k as i32;
                dk = intcgto(k, &bas) as usize;

                lam = 0;
                for l in 0..nbas {
                    shls[3] = l as i32;
                    dl = intcgto(l, &bas) as usize;

                    buf = vec![0.0; di * dj * dk * dl];

                    func(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
                    let mut c: usize = 0;
                    for laml in lam..(lam + dl) {
                        for sigk in sig..(sig + dk) {
                            for nuj in nu..(nu + dj) {
                                for mui in mu..(mu + di) {
                                    R[mui * nshells.pow(3)
                                        + nuj * nshells.pow(2)
                                        + sigk * nshells
                                        + laml] = buf[c];
                                    c += 1;
                                }
                            }
                        }
                    }

                    lam += dl;
                }
                sig += dk;
            }
            nu += dj;
        }
        mu += di;
    }

    return R;
}

fn integrals(
    natm: usize,
    nbas: usize,
    nshells: usize,
    atm: &mut [i32],
    bas: &mut [i32],
    env: &mut [f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut S = vec![0.0; nshells * nshells];
    let mut H = vec![0.0; nshells * nshells];
    let mut two = vec![0.0; nshells * nshells * nshells * nshells];

    let mut buf: Vec<f64>;

    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut mu;
    let mut nu;
    let mut sig;
    let mut lam;

    let mut di;
    let mut dj;
    let mut dk;
    let mut dl;

    let mut T = vec![0.0; nshells * nshells];
    let mut V = vec![0.0; nshells * nshells];

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32;
        di = CINTcgto_cart(i, &bas) as usize;

        nu = 0;
        for j in 0..nbas {
            sig = 0;

            shls[1] = j as i32;
            dj = CINTcgto_cart(j, &bas) as usize;

            buf = vec![0.0; di * dj];

            cint1e_ovlp_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    S[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            cint1e_kin_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    T[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            cint1e_nuc_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    V[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    H[mui * nshells + nuj] = T[mui * nshells + nuj] + V[mui * nshells + nuj];
                }
            }

            for k in 0..nbas {
                shls[2] = k as i32;
                dk = CINTcgto_cart(k, &bas) as usize;

                lam = 0;
                for l in 0..nbas {
                    shls[3] = l as i32;
                    dl = CINTcgto_cart(l, &bas) as usize;

                    buf = vec![0.0; di * dj * dk * dl];

                    cint2e_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
                    let mut c: usize = 0;
                    for laml in lam..(lam + dl) {
                        for sigk in sig..(sig + dk) {
                            for nuj in nu..(nu + dj) {
                                for mui in mu..(mu + di) {
                                    two[mui * nshells.pow(3)
                                        + nuj * nshells.pow(2)
                                        + sigk * nshells
                                        + laml] = buf[c];
                                    c += 1;
                                }
                            }
                        }
                    }

                    lam += dl;
                }
                sig += dk;
            }
            nu += dj;
        }
        mu += di;
    }

    return (S, H, two);
}

fn find_X(n: usize, S: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let s_mat = mat::from_column_major_slice::<f64>(&S, n, n);
    let eig_decomp: Eigendecomposition<c64> = Eigendecomposition::new_from_real(s_mat);
    let eigenvalues = eig_decomp.s();
    let eigenvectors = eig_decomp.u();

    let mut U = vec![0.0; n * n];

    for i in 0..eigenvectors.nrows() {
        for j in 0..eigenvectors.ncols() {
            U[i * n + j] = eigenvectors.read(i, j).re();
        }
    }

    let eign = eigenvalues.column_vector();
    let mut eig = vec![0.0; n];
    for i in 0..n {
        eig[i] = eign.read(i).re();
    }

    sort(n, &mut eig, &mut U);

    let mut lamb = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                let li: f64 = eig[i];
                lamb[i * n + j] = li.powf(-0.5);
            }
        }
    }

    let X = matmult(n, &U, &lamb);
    let Xdag = transpose(n, &X);

    return (X, Xdag);
}

pub fn calc_F(n: usize, P: &[f64], two: &[f64], H: &[f64]) -> Vec<f64> {
    let mut G = vec![0.0; n * n];
    let mut F = vec![0.0; n * n];

    for mu in 0..n {
        for nu in 0..n {
            for la in 0..n {
                for sig in 0..n {
                    G[mu * n + nu] += P[la * n + sig]
                        * (two[mu * n.pow(3) + nu * n.pow(2) + sig * n + la]
                            - 0.5 * two[mu * n.pow(3) + la * n.pow(2) + sig * n + nu]);
                }
            }

            F[mu * n + nu] += G[mu * n + nu] + H[mu * n + nu];
        }
    }

    return F;
}

fn calc_Fprime(n: usize, F: &[f64], X: &[f64], Xdag: &[f64]) -> Vec<f64> {
    let inter = matmult(n, &Xdag, &F);
    let Fprime = matmult(n, &inter, X);
    return Fprime;
}

fn diag_F(n: usize, Fprime: &[f64], X: &[f64]) -> Vec<f64> {
    let fprime_mat = mat::from_column_major_slice::<f64>(&Fprime, n, n);
    let eig_decomp: Eigendecomposition<c64> = Eigendecomposition::new_from_real(fprime_mat);
    let eigenvalues = eig_decomp.s();
    let eigenvectors = eig_decomp.u();

    let mut U = vec![0.0; n * n];
    for i in 0..eigenvectors.nrows() {
        for j in 0..eigenvectors.ncols() {
            U[i * n + j] = eigenvectors.read(i, j).re();
        }
    }

    let eign = eigenvalues.column_vector();
    let mut eig = vec![0.0; n];
    for i in 0..n {
        eig[i] = eign.read(i).re();
    }

    sort(n, &mut eig, &mut U);

    let Udag = transpose(n, &U);
    let inter = matmult(n, &Udag, Fprime);
    let f = matmult(n, &inter, &U);
    let Cprime = dcopya(n, &U);

    let _epsilon = dcopya(n, &f);
    let C = matmult(n, X, &Cprime);

    return C;
}

fn calc_P(n: usize, nelec: usize, C: &mut [f64]) -> Vec<f64> {
    let mut P = vec![0.0; n * n];
    for mu in 0..n {
        for nu in 0..n {
            for i in 0..(nelec / 2) {
                P[mu * n + nu] += 2.0 * C[mu * n + i] * C[nu * n + i];
            }
        }
    }
    return P;
}

fn f_delta(n: usize, P: &mut [f64], Pold: &mut [f64]) -> f64 {
    let mut delta: f64 = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            delta += (P[mu * n + nu] - Pold[mu * n + nu]).powf(2.0);
        }
    }
    delta = delta.powf(0.5) / 2.0;
    return delta;
}

pub fn norm(atm: &mut [i32], env: &mut [f64], i: usize, j: usize) -> f64 {
    let xi: f64 = env[(atm[i * 6 + 1]) as usize];
    let xj: f64 = env[(atm[j * 6 + 1]) as usize];

    let yi: f64 = env[(atm[i * 6 + 1] + 1) as usize];
    let yj: f64 = env[(atm[j * 6 + 1] + 1) as usize];

    let zi: f64 = env[(atm[i * 6 + 1] + 2) as usize];
    let zj: f64 = env[(atm[j * 6 + 1] + 2) as usize];

    return ((xi - xj).powf(2.0) + (yi - yj).powf(2.0) + (zi - zj).powf(2.0)).powf(0.5);
}

#[no_mangle]
pub fn density(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env: &mut Vec<f64>,
    nelec: usize,
    imax: i32,
    conv: f64,
) -> Vec<f64> {
    let (natm, nbas) = nmol(atm, bas);
    let nshells = angl(bas, 0);

    let (S, H, two) = integrals(natm, nbas, nshells, atm, bas, env);
    let (X, Xdag) = find_X(nshells, &S);

    let mut P = vec![0.0; nshells * nshells];
    for i in 0..nshells {
        for j in 0..nshells {
            if i == j {
                P[i * nshells + j] = 1.0;
            }
        }
    }

    let mut i: i32 = 0;
    let mut delta: f64 = 1.0;

    let mut Pold;
    let mut F;
    let mut Fprime;
    let mut C;

    while delta > conv && i < imax {
        Pold = dcopya(nshells, &P);

        F = calc_F(nshells, &P, &two, &H);
        Fprime = calc_Fprime(nshells, &F, &X, &Xdag);
        C = diag_F(nshells, &Fprime, &X);

        P = calc_P(nshells, nelec, &mut C);
        delta = f_delta(nshells, &mut P, &mut Pold);
        i += 1;
    }

    if delta > conv && i == imax {
        println!("did not converge");
        for i in 0..nshells {
            for j in 0..nshells {
                P[i * nshells + j] = 0.0;
            }
        }
    }

    return P;
}

#[no_mangle]
pub fn energy(atm: &mut Vec<i32>, bas: &mut Vec<i32>, env: &mut Vec<f64>, P: &mut Vec<f64>) -> f64 {
    let (natm, nbas) = nmol(atm, bas);
    let nshells = angl(bas, 0);

    let (_, H, two) = integrals(natm, nbas, nshells, atm, bas, env);
    let F = calc_F(nshells, P, &two, &H);

    let mut E0: f64 = 0.0;
    for mu in 0..nshells {
        for nu in 0..nshells {
            E0 += 0.5 * P[mu * nshells + nu] * (H[mu * nshells + nu] + F[mu * nshells + nu]);
        }
    }

    let mut Enuc: f64 = 0.0;
    for i in 0..natm {
        for j in 0..natm {
            if i > j {
                Enuc += (atm[i * 6 + 0] * atm[j * 6 + 0]) as f64 / (norm(atm, env, i, j));
            }
        }
    }

    return E0 + Enuc;
}

#[no_mangle]
pub fn energyfast(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env: &mut Vec<f64>,
    P: &mut Vec<f64>,
) -> f64 {
    let (natm, nbas) = nmol(atm, bas);
    let nshells = angl(bas, 0);

    let mut buf;
    let mut shls = [0; 4];

    let mut mu;
    let mut nu;
    let mut sig;
    let mut lam;

    let mut E0: f64 = 0.0;

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32;
        let di = CINTcgto_cart(i, &bas) as usize;
        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32;
            let dj = CINTcgto_cart(j, &bas) as usize;

            buf = vec![0.0; di * dj];

            cint1e_kin_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    E0 += P[mui * nshells + nuj] * buf[c];
                    c += 1;
                }
            }

            cint1e_nuc_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    E0 += P[mui * nshells + nuj] * buf[c];
                    c += 1;
                }
            }
            nu += dj;
        }
        mu += di;
    }

    mu = 0;
    for i in 0..nbas {
        shls[0] = i as i32;
        let di = CINTcgto_cart(i, &bas) as usize;
        nu = 0;
        for j in 0..nbas {
            shls[1] = j as i32;
            let dj = CINTcgto_cart(j, &bas) as usize;
            sig = 0;
            for k in 0..nbas {
                shls[2] = k as i32;
                let dk = CINTcgto_cart(k, &bas) as usize;
                lam = 0;
                for l in 0..nbas {
                    shls[3] = l as i32;
                    let dl = CINTcgto_cart(l, &bas) as usize;

                    buf = vec![0.0; di * dj * dk * dl];

                    cint2e_cart(&mut buf, shls, atm, natm as i32, bas, nbas as i32, env);
                    let mut c: usize = 0;
                    for laml in lam..(lam + dl) {
                        for sigk in sig..(sig + dk) {
                            for nuj in nu..(nu + dj) {
                                for mui in mu..(mu + di) {
                                    E0 += 0.5
                                        * (P[mui * nshells + nuj] * P[sigk * nshells + laml]
                                            - 0.5
                                                * P[mui * nshells + sigk]
                                                * P[nuj * nshells + laml])
                                        * buf[c];
                                    c += 1;
                                }
                            }
                        }
                    }
                    lam += dl;
                }
                sig += dk;
            }
            nu += dj;
        }
        mu += di;
    }

    let mut Enuc: f64 = 0.0;
    for i in 0..natm {
        for j in 0..natm {
            if i > j {
                Enuc += (atm[i * 6 + 0] * atm[j * 6 + 0]) as f64 / (norm(atm, env, i, j));
            }
        }
    }

    return E0 + Enuc;
}

#[no_mangle]
pub fn scf(
    atm: &mut Vec<i32>,
    bas: &mut Vec<i32>,
    env: &mut Vec<f64>,
    nelec: usize,
    imax: i32,
    conv: f64,
) -> f64 {
    let mut P = density(atm, bas, env, nelec, imax, conv);
    let Etot = energyfast(atm, bas, env, &mut P);
    return Etot;
}
