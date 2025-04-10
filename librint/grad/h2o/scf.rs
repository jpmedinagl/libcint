#![allow(non_snake_case, non_upper_case_globals,unused_variables,improper_ctypes_definitions,static_mut_refs)]

use librint::utils::print_arr;
use librint::utils::read_basis;

use librint::cint_bas::CINTcgto_cart;
use librint::cint_bas::CINTcgto_spheric;
use librint::cint1e::cint1e_ovlp_cart;
use librint::cint1e::cint1e_ovlp_sph;
use librint::cint1e::cint1e_nuc_cart;
use librint::cint1e::cint1e_nuc_sph;
use librint::intor1::cint1e_kin_cart;
use librint::intor1::cint1e_kin_sph;
use librint::cint2e::cint2e_cart;
use librint::cint2e::cint2e_sph;

use faer::mat;
use faer::linalg::solvers::Eigendecomposition;
use faer::complex_native::c64;

pub const ATM_SLOTS: usize = 6;
pub const BAS_SLOTS: usize = 8;

pub const N_STO: i32 = 3;

fn matmult(
    n: usize,
    A: &mut [f64],
    B: &mut [f64],
    C: &mut [f64],
) {
    // cij = aik bkj
    C.fill(0.0);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

fn transpose(
    n: usize,
    C: &mut [f64],
    Ct: &mut [f64],
) {
    for i in 0..n {
        for j in 0..n {
            Ct[i * n + j] = C[j * n + i]
        }
    }
}

fn dcopy(
    n: usize,
    A: &mut [f64],
    Ac: &mut [f64],
) {
    for i in 0..n {
        for j in 0..n {
            Ac[i * n + j] = A[i * n + j];
        }
    }
}

fn sort(
    n: usize,
    eig: &mut [f64],
    U: &mut [f64],
) {
    for i in (0..n).rev() {
        for j in 0..i {
            if (eig[j] > eig[j+1]) {
                let tmp = eig[j];
                eig[j] = eig[j+1];
                eig[j+1] = tmp;

                for k in 0..n {
                    let tmp = U[k*n + j];
                    U[k*n + j] = U[k*n + j+1];
                    U[k*n + j+1] = tmp;
                }
            }
        }
    }
}

fn int_cart(
    natm: usize,
    nbas: usize,
    nshells: usize,
    atm: &mut [i32],
    bas: &mut [i32],
    env: &mut [f64],
    S: &mut [f64],
    T: &mut [f64],
    V: &mut [f64],
    H: &mut [f64],
    two: &mut [f64],
) {
    let mut buf: Vec<f64>;

    let mut shls: [i32; 4] = [0, 0, 0, 0];

    let mut mu;
    let mut nu;
    let mut sig;
    let mut lam;

    let mut di = 0;
    let mut dj = 0;
    let mut dk = 0;
    let mut dl = 0;

    mu = 0;
    for i in 0..nbas {
        nu = 0;
        for j in 0..nbas {
            sig = 0;

            shls[0] = i as i32;
            shls[1] = j as i32;

            di = CINTcgto_cart(i, &bas) as usize;
            dj = CINTcgto_cart(j, &bas) as usize;

            buf = vec![0.0; di * dj];

            cint1e_ovlp_cart(&mut buf, &mut shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    S[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            cint1e_kin_cart(&mut buf, &mut shls, atm, natm as i32, bas, nbas as i32, env);
            let mut c: usize = 0;
            for nuj in nu..(nu + dj) {
                for mui in mu..(mu + di) {
                    T[mui * nshells + nuj] = buf[c];
                    c += 1;
                }
            }

            cint1e_nuc_cart(&mut buf, &mut shls, atm, natm as i32, bas, nbas as i32, env);
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
                lam = 0;
                for l in 0..nbas {
                    shls[2] = k as i32;
                    shls[3] = l as i32;

                    dk = CINTcgto_cart(k, &bas) as usize;
                    dl = CINTcgto_cart(l, &bas) as usize;
                    
                    buf = vec![0.0; di * dj * dk * dl];
        
                    cint2e_cart(&mut buf, &mut shls, atm, natm as i32, bas, nbas as i32, env);
                    let mut c: usize = 0;
                    for laml in lam..(lam + dl) {
                        for sigk in sig..(sig + dk) {
                            for nuj in nu..(nu + dj) {
                                for mui in mu..(mu + di) {
                                    two[mui*nshells.pow(3) + nuj*nshells.pow(2) + sigk*nshells + laml] = buf[c];
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
}

fn find_X(
    n: usize,
    S: &mut [f64],
    X: &mut [f64],
    Xdag: &mut [f64],
) {
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

    matmult(n, &mut U, &mut lamb, X);
    transpose(n, X, Xdag);
}

fn calc_F(
    n: usize,
    P: &mut [f64],
    two: &mut [f64],
    H: &mut [f64],
    G: &mut [f64],
    F: &mut [f64],
) {
    G.fill(0.0);
    F.fill(0.0);

    for mu in 0..n {
        for nu in 0..n {
            for la in 0..n {
                for sig in 0..n {
                    // if (P[la*n + sig] != 0.0) {
                    //     println!("{:.6} {:.6}", P[la*n + sig], two[mu*n.pow(3) + nu*n.pow(2) + sig*n + la]
                    //     - 0.5 * two[mu*n.pow(3) + la*n.pow(2) + sig*n + nu]);
                    // }
                    G[mu*n + nu] += P[la*n + sig]
                        * (two[mu*n.pow(3) + nu*n.pow(2) + sig*n + la]
                        - 0.5 * two[mu*n.pow(3) + la*n.pow(2) + sig*n + nu]);
                }
            }

            F[mu*n + nu] += G[mu*n + nu] + H[mu*n + nu];
        }
    }
}

fn calc_Fprime(
    n: usize,
    F: &mut [f64],
    X: &mut [f64],
    Xdag: &mut [f64],
    mut Fprime: &mut [f64],
) {
    let mut inter = vec![0.0; n * n];
    matmult(n, Xdag, F, &mut inter);
    matmult(n, &mut inter, X, &mut Fprime);
}

fn diag_F(
    n: usize,
    Fprime: &mut [f64],
    X: &mut [f64],
    C: &mut [f64],
    epsilon: &mut [f64],
) {
    // diagonalize Fprime
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

    let mut Udag = vec![0.0; n * n];
    transpose(n, &mut U, &mut Udag);

    let mut inter = vec![0.0; n * n];
    matmult(n, &mut Udag, Fprime, &mut inter); 

    let mut f = vec![0.0; n * n];
    matmult(n, &mut inter, &mut U, &mut f);

    let mut Cprime = vec![0.0; n * n];
    dcopy(n, &mut U, &mut Cprime);

    dcopy(n, &mut f, epsilon);

    matmult(n, X, &mut Cprime, C);
}

fn calc_P(
    n: usize,
    nelec: usize,
    C: &mut [f64],
) -> Vec<f64> {
    let mut P = vec![0.0; n * n];
    for mu in 0..n {
        for nu in 0..n {
            for i in 0..(nelec/2) {
                P[mu * n + nu] += 2.0 * C[mu * n + i] * C[nu * n + i];
            }
        }
    }
    return P;
}

fn f_delta(
    n: usize,
    P: &mut [f64],
    Pold: &mut [f64],
) -> f64 {
    let mut delta: f64 = 0.0;
    for mu in 0..n {
        for nu in 0..n {
            delta += (P[mu * n + nu] - Pold[mu * n + nu]).powf(2.0);
        }
    }
    delta = delta.powf(0.5) / 2.0;
    return delta;
}

fn norm(
    atm: &mut [i32],
    env: &mut [f64],
    i: usize,
    j: usize,
) -> f64 {
    let xi: f64 = env[(atm[i*6 + 1]) as usize];
    let xj: f64 = env[(atm[j*6 + 1]) as usize];

    let yi: f64 = env[(atm[i*6 + 1] + 1) as usize];
    let yj: f64 = env[(atm[j*6 + 1] + 1) as usize];

    let zi: f64 = env[(atm[i*6 + 1] + 2) as usize];
    let zj: f64 = env[(atm[j*6 + 1] + 2) as usize];

    return ((xi - xj).powf(2.0) + (yi - yj).powf(2.0) + (zi - zj).powf(2.0)).powf(0.5);
}

fn RHF(
    natm: usize,
    nbas: usize,
    nelec: usize,
    nshells: usize,
    atm: &mut [i32],
    bas: &mut [i32],
    env: &mut [f64],
    imax: i32,
    conv: f64,
) -> Vec<f64> {
    let mut S = vec![0.0; nshells * nshells];
    let mut T = vec![0.0; nshells * nshells];
    let mut V = vec![0.0; nshells * nshells];
    let mut H = vec![0.0; nshells * nshells];
    let mut two = vec![0.0; nshells * nshells * nshells * nshells];
    int_cart(natm, nbas, nshells, atm, bas, env, &mut S, &mut T, &mut V, &mut H, &mut two);

    let mut X = vec![0.0; nshells * nshells];
    let mut Xdag = vec![0.0; nshells * nshells];
    find_X(nshells, &mut S, &mut X, &mut Xdag);

    let mut P = vec![0.0; nshells * nshells];
    for i in 0..nshells {
        for j in 0..nshells {
            if i == j {
                P[i*nshells + j] = 1.0;
            }
        }
    }

    let mut i: i32 = 0;
    let mut delta: f64 = 1.0;

    let mut Pold = vec![0.0; nshells * nshells];
    let mut G = vec![0.0; nshells * nshells];
    let mut F = vec![0.0; nshells * nshells];
    let mut Fprime = vec![0.0; nshells * nshells];
    let mut C = vec![0.0; nshells * nshells];
    let mut epsilon = vec![0.0; nshells * nshells];

    while delta > conv && i < imax {
        dcopy(nshells, &mut P, &mut Pold);

        calc_F(nshells, &mut P, &mut two, &mut H, &mut G, &mut F);
        calc_Fprime(nshells, &mut F, &mut X, &mut Xdag, &mut Fprime);
        diag_F(nshells, &mut Fprime, &mut X, &mut C, &mut epsilon);

        P = calc_P(nshells, nelec, &mut C);
        delta = f_delta(nshells, &mut P, &mut Pold);
        i += 1;

        // println!("i: {}", i);
        // print_arr(nshells, 2, &mut P);
    }

    println!("Conv P: ({}/{})", i, imax);
    print_arr(nshells, 2, &mut P);

    if delta > conv && i == imax {
        println!("did not converge");
    }

    return P;
}

fn energy(
    natm: usize,
    nbas: usize,
    nelec: usize,
    nshells: usize,
    atm: &mut [i32],
    bas: &mut [i32],
    env: &mut [f64],
    mut P: Vec<f64>,
) -> f64 {
    let mut S = vec![0.0; nshells * nshells];
    let mut T = vec![0.0; nshells * nshells];
    let mut V = vec![0.0; nshells * nshells];
    let mut H = vec![0.0; nshells * nshells];
    let mut two = vec![0.0; nshells * nshells * nshells * nshells];
    int_cart(natm, nbas, nshells, atm, bas, env, &mut S, &mut T, &mut V, &mut H, &mut two);

    let mut G = vec![0.0; nshells * nshells];
    let mut F = vec![0.0; nshells * nshells];
    calc_F(nshells, &mut P, &mut two, &mut H, &mut G, &mut F);

    let mut E0: f64 = 0.0;
    for mu in 0..nshells {
        for nu in 0..nshells {
            E0 += 0.5 * P[mu*nshells + nu] * (H[mu*nshells + nu]  + F[mu*nshells + nu]);
        }
    }

    let mut Enuc: f64 = 0.0;
    for i in 0..natm {
        for j in 0..natm {
            if i > j {
                Enuc += (atm[i*6 + 0] * atm[j*6 + 0]) as f64 / (norm(atm, env, i, j));
            }
        }
    }

    return E0 + Enuc;
}

fn main() {
    let mut atm: Vec<i32> = Vec::new();
    let mut bas: Vec<i32> = Vec::new();
    let mut env: Vec<f64> = Vec::new();

    let path = librint::get_path();
    read_basis(&path, &mut atm, &mut bas, &mut env);

    const natm: usize = 3;
    const nelec: usize = 10;

    const nbas: usize = 5;
    const nshells: usize = 7;

    const imax: i32 = 20;
    const conv: f64 = 0.000001;

    let P = RHF(natm, nbas, nelec, nshells, &mut atm, &mut bas, &mut env, imax, conv);

    let Etot: f64 = energy(natm, nbas, nelec, nshells, &mut atm, &mut bas, &mut env, P);
    println!("Etot: {}", Etot);
}
