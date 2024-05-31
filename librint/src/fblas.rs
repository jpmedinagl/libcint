#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]
#[no_mangle]
pub unsafe extern "C" fn CINTdset0(mut n: libc::c_int, mut x: *mut libc::c_double) {
    let mut i: libc::c_int = 0;
    i = 0 as libc::c_int;
    while i < n {
        *x.offset(i as isize) = 0 as libc::c_int as libc::c_double;
        i += 1;
        i;
    }
}
#[no_mangle]
pub unsafe extern "C" fn CINTdaxpy2v(
    mut n: libc::c_int,
    mut a: libc::c_double,
    mut x: *mut libc::c_double,
    mut y: *mut libc::c_double,
    mut v: *mut libc::c_double,
) {
    let mut i: libc::c_int = 0;
    i = 0 as libc::c_int;
    while i < n {
        *v.offset(i as isize) = a * *x.offset(i as isize) + *y.offset(i as isize);
        i += 1;
        i;
    }
}
#[no_mangle]
pub unsafe extern "C" fn CINTdmat_transpose(
    mut a_t: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut m: libc::c_int,
    mut n: libc::c_int,
) {
    let mut i: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    j = 0 as libc::c_int;
    while j < n - 3 as libc::c_int {
        i = 0 as libc::c_int;
        while i < m {
            *a_t
                .offset(
                    ((j + 0 as libc::c_int) * m + i) as isize,
                ) = *a.offset((i * n + j + 0 as libc::c_int) as isize);
            *a_t
                .offset(
                    ((j + 1 as libc::c_int) * m + i) as isize,
                ) = *a.offset((i * n + j + 1 as libc::c_int) as isize);
            *a_t
                .offset(
                    ((j + 2 as libc::c_int) * m + i) as isize,
                ) = *a.offset((i * n + j + 2 as libc::c_int) as isize);
            *a_t
                .offset(
                    ((j + 3 as libc::c_int) * m + i) as isize,
                ) = *a.offset((i * n + j + 3 as libc::c_int) as isize);
            i += 1;
            i;
        }
        j += 4 as libc::c_int;
    }
    match n - j {
        1 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t.offset((j * m + i) as isize) = *a.offset((i * n + j) as isize);
                i += 1;
                i;
            }
        }
        2 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t
                    .offset(
                        ((j + 0 as libc::c_int) * m + i) as isize,
                    ) = *a.offset((i * n + j + 0 as libc::c_int) as isize);
                *a_t
                    .offset(
                        ((j + 1 as libc::c_int) * m + i) as isize,
                    ) = *a.offset((i * n + j + 1 as libc::c_int) as isize);
                i += 1;
                i;
            }
        }
        3 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t
                    .offset(
                        ((j + 0 as libc::c_int) * m + i) as isize,
                    ) = *a.offset((i * n + j + 0 as libc::c_int) as isize);
                *a_t
                    .offset(
                        ((j + 1 as libc::c_int) * m + i) as isize,
                    ) = *a.offset((i * n + j + 1 as libc::c_int) as isize);
                *a_t
                    .offset(
                        ((j + 2 as libc::c_int) * m + i) as isize,
                    ) = *a.offset((i * n + j + 2 as libc::c_int) as isize);
                i += 1;
                i;
            }
        }
        _ => {}
    };
}
#[no_mangle]
pub unsafe extern "C" fn CINTdplus_transpose(
    mut a_t: *mut libc::c_double,
    mut a: *mut libc::c_double,
    mut m: libc::c_int,
    mut n: libc::c_int,
) {
    let mut i: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut k: libc::c_int = 0;
    j = 0 as libc::c_int;
    while j < n - 3 as libc::c_int {
        i = 0 as libc::c_int;
        while i < m {
            *a_t.offset(((j + 0 as libc::c_int) * m + i) as isize)
                += *a.offset((i * n + j + 0 as libc::c_int) as isize);
            *a_t.offset(((j + 1 as libc::c_int) * m + i) as isize)
                += *a.offset((i * n + j + 1 as libc::c_int) as isize);
            *a_t.offset(((j + 2 as libc::c_int) * m + i) as isize)
                += *a.offset((i * n + j + 2 as libc::c_int) as isize);
            *a_t.offset(((j + 3 as libc::c_int) * m + i) as isize)
                += *a.offset((i * n + j + 3 as libc::c_int) as isize);
            i += 1;
            i;
        }
        j += 4 as libc::c_int;
    }
    match n - j {
        1 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t.offset((j * m + i) as isize) += *a.offset((i * n + j) as isize);
                i += 1;
                i;
            }
        }
        2 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t.offset(((j + 0 as libc::c_int) * m + i) as isize)
                    += *a.offset((i * n + j + 0 as libc::c_int) as isize);
                *a_t.offset(((j + 1 as libc::c_int) * m + i) as isize)
                    += *a.offset((i * n + j + 1 as libc::c_int) as isize);
                i += 1;
                i;
            }
        }
        3 => {
            i = 0 as libc::c_int;
            while i < m {
                *a_t.offset(((j + 0 as libc::c_int) * m + i) as isize)
                    += *a.offset((i * n + j + 0 as libc::c_int) as isize);
                *a_t.offset(((j + 1 as libc::c_int) * m + i) as isize)
                    += *a.offset((i * n + j + 1 as libc::c_int) as isize);
                *a_t.offset(((j + 2 as libc::c_int) * m + i) as isize)
                    += *a.offset((i * n + j + 2 as libc::c_int) as isize);
                i += 1;
                i;
            }
        }
        _ => {}
    };
}
#[no_mangle]
pub unsafe extern "C" fn CINTdgemm_NN1(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut k: libc::c_int,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut c: *mut libc::c_double,
    mut ldc: libc::c_int,
) {
    let mut i: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut kp: libc::c_int = 0;
    let mut bi: libc::c_double = 0.;
    j = 0 as libc::c_int;
    while j < n {
        i = 0 as libc::c_int;
        while i < m {
            *c.offset((i + ldc * j) as isize) = 0 as libc::c_int as libc::c_double;
            i += 1;
            i;
        }
        kp = 0 as libc::c_int;
        while kp < k {
            bi = *b.offset((kp + k * j) as isize);
            i = 0 as libc::c_int;
            while i < m {
                *c.offset((i + ldc * j) as isize)
                    += *a.offset((i + m * kp) as isize) * bi;
                i += 1;
                i;
            }
            kp += 1;
            kp;
        }
        j += 1;
        j;
    }
}
#[no_mangle]
pub unsafe extern "C" fn CINTdgemm_NN(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut k: libc::c_int,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut c: *mut libc::c_double,
) {
    CINTdgemm_NN1(m, n, k, a, b, c, m);
}
#[no_mangle]
pub unsafe extern "C" fn CINTdgemm_TN(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut k: libc::c_int,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut c: *mut libc::c_double,
) {
    let mut i: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut kp: libc::c_int = 0;
    let mut ci: libc::c_double = 0.;
    j = 0 as libc::c_int;
    while j < n {
        i = 0 as libc::c_int;
        while i < m {
            ci = 0 as libc::c_int as libc::c_double;
            kp = 0 as libc::c_int;
            while kp < k {
                ci
                    += *a.offset((kp + k * i) as isize)
                        * *b.offset((kp + k * j) as isize);
                kp += 1;
                kp;
            }
            *c.offset((i + m * j) as isize) = ci;
            i += 1;
            i;
        }
        j += 1;
        j;
    }
}
#[no_mangle]
pub unsafe extern "C" fn CINTdgemm_NT(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut k: libc::c_int,
    mut a: *mut libc::c_double,
    mut b: *mut libc::c_double,
    mut c: *mut libc::c_double,
) {
    let mut i: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut kp: libc::c_int = 0;
    let mut bi: libc::c_double = 0.;
    j = 0 as libc::c_int;
    while j < n {
        i = 0 as libc::c_int;
        while i < m {
            *c.offset((i + m * j) as isize) = 0 as libc::c_int as libc::c_double;
            i += 1;
            i;
        }
        kp = 0 as libc::c_int;
        while kp < k {
            bi = *b.offset((j + n * kp) as isize);
            i = 0 as libc::c_int;
            while i < m {
                *c.offset((i + m * j) as isize) += *a.offset((i + m * kp) as isize) * bi;
                i += 1;
                i;
            }
            kp += 1;
            kp;
        }
        j += 1;
        j;
    }
}