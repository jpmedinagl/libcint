#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]

use crate::optimizer::CINTOpt_log_max_pgto_coeff;
use crate::optimizer::CINTOpt_non0coeff_byshell;
use crate::optimizer::CINTset_pairdata;
use crate::g1e::CINTg1e_index_xyz;
use crate::g1e::CINTprim_to_ctr_0;
use crate::g1e::CINTprim_to_ctr_1;
use crate::g1e_grids::CINTinit_int1e_grids_EnvVars;
use crate::g1e_grids::CINTg0_1e_grids;
use crate::g1e_grids::CINTg0_1e_grids;
use crate::g1e_grids::CINTgout1e_grids;
use crate::cart2sph::c2s_sph_1e_grids;
use crate::cart2sph::c2s_cart_1e_grids;
use crate::cart2sph::c2s_grids_dset0;

use crate::cint::PairData;
use crate::cint::CINTOPt;
use crate::cint::CINTEnvVars;

pub type size_t = libc::c_ulong;
pub type uintptr_t = libc::c_ulong;

extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
}

#[no_mangle]
pub unsafe extern "C" fn CINT1e_grids_loop(
    mut gctr: *mut libc::c_double,
    mut envs: *mut CINTEnvVars,
    mut cache: *mut libc::c_double,
) -> libc::c_int {
    let mut shls: *mut libc::c_int = (*envs).shls;
    let mut bas: *mut libc::c_int = (*envs).bas;
    let mut env: *mut libc::c_double = (*envs).env;
    let mut i_sh: libc::c_int = *shls.offset(0 as libc::c_int as isize);
    let mut j_sh: libc::c_int = *shls.offset(1 as libc::c_int as isize);
    let mut i_ctr: libc::c_int = (*envs).x_ctr[0 as libc::c_int as usize];
    let mut j_ctr: libc::c_int = (*envs).x_ctr[1 as libc::c_int as usize];
    let mut i_prim: libc::c_int = *bas
        .offset((8 as libc::c_int * i_sh + 2 as libc::c_int) as isize);
    let mut j_prim: libc::c_int = *bas
        .offset((8 as libc::c_int * j_sh + 2 as libc::c_int) as isize);
    let mut nf: libc::c_int = (*envs).nf;
    let mut n_comp: libc::c_int = (*envs).ncomp_e1 * (*envs).ncomp_tensor;
    let mut ngrids: libc::c_int = (*envs).c2rust_unnamed_0.ngrids;
    let mut grids: *mut libc::c_double = (*envs).c2rust_unnamed_1.grids;
    let mut ai: *mut libc::c_double = env
        .offset(
            *bas.offset((8 as libc::c_int * i_sh + 5 as libc::c_int) as isize) as isize,
        );
    let mut aj: *mut libc::c_double = env
        .offset(
            *bas.offset((8 as libc::c_int * j_sh + 5 as libc::c_int) as isize) as isize,
        );
    let mut ci: *mut libc::c_double = env
        .offset(
            *bas.offset((8 as libc::c_int * i_sh + 6 as libc::c_int) as isize) as isize,
        );
    let mut cj: *mut libc::c_double = env
        .offset(
            *bas.offset((8 as libc::c_int * j_sh + 6 as libc::c_int) as isize) as isize,
        );
    let mut expcutoff: libc::c_double = (*envs).expcutoff;
    let mut log_maxci: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut log_maxcj: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut pdata_base: *mut PairData = 0 as *mut PairData;
    let mut pdata_ij: *mut PairData = 0 as *mut PairData;
    log_maxci = ((cache as uintptr_t).wrapping_add(7 as libc::c_int as libc::c_ulong)
        & (8 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut libc::c_double;
    cache = log_maxci.offset((i_prim + j_prim) as isize);
    pdata_base = ((cache as uintptr_t).wrapping_add(7 as libc::c_int as libc::c_ulong)
        & (8 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut PairData;
    cache = pdata_base.offset((i_prim * j_prim) as isize) as *mut libc::c_double;
    log_maxcj = log_maxci.offset(i_prim as isize);
    CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    if CINTset_pairdata(
        pdata_base,
        ai,
        aj,
        (*envs).ri,
        (*envs).rj,
        log_maxci,
        log_maxcj,
        (*envs).li_ceil,
        (*envs).lj_ceil,
        i_prim,
        j_prim,
        (*envs).rirj[0 as libc::c_int as usize] * (*envs).rirj[0 as libc::c_int as usize]
            + (*envs).rirj[1 as libc::c_int as usize]
                * (*envs).rirj[1 as libc::c_int as usize]
            + (*envs).rirj[2 as libc::c_int as usize]
                * (*envs).rirj[2 as libc::c_int as usize],
        expcutoff,
        env,
    ) != 0
    {
        return 0 as libc::c_int;
    }
    let mut fac1i: libc::c_double = 0.;
    let mut fac1j: libc::c_double = 0.;
    let mut expij: libc::c_double = 0.;
    let mut cutoff: libc::c_double = 0.;
    let mut rij: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut ip: libc::c_int = 0;
    let mut jp: libc::c_int = 0;
    let mut i: libc::c_int = 0;
    let mut grids_offset: libc::c_int = 0;
    let mut bgrids: libc::c_int = 0;
    let mut empty: [libc::c_int; 4] = [
        1 as libc::c_int,
        1 as libc::c_int,
        1 as libc::c_int,
        1 as libc::c_int,
    ];
    let mut gempty: *mut libc::c_int = empty
        .as_mut_ptr()
        .offset(0 as libc::c_int as isize);
    let mut iempty: *mut libc::c_int = empty
        .as_mut_ptr()
        .offset(1 as libc::c_int as isize);
    let mut jempty: *mut libc::c_int = empty
        .as_mut_ptr()
        .offset(2 as libc::c_int as isize);
    let mut all_empty: libc::c_int = 1 as libc::c_int;
    let mut idx: *mut libc::c_int = 0 as *mut libc::c_int;
    idx = ((cache as uintptr_t).wrapping_add(7 as libc::c_int as libc::c_ulong)
        & (8 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut libc::c_int;
    cache = idx.offset((nf * 3 as libc::c_int) as isize) as *mut libc::c_double;
    CINTg1e_index_xyz(idx, envs);
    let mut non0ctri: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut non0ctrj: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut non0idxi: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut non0idxj: *mut libc::c_int = 0 as *mut libc::c_int;
    non0ctri = ((cache as uintptr_t).wrapping_add(7 as libc::c_int as libc::c_ulong)
        & (8 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut libc::c_int;
    cache = non0ctri.offset((i_prim + j_prim + i_prim * i_ctr + j_prim * j_ctr) as isize)
        as *mut libc::c_double;
    non0ctrj = non0ctri.offset(i_prim as isize);
    non0idxi = non0ctrj.offset(j_prim as isize);
    non0idxj = non0idxi.offset((i_prim * i_ctr) as isize);
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    let nc: libc::c_int = i_ctr * j_ctr;
    let leng: libc::c_int = (*envs).g_size * 3 as libc::c_int
        * (((1 as libc::c_int) << (*envs).gbits) + 1 as libc::c_int);
    let lenj: libc::c_int = 104 as libc::c_int * nf * nc * n_comp;
    let leni: libc::c_int = 104 as libc::c_int * nf * i_ctr * n_comp;
    let len0: libc::c_int = 104 as libc::c_int * nf * n_comp;
    let len: libc::c_int = leng + lenj + leni + len0;
    let mut gridsT: *mut libc::c_double = 0 as *mut libc::c_double;
    gridsT = ((cache as uintptr_t).wrapping_add(63 as libc::c_int as libc::c_ulong)
        & (64 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut libc::c_double;
    cache = gridsT.offset((len + 104 as libc::c_int * 3 as libc::c_int) as isize);
    let mut g: *mut libc::c_double = gridsT
        .offset((104 as libc::c_int * 3 as libc::c_int) as isize);
    let mut g1: *mut libc::c_double = g.offset(leng as isize);
    let mut gout: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut gctri: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut gctrj: *mut libc::c_double = 0 as *mut libc::c_double;
    if n_comp == 1 as libc::c_int {
        gctrj = gctr;
    } else {
        gctrj = g1;
        g1 = g1.offset(lenj as isize);
    }
    if j_ctr == 1 as libc::c_int {
        gctri = gctrj;
        iempty = jempty;
    } else {
        gctri = g1;
        g1 = g1.offset(leni as isize);
    }
    if i_ctr == 1 as libc::c_int {
        gout = gctri;
        gempty = iempty;
    } else {
        gout = g1;
    }
    grids_offset = 0 as libc::c_int;
    while grids_offset < ngrids {
        (*envs).c2rust_unnamed.grids_offset = grids_offset;
        bgrids = if ngrids - grids_offset < 104 as libc::c_int {
            ngrids - grids_offset
        } else {
            104 as libc::c_int
        };
        i = 0 as libc::c_int;
        while i < bgrids {
            *gridsT
                .offset(
                    (i + 104 as libc::c_int * 0 as libc::c_int) as isize,
                ) = *grids
                .offset(
                    ((grids_offset + i) * 3 as libc::c_int + 0 as libc::c_int) as isize,
                );
            *gridsT
                .offset(
                    (i + 104 as libc::c_int * 1 as libc::c_int) as isize,
                ) = *grids
                .offset(
                    ((grids_offset + i) * 3 as libc::c_int + 1 as libc::c_int) as isize,
                );
            *gridsT
                .offset(
                    (i + 104 as libc::c_int * 2 as libc::c_int) as isize,
                ) = *grids
                .offset(
                    ((grids_offset + i) * 3 as libc::c_int + 2 as libc::c_int) as isize,
                );
            i += 1;
            i;
        }
        empty[0 as libc::c_int as usize] = 1 as libc::c_int;
        empty[1 as libc::c_int as usize] = 1 as libc::c_int;
        empty[2 as libc::c_int as usize] = 1 as libc::c_int;
        if n_comp == 1 as libc::c_int {
            gctrj = gctr.offset((grids_offset * nf * nc) as isize);
            if j_ctr == 1 as libc::c_int {
                gctri = gctrj;
            }
            if i_ctr == 1 as libc::c_int {
                gout = gctri;
            }
        }
        pdata_ij = pdata_base;
        jp = 0 as libc::c_int;
        while jp < j_prim {
            (*envs).aj[0 as libc::c_int as usize] = *aj.offset(jp as isize);
            if j_ctr == 1 as libc::c_int {
                fac1j = (*envs).common_factor * *cj.offset(jp as isize);
            } else {
                fac1j = (*envs).common_factor;
                *iempty = 1 as libc::c_int;
            }
            ip = 0 as libc::c_int;
            while ip < i_prim {
                if !((*pdata_ij).cceij > expcutoff) {
                    (*envs).ai[0 as libc::c_int as usize] = *ai.offset(ip as isize);
                    expij = (*pdata_ij).eij;
                    cutoff = expcutoff - (*pdata_ij).cceij;
                    rij = ((*pdata_ij).rij).as_mut_ptr();
                    (*envs)
                        .rij[0 as libc::c_int
                        as usize] = *rij.offset(0 as libc::c_int as isize);
                    (*envs)
                        .rij[1 as libc::c_int
                        as usize] = *rij.offset(1 as libc::c_int as isize);
                    (*envs)
                        .rij[2 as libc::c_int
                        as usize] = *rij.offset(2 as libc::c_int as isize);
                    if i_ctr == 1 as libc::c_int {
                        fac1i = fac1j * *ci.offset(ip as isize) * expij;
                    } else {
                        fac1i = fac1j * expij;
                    }
                    (*envs).fac[0 as libc::c_int as usize] = fac1i;
                    CINTg0_1e_grids(g, cutoff, envs, cache, gridsT);
                    ::core::mem::transmute::<
                        _,
                        fn(_, _, _, _, _),
                    >(
                        (Some(((*envs).f_gout).expect("non-null function pointer")))
                            .expect("non-null function pointer"),
                    )(gout, g, idx, envs, *gempty);
                    if i_ctr > 1 as libc::c_int {
                        if *iempty != 0 {
                            CINTprim_to_ctr_0(
                                gctri,
                                gout,
                                ci.offset(ip as isize),
                                (bgrids * nf * n_comp) as size_t,
                                i_prim,
                                i_ctr,
                                *non0ctri.offset(ip as isize),
                                non0idxi.offset((ip * i_ctr) as isize),
                            );
                        } else {
                            CINTprim_to_ctr_1(
                                gctri,
                                gout,
                                ci.offset(ip as isize),
                                (bgrids * nf * n_comp) as size_t,
                                i_prim,
                                i_ctr,
                                *non0ctri.offset(ip as isize),
                                non0idxi.offset((ip * i_ctr) as isize),
                            );
                        }
                    }
                    *iempty = 0 as libc::c_int;
                }
                ip += 1;
                ip;
                pdata_ij = pdata_ij.offset(1);
                pdata_ij;
            }
            if *iempty == 0 {
                if j_ctr > 1 as libc::c_int {
                    if *jempty != 0 {
                        CINTprim_to_ctr_0(
                            gctrj,
                            gctri,
                            cj.offset(jp as isize),
                            (bgrids * nf * i_ctr * n_comp) as size_t,
                            j_prim,
                            j_ctr,
                            *non0ctrj.offset(jp as isize),
                            non0idxj.offset((jp * j_ctr) as isize),
                        );
                    } else {
                        CINTprim_to_ctr_1(
                            gctrj,
                            gctri,
                            cj.offset(jp as isize),
                            (bgrids * nf * i_ctr * n_comp) as size_t,
                            j_prim,
                            j_ctr,
                            *non0ctrj.offset(jp as isize),
                            non0idxj.offset((jp * j_ctr) as isize),
                        );
                    }
                }
                *jempty = 0 as libc::c_int;
            }
            jp += 1;
            jp;
        }
        if n_comp > 1 as libc::c_int && *jempty == 0 {
            _transpose_comps(
                gctr.offset((grids_offset * nf * nc) as isize),
                gctrj,
                bgrids,
                nf * nc,
                ngrids,
                n_comp,
            );
        }
        all_empty &= *jempty;
        grids_offset += 104 as libc::c_int;
    }
    return (all_empty == 0) as libc::c_int;
}
unsafe extern "C" fn _transpose_comps(
    mut gctr: *mut libc::c_double,
    mut gctrj: *mut libc::c_double,
    mut bgrids: libc::c_int,
    mut dij: libc::c_int,
    mut ngrids: libc::c_int,
    mut n_comp: libc::c_int,
) {
    let mut n: libc::c_int = 0;
    let mut ic: libc::c_int = 0;
    let mut ig: libc::c_int = 0;
    let mut pgctr: *mut libc::c_double = 0 as *mut libc::c_double;
    let mut pgctrj: *mut libc::c_double = 0 as *mut libc::c_double;
    ic = 0 as libc::c_int;
    while ic < n_comp {
        pgctr = gctr.offset((ic * dij * ngrids) as isize);
        n = 0 as libc::c_int;
        while n < dij {
            pgctrj = gctrj.offset(((n * n_comp + ic) * bgrids) as isize);
            ig = 0 as libc::c_int;
            while ig < bgrids {
                *pgctr.offset((ig + n * bgrids) as isize) = *pgctrj.offset(ig as isize);
                ig += 1;
                ig;
            }
            n += 1;
            n;
        }
        ic += 1;
        ic;
    }
}
#[no_mangle]
pub unsafe extern "C" fn int1e_grids_cache_size(mut envs: *mut CINTEnvVars) -> size_t {
    let mut bas: *mut libc::c_int = (*envs).bas;
    let mut shls: *mut libc::c_int = (*envs).shls;
    let mut x_ctr: *mut libc::c_int = ((*envs).x_ctr).as_mut_ptr();
    let mut ngrids: libc::c_int = (*envs).c2rust_unnamed_0.ngrids;
    let mut nroots: libc::c_int = (*envs).nrys_roots;
    let mut nf: libc::c_int = (*envs).nf;
    let mut nc: libc::c_int = ngrids * nf * *x_ctr.offset(0 as libc::c_int as isize)
        * *x_ctr.offset(1 as libc::c_int as isize);
    let mut n_comp: libc::c_int = (*envs).ncomp_e1 * (*envs).ncomp_tensor;
    let mut i_prim: libc::c_int = *bas
        .offset(
            (8 as libc::c_int * *shls.offset(0 as libc::c_int as isize)
                + 2 as libc::c_int) as isize,
        );
    let mut j_prim: libc::c_int = *bas
        .offset(
            (8 as libc::c_int * *shls.offset(1 as libc::c_int as isize)
                + 2 as libc::c_int) as isize,
        );
    let mut pdata_size: libc::c_int = i_prim * j_prim * 5 as libc::c_int
        + i_prim * *x_ctr.offset(0 as libc::c_int as isize)
        + j_prim * *x_ctr.offset(1 as libc::c_int as isize)
        + (i_prim + j_prim) * 2 as libc::c_int + (*envs).nf * 3 as libc::c_int;
    let mut leng: size_t = ((*envs).g_size * 3 as libc::c_int
        * (((1 as libc::c_int) << (*envs).gbits) + 1 as libc::c_int)) as size_t;
    let mut len0: size_t = (104 as libc::c_int * nf * n_comp) as size_t;
    let mut leni: size_t = len0
        .wrapping_mul(*x_ctr.offset(0 as libc::c_int as isize) as libc::c_ulong);
    let mut lenj: size_t = leni
        .wrapping_mul(*x_ctr.offset(1 as libc::c_int as isize) as libc::c_ulong);
    let mut cache_size: size_t = if ((nc * n_comp) as libc::c_ulong)
        .wrapping_add(leng)
        .wrapping_add(len0)
        .wrapping_add(leni)
        .wrapping_add(lenj)
        .wrapping_add(pdata_size as libc::c_ulong)
        .wrapping_add(
            (104 as libc::c_int
                * (if n_comp > nroots + 10 as libc::c_int {
                    n_comp
                } else {
                    nroots + 10 as libc::c_int
                })) as libc::c_ulong,
        )
        > (nc * n_comp + 104 as libc::c_int * nf * 8 as libc::c_int * 2 as libc::c_int)
            as libc::c_ulong
    {
        ((nc * n_comp) as libc::c_ulong)
            .wrapping_add(leng)
            .wrapping_add(len0)
            .wrapping_add(leni)
            .wrapping_add(lenj)
            .wrapping_add(pdata_size as libc::c_ulong)
            .wrapping_add(
                (104 as libc::c_int
                    * (if n_comp > nroots + 10 as libc::c_int {
                        n_comp
                    } else {
                        nroots + 10 as libc::c_int
                    })) as libc::c_ulong,
            )
    } else {
        (nc * n_comp + 104 as libc::c_int * nf * 8 as libc::c_int * 2 as libc::c_int)
            as libc::c_ulong
    };
    return cache_size.wrapping_add(32 as libc::c_int as libc::c_ulong);
}
#[no_mangle]
pub unsafe extern "C" fn CINT1e_grids_drv(
    mut out: *mut libc::c_double,
    mut dims: *mut libc::c_int,
    mut envs: *mut CINTEnvVars,
    mut cache: *mut libc::c_double,
    mut f_c2s: Option::<unsafe extern "C" fn() -> ()>,
) -> libc::c_int {
    if out.is_null() {
        return int1e_grids_cache_size(envs) as libc::c_int;
    }
    let mut x_ctr: *mut libc::c_int = ((*envs).x_ctr).as_mut_ptr();
    let mut ngrids_nf: libc::c_int = (*envs).c2rust_unnamed_0.ngrids * (*envs).nf;
    let mut nc: libc::c_int = ngrids_nf * *x_ctr.offset(0 as libc::c_int as isize)
        * *x_ctr.offset(1 as libc::c_int as isize);
    let mut n_comp: libc::c_int = (*envs).ncomp_e1 * (*envs).ncomp_tensor;
    let mut stack: *mut libc::c_double = 0 as *mut libc::c_double;
    if cache.is_null() {
        let mut cache_size: size_t = int1e_grids_cache_size(envs);
        stack = malloc(
            (::core::mem::size_of::<libc::c_double>() as libc::c_ulong)
                .wrapping_mul(cache_size),
        ) as *mut libc::c_double;
        cache = stack;
    }
    let mut gctr: *mut libc::c_double = 0 as *mut libc::c_double;
    gctr = ((cache as uintptr_t).wrapping_add(63 as libc::c_int as libc::c_ulong)
        & (64 as libc::c_int as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut libc::c_double;
    cache = gctr.offset((nc * n_comp) as isize);
    let mut has_value: libc::c_int = CINT1e_grids_loop(gctr, envs, cache);
    let mut counts: [libc::c_int; 4] = [0; 4];
    if dims.is_null() {
        dims = counts.as_mut_ptr();
    }
    if f_c2s
        == ::core::mem::transmute::<
            Option::<
                unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    *mut libc::c_double,
                ) -> (),
            >,
            Option::<unsafe extern "C" fn() -> ()>,
        >(
            Some(
                c2s_sph_1e_grids
                    as unsafe extern "C" fn(
                        *mut libc::c_double,
                        *mut libc::c_double,
                        *mut libc::c_int,
                        *mut CINTEnvVars,
                        *mut libc::c_double,
                    ) -> (),
            ),
        )
    {
        counts[0 as libc::c_int
            as usize] = ((*envs).i_l * 2 as libc::c_int + 1 as libc::c_int)
            * *x_ctr.offset(0 as libc::c_int as isize);
        counts[1 as libc::c_int
            as usize] = ((*envs).j_l * 2 as libc::c_int + 1 as libc::c_int)
            * *x_ctr.offset(1 as libc::c_int as isize);
        counts[2 as libc::c_int as usize] = (*envs).c2rust_unnamed_0.ngrids;
        counts[3 as libc::c_int as usize] = 1 as libc::c_int;
    } else if f_c2s
        == ::core::mem::transmute::<
            Option::<
                unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    *mut libc::c_double,
                ) -> (),
            >,
            Option::<unsafe extern "C" fn() -> ()>,
        >(
            Some(
                c2s_cart_1e_grids
                    as unsafe extern "C" fn(
                        *mut libc::c_double,
                        *mut libc::c_double,
                        *mut libc::c_int,
                        *mut CINTEnvVars,
                        *mut libc::c_double,
                    ) -> (),
            ),
        )
    {
        counts[0 as libc::c_int
            as usize] = (*envs).nfi * *x_ctr.offset(0 as libc::c_int as isize);
        counts[1 as libc::c_int
            as usize] = (*envs).nfj * *x_ctr.offset(1 as libc::c_int as isize);
        counts[2 as libc::c_int as usize] = (*envs).c2rust_unnamed_0.ngrids;
        counts[3 as libc::c_int as usize] = 1 as libc::c_int;
    }
    let mut nout: libc::c_int = *dims.offset(0 as libc::c_int as isize)
        * *dims.offset(1 as libc::c_int as isize)
        * *dims.offset(2 as libc::c_int as isize);
    let mut n: libc::c_int = 0;
    if has_value != 0 {
        n = 0 as libc::c_int;
        while n < n_comp {
            ::core::mem::transmute::<
                _,
                fn(_, _, _, _, _),
            >(
                (Some(f_c2s.expect("non-null function pointer")))
                    .expect("non-null function pointer"),
            )(
                out.offset((nout * n) as isize),
                gctr.offset((nc * n) as isize),
                dims,
                envs,
                cache,
            );
            n += 1;
            n;
        }
    } else {
        n = 0 as libc::c_int;
        while n < n_comp {
            c2s_grids_dset0(out.offset((nout * n) as isize), dims, counts.as_mut_ptr());
            n += 1;
            n;
        }
    }
    if !stack.is_null() {
        free(stack as *mut libc::c_void);
    }
    return has_value;
}
#[no_mangle]
pub unsafe extern "C" fn int1e_grids_sph(
    mut out: *mut libc::c_double,
    mut dims: *mut libc::c_int,
    mut shls: *mut libc::c_int,
    mut atm: *mut libc::c_int,
    mut natm: libc::c_int,
    mut bas: *mut libc::c_int,
    mut nbas: libc::c_int,
    mut env: *mut libc::c_double,
    mut opt: *mut CINTOpt,
    mut cache: *mut libc::c_double,
) -> libc::c_int {
    let mut ng: [libc::c_int; 8] = [
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        1 as libc::c_int,
        0 as libc::c_int,
        1 as libc::c_int,
    ];
    let mut envs: CINTEnvVars = CINTEnvVars::new();
    CINTinit_int1e_grids_EnvVars(
        &mut envs,
        ng.as_mut_ptr(),
        shls,
        atm,
        natm,
        bas,
        nbas,
        env,
    );
    envs
        .f_gout = ::core::mem::transmute::<
        Option::<
            unsafe extern "C" fn(
                *mut libc::c_double,
                *mut libc::c_double,
                *mut libc::c_int,
                *mut CINTEnvVars,
                libc::c_int,
            ) -> (),
        >,
        Option::<unsafe extern "C" fn() -> ()>,
    >(
        Some(
            CINTgout1e_grids
                as unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    libc::c_int,
                ) -> (),
        ),
    );
    return CINT1e_grids_drv(
        out,
        dims,
        &mut envs,
        cache,
        ::core::mem::transmute::<
            Option::<
                unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    *mut libc::c_double,
                ) -> (),
            >,
            Option::<unsafe extern "C" fn() -> ()>,
        >(
            Some(
                c2s_sph_1e_grids
                    as unsafe extern "C" fn(
                        *mut libc::c_double,
                        *mut libc::c_double,
                        *mut libc::c_int,
                        *mut CINTEnvVars,
                        *mut libc::c_double,
                    ) -> (),
            ),
        ),
    );
}
#[no_mangle]
pub unsafe extern "C" fn int1e_grids_optimizer(
    mut opt: *mut *mut CINTOpt,
    mut atm: *mut libc::c_int,
    mut natm: libc::c_int,
    mut bas: *mut libc::c_int,
    mut nbas: libc::c_int,
    mut env: *mut libc::c_double,
) {
    *opt = 0 as *mut CINTOpt;
}
#[no_mangle]
pub unsafe extern "C" fn int1e_grids_cart(
    mut out: *mut libc::c_double,
    mut dims: *mut libc::c_int,
    mut shls: *mut libc::c_int,
    mut atm: *mut libc::c_int,
    mut natm: libc::c_int,
    mut bas: *mut libc::c_int,
    mut nbas: libc::c_int,
    mut env: *mut libc::c_double,
    mut opt: *mut CINTOpt,
    mut cache: *mut libc::c_double,
) -> libc::c_int {
    let mut ng: [libc::c_int; 8] = [
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        0 as libc::c_int,
        1 as libc::c_int,
        0 as libc::c_int,
        1 as libc::c_int,
    ];
    let mut envs: CINTEnvVars = CINTEnvVars::new();
    CINTinit_int1e_grids_EnvVars(
        &mut envs,
        ng.as_mut_ptr(),
        shls,
        atm,
        natm,
        bas,
        nbas,
        env,
    );
    envs
        .f_gout = ::core::mem::transmute::<
        Option::<
            unsafe extern "C" fn(
                *mut libc::c_double,
                *mut libc::c_double,
                *mut libc::c_int,
                *mut CINTEnvVars,
                libc::c_int,
            ) -> (),
        >,
        Option::<unsafe extern "C" fn() -> ()>,
    >(
        Some(
            CINTgout1e_grids
                as unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    libc::c_int,
                ) -> (),
        ),
    );
    return CINT1e_grids_drv(
        out,
        dims,
        &mut envs,
        cache,
        ::core::mem::transmute::<
            Option::<
                unsafe extern "C" fn(
                    *mut libc::c_double,
                    *mut libc::c_double,
                    *mut libc::c_int,
                    *mut CINTEnvVars,
                    *mut libc::c_double,
                ) -> (),
            >,
            Option::<unsafe extern "C" fn() -> ()>,
        >(
            Some(
                c2s_cart_1e_grids
                    as unsafe extern "C" fn(
                        *mut libc::c_double,
                        *mut libc::c_double,
                        *mut libc::c_int,
                        *mut CINTEnvVars,
                        *mut libc::c_double,
                    ) -> (),
            ),
        ),
    );
}
