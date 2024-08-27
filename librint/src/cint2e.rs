#![allow(dead_code, mutable_transmutes, non_camel_case_types, non_snake_case, non_upper_case_globals, unused_assignments, unused_mut)]

use crate::g1e::CINTprim_to_ctr_0;
use crate::g1e::CINTprim_to_ctr_1;
use crate::g2e::CINTg2e_index_xyz;
use crate::g2e::CINTinit_int2e_EnvVars;
use crate::optimizer::CINTOpt_log_max_pgto_coeff;
use crate::optimizer::CINTOpt_non0coeff_byshell;
use crate::optimizer::CINTset_pairdata;
// use crate::optimizer::CINTall_2e_optimizer;
use crate::fblas::CINTdplus_transpose;
use crate::fblas::CINTdmat_transpose;
// use crate::cart2sph::c2s_sph_2e1;
// use crate::cart2sph::c2s_cart_2e1;
// use crate::cart2sph::c2s_dset0;

use crate::cart2sph::c2s_dset0_cpy;

use crate::cart2sph::c2s_sph_2e1_cpy;
use crate::cart2sph::c2s_cart_2e1_cpy;

use crate::cint::PairData;
// use crate::cint::CINTOpt;
use crate::cint::CINTEnvVars;

use crate::cint::F_FC2S;

extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
}

pub type size_t = libc::c_ulong;
pub type uintptr_t = libc::c_ulong;

#[no_mangle]
pub unsafe extern "C" fn CINT2e_loop_nopt_cpy(
    gctr: &mut [f64],
    envs: &CINTEnvVars,
    cache: &mut [f64],
    empty: &mut i32,
) -> i32 {
    let shls: [i32; 4] = envs.shls;
    let bas: &[i32] = &envs.bas;
    let env: &[f64] = &envs.env;
    let i_sh: usize = shls[0] as usize;
    let j_sh: usize = shls[1] as usize;
    let k_sh: usize = shls[2] as usize;
    let l_sh: usize = shls[3] as usize;
    let i_ctr: usize = envs.x_ctr[0] as usize;
    let j_ctr: usize = envs.x_ctr[1] as usize;
    let k_ctr: usize = envs.x_ctr[2] as usize;
    let l_ctr: usize = envs.x_ctr[3] as usize;
    let i_prim: i32 = bas[8 * i_sh + 2];
    let j_prim: i32 = bas[8 * j_sh + 2];
    let k_prim: i32 = bas[8 * k_sh + 2];
    let l_prim: i32 = bas[8 * l_sh + 2];
    let rk: [f64; 3] = envs.rk;
    let rl: [f64; 3] = envs.c2rust_unnamed_1.rl;
    
    let mut ai: *mut f64 = env.offset(*bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,);
    let mut aj: *mut f64 = env.offset(*bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,);
    let mut ak: *mut f64 = env.offset(*bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,);
    let mut al: *mut f64 = env.offset(*bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,);

    let mut ci: *mut f64 = env.offset(*bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,);
    let mut cj: *mut f64 = env.offset(*bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,);
    let mut ck: *mut f64 = env.offset(*bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,);
    let mut cl: *mut f64 = env.offset(*bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,);
    
    let expcutoff: f64 = envs.expcutoff;
    let mut rr_ij: f64 = envs.rirj[0] * envs.rirj[0] + envs.rirj[1] * envs.rirj[1] + envs.rirj[2] * envs.rirj[2];
    let mut rr_kl: f64 = envs.rkrl[0] * envs.rkrl[0] + envs.rkrl[1] * envs.rkrl[1] + envs.rkrl[2] * envs.rkrl[2];

    // let mut log_maxci: *mut f64 = 0 as *mut f64;
    // let mut log_maxcj: *mut f64 = 0 as *mut f64;
    // let mut log_maxck: *mut f64 = 0 as *mut f64;
    // let mut log_maxcl: *mut f64 = 0 as *mut f64;
    // let mut pdata_base: *mut PairData = 0 as *mut PairData;
    // let mut pdata_ij: *mut PairData = 0 as *mut PairData;
    // log_maxci = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
    //     & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
    //     as *mut f64;
    // cache = log_maxci.offset((i_prim + j_prim + k_prim + l_prim) as isize);
    // pdata_base = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
    //     & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
    //     as *mut PairData;
    // cache = pdata_base.offset((i_prim * j_prim) as isize) as *mut f64;
    // log_maxcj = log_maxci.offset(i_prim as isize);
    // log_maxck = log_maxcj.offset(j_prim as isize);
    // log_maxcl = log_maxck.offset(k_prim as isize);

    let mut log_maxci = vec![0.0; i_prim as usize].into_boxed_slice();
    let mut log_maxcj = vec![0.0; j_prim as usize].into_boxed_slice();
    let mut log_maxck = vec![0.0; k_prim as usize].into_boxed_slice();
    let mut log_maxcl = vec![0.0; l_prim as usize].into_boxed_slice();

    let pdata_base = vec![0.0; (i_prim * j_prim) as usize].into_boxed_slice();


    CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxck, ck, k_prim, k_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcl, cl, l_prim, l_ctr);

    if CINTset_pairdata(
        pdata_base,
        ai,
        aj,
        (*envs).ri.as_mut_ptr(),
        (*envs).rj.as_mut_ptr(),
        log_maxci,
        log_maxcj,
        (*envs).li_ceil,
        (*envs).lj_ceil,
        i_prim,
        j_prim,
        rr_ij,
        expcutoff,
        env,
    ) != 0
    {
        return 0 as i32;
    }

    let n_comp: i32 = envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor;
    let nf: size_t = envs.nf as size_t;
    let fac1i: f64 = 0.;
    let fac1j: f64 = 0.;
    let fac1k: f64 = 0.;
    let fac1l: f64 = 0.;
    let mut ip: i32 = 0;
    let mut jp: i32 = 0;
    let mut kp: i32 = 0;
    let mut lp: i32 = 0;

    let mut empty: [i32; 5] = [1, 1, 1, 1, 1];

    let mut iempty_idx = 0;
    let mut jempty_idx = 1;
    let mut kempty_idx = 2;
    let mut lempty_idx = 3;
    let mut gempty_idx = 4;

    // let mut iempty: *mut i32 = _empty
    //     .as_mut_ptr()
    //     .offset(0 as i32 as isize);
    // let mut jempty: *mut i32 = _empty
    //     .as_mut_ptr()
    //     .offset(1 as i32 as isize);
    // let mut kempty: *mut i32 = _empty
    //     .as_mut_ptr()
    //     .offset(2 as i32 as isize);
    // let mut lempty: *mut i32 = _empty
    //     .as_mut_ptr()
    //     .offset(3 as i32 as isize);
    // let mut gempty: *mut i32 = _empty
    //     .as_mut_ptr()
    //     .offset(4 as i32 as isize);

    let lkl: i32 = envs.lk_ceil + envs.ll_ceil;
    let akl: f64 = 0.;
    let ekl: f64 = 0.;
    let expijkl: f64 = 0.;
    let ccekl: f64 = 0.;
    let log_rr_kl: f64 = 0.;
    let eijcutoff: f64 = 0.;
    let cutoff: f64 = 0.;
    let rkl: [f64; 3] = [0.; 3];
    let rij: *mut f64 = 0 as *mut f64;

    akl = ak[k_prim - 1] + al[l_prim - 1];

    log_rr_kl = 1.7f64 - 1.5f64 * akl.ln(); //log(akl);
    let mut omega: f64 = env[8];

    if omega < 0.0 {
        if envs.rys_order > 1 {
            let r_guess: f64 = 8.0f64;
            let omega2: f64 = omega * omega;
            let lij: i32 = envs.li_ceil + envs.lj_ceil;
            if lij > 0 {
                let aij: f64 = ai[i_prim - 1] + aj[j_prim - 1];
                let dist_ij: f64 = rr_ij.sqrt();
                let theta: f64 = omega2 / (omega2 + aij);
                expcutoff += lij as f64 * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
            }
            if lkl > 0 as i32 {
                let theta_0: f64 = omega2 / (omega2 + akl);
                log_rr_kl += lkl as f64 * (rr_kl.sqrt() + theta_0 * r_guess + 1.0f64).ln(); //log(sqrt(rr_kl) + theta_0 * r_guess + 1.0f64);
            }
        }
    } else if lkl > 0 {
        log_rr_kl += lkl as f64 * (rr_kl.sqrt() + 1.0f64).ln(); //log(sqrt(rr_kl) + 1.0f64);
    }

    // let mut idx: *mut i32 = 0 as *mut i32;
    // idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
    //     & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
    //     as *mut i32;
    // cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
    //     as *mut f64;

    let mut idx = vec![0; nf as usize* 3];
    
    CINTg2e_index_xyz(idx.as_mut_ptr(), envs);
    let mut non0ctri: *mut i32 = 0 as *mut i32;
    let mut non0ctrj: *mut i32 = 0 as *mut i32;
    let mut non0ctrk: *mut i32 = 0 as *mut i32;
    let mut non0ctrl: *mut i32 = 0 as *mut i32;
    let mut non0idxi: *mut i32 = 0 as *mut i32;
    let mut non0idxj: *mut i32 = 0 as *mut i32;
    let mut non0idxk: *mut i32 = 0 as *mut i32;
    let mut non0idxl: *mut i32 = 0 as *mut i32;
    non0ctri = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut i32;
    cache = non0ctri
        .offset(
            (i_prim + j_prim + k_prim + l_prim + i_prim * i_ctr + j_prim * j_ctr
                + k_prim * k_ctr + l_prim * l_ctr) as isize,
        ) as *mut f64;
    non0ctrj = non0ctri.offset(i_prim as isize);
    non0ctrk = non0ctrj.offset(j_prim as isize);
    non0ctrl = non0ctrk.offset(k_prim as isize);
    non0idxi = non0ctrl.offset(l_prim as isize);
    non0idxj = non0idxi.offset((i_prim * i_ctr) as isize);
    non0idxk = non0idxj.offset((j_prim * j_ctr) as isize);
    non0idxl = non0idxk.offset((k_prim * k_ctr) as isize);
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, cl, l_prim, l_ctr);
    let mut nc: i32 = i_ctr * j_ctr * k_ctr * l_ctr;
    let mut leng: size_t = ((*envs).g_size * 3 as i32
        * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
    let mut lenl: size_t = nf
        .wrapping_mul(nc as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut lenk: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(j_ctr as libc::c_ulong)
        .wrapping_mul(k_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut lenj: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(j_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut leni: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
    let mut len: size_t = leng
        .wrapping_add(lenl)
        .wrapping_add(lenk)
        .wrapping_add(lenj)
        .wrapping_add(leni)
        .wrapping_add(len0);
    let mut g: *mut f64 = 0 as *mut f64;
    g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut f64;
    cache = g.offset(len as isize);
    let mut g1: *mut f64 = g.offset(leng as isize);
    let mut gout: *mut f64 = 0 as *mut f64;
    let mut gctri: *mut f64 = 0 as *mut f64;
    let mut gctrj: *mut f64 = 0 as *mut f64;
    let mut gctrk: *mut f64 = 0 as *mut f64;
    let mut gctrl: *mut f64 = 0 as *mut f64;
    if n_comp == 1 as i32 {
        gctrl = gctr;
        lempty = empty;
    } else {
        gctrl = g1;
        g1 = g1.offset(lenl as isize);
    }
    if l_ctr == 1 as i32 {
        gctrk = gctrl;
        kempty = lempty;
    } else {
        gctrk = g1;
        g1 = g1.offset(lenk as isize);
    }
    if k_ctr == 1 as i32 {
        gctrj = gctrk;
        jempty = kempty;
    } else {
        gctrj = g1;
        g1 = g1.offset(lenj as isize);
    }
    if j_ctr == 1 as i32 {
        gctri = gctrj;
        iempty = jempty;
    } else {
        gctri = g1;
        g1 = g1.offset(leni as isize);
    }
    if i_ctr == 1 as i32 {
        gout = gctri;
        gempty = iempty;
    } else {
        gout = g1;
        g1 = g1.offset(leng as isize);
    }
    lp = 0 as i32;
    while lp < l_prim {
        (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
        if l_ctr == 1 as i32 {
            fac1l = (*envs).common_factor * *cl.offset(lp as isize);
        } else {
            fac1l = (*envs).common_factor;
            *kempty = 1 as i32;
        }
        kp = 0 as i32;
        while kp < k_prim {
            akl = *ak.offset(kp as isize) + *al.offset(lp as isize);
            ekl = rr_kl * *ak.offset(kp as isize) * *al.offset(lp as isize) / akl;
            ccekl = ekl - log_rr_kl - *log_maxck.offset(kp as isize)
                - *log_maxcl.offset(lp as isize);
            if !(ccekl > expcutoff) {
                (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
                rkl[0 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(0 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(0 as i32 as isize))
                    / akl;
                rkl[1 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(1 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(1 as i32 as isize))
                    / akl;
                rkl[2 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(2 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(2 as i32 as isize))
                    / akl;
                eijcutoff = expcutoff - ccekl;
                ekl = (-ekl).exp(); // exp(-ekl);
                if k_ctr == 1 as i32 {
                    fac1k = fac1l * *ck.offset(kp as isize);
                } else {
                    fac1k = fac1l;
                    *jempty = 1 as i32;
                }
                pdata_ij = pdata_base;
                jp = 0 as i32;
                while jp < j_prim {
                    (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
                    if j_ctr == 1 as i32 {
                        fac1j = fac1k * *cj.offset(jp as isize);
                    } else {
                        fac1j = fac1k;
                        *iempty = 1 as i32;
                    }
                    ip = 0 as i32;
                    while ip < i_prim {
                        if !((*pdata_ij).cceij > eijcutoff) {
                            (*envs)
                                .ai[0 as i32 as usize] = *ai.offset(ip as isize);
                            rij = ((*pdata_ij).rij).as_mut_ptr();
                            cutoff = eijcutoff - (*pdata_ij).cceij;
                            expijkl = (*pdata_ij).eij * ekl;
                            if i_ctr == 1 as i32 {
                                fac1i = fac1j * *ci.offset(ip as isize) * expijkl;
                            } else {
                                fac1i = fac1j * expijkl;
                            }
                            (*envs).fac[0 as i32 as usize] = fac1i;
                            if ::core::mem::transmute::<
                                _,
                                fn(_, _, _, _, _) -> i32,
                            >(
                                (Some(
                                    ((*envs).f_g0_2e).expect("non-null function pointer"),
                                ))
                                    .expect("non-null function pointer"),
                            )(g, rij, rkl.as_mut_ptr(), cutoff, envs) != 0
                            {
                                ::core::mem::transmute::<
                                    _,
                                    fn(_, _, _, _, _),
                                >(
                                    (Some(((*envs).f_gout).expect("non-null function pointer")))
                                        .expect("non-null function pointer"),
                                )(gout, g, idx, envs, *gempty);
                                if i_ctr > 1 as i32 {
                                    if *iempty != 0 {
                                        CINTprim_to_ctr_0(
                                            gctri,
                                            gout,
                                            ci.offset(ip as isize),
                                            len0,
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
                                            len0,
                                            i_prim,
                                            i_ctr,
                                            *non0ctri.offset(ip as isize),
                                            non0idxi.offset((ip * i_ctr) as isize),
                                        );
                                    }
                                }
                                *iempty = 0 as i32;
                            }
                        }
                        ip += 1;
                        pdata_ij = pdata_ij.offset(1);
                    }
                    if *iempty == 0 {
                        if j_ctr > 1 as i32 {
                            if *jempty != 0 {
                                CINTprim_to_ctr_0(
                                    gctrj,
                                    gctri,
                                    cj.offset(jp as isize),
                                    leni,
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
                                    leni,
                                    j_prim,
                                    j_ctr,
                                    *non0ctrj.offset(jp as isize),
                                    non0idxj.offset((jp * j_ctr) as isize),
                                );
                            }
                        }
                        *jempty = 0 as i32;
                    }
                    jp += 1;
                }
                if *jempty == 0 {
                    if k_ctr > 1 as i32 {
                        if *kempty != 0 {
                            CINTprim_to_ctr_0(
                                gctrk,
                                gctrj,
                                ck.offset(kp as isize),
                                lenj,
                                k_prim,
                                k_ctr,
                                *non0ctrk.offset(kp as isize),
                                non0idxk.offset((kp * k_ctr) as isize),
                            );
                        } else {
                            CINTprim_to_ctr_1(
                                gctrk,
                                gctrj,
                                ck.offset(kp as isize),
                                lenj,
                                k_prim,
                                k_ctr,
                                *non0ctrk.offset(kp as isize),
                                non0idxk.offset((kp * k_ctr) as isize),
                            );
                        }
                    }
                    *kempty = 0 as i32;
                }
            }
            kp += 1;
        }
        if *kempty == 0 {
            if l_ctr > 1 as i32 {
                if *lempty != 0 {
                    CINTprim_to_ctr_0(
                        gctrl,
                        gctrk,
                        cl.offset(lp as isize),
                        lenk,
                        l_prim,
                        l_ctr,
                        *non0ctrl.offset(lp as isize),
                        non0idxl.offset((lp * l_ctr) as isize),
                    );
                } else {
                    CINTprim_to_ctr_1(
                        gctrl,
                        gctrk,
                        cl.offset(lp as isize),
                        lenk,
                        l_prim,
                        l_ctr,
                        *non0ctrl.offset(lp as isize),
                        non0idxl.offset((lp * l_ctr) as isize),
                    );
                }
            }
            *lempty = 0 as i32;
        }
        lp += 1;
    }
    if n_comp > 1 as i32 && *lempty == 0 {
        if *empty != 0 {
            CINTdmat_transpose(
                gctr,
                gctrl,
                nf.wrapping_mul(nc as libc::c_ulong) as i32,
                n_comp,
            );
            *empty = 0 as i32;
        } else {
            CINTdplus_transpose(
                gctr,
                gctrl,
                nf.wrapping_mul(nc as libc::c_ulong) as i32,
                n_comp,
            );
        }
    }
    return (*empty == 0) as i32;
}

#[no_mangle]
pub unsafe extern "C" fn CINT2e_loop_nopt(
    mut gctr: *mut f64,
    mut envs: *mut CINTEnvVars,
    mut cache: *mut f64,
    mut empty: *mut i32,
) -> i32 {
    let mut shls: *mut i32 = (*envs).shls.as_mut_ptr();
    let mut bas: *mut i32 = (*envs).bas.as_mut_ptr();
    let mut env: *mut f64 = (*envs).env.as_mut_ptr();
    let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
    let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
    let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
    let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
    let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
    let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
    let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
    let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
    let mut i_prim: i32 = *bas
        .offset((8 as i32 * i_sh + 2 as i32) as isize);
    let mut j_prim: i32 = *bas
        .offset((8 as i32 * j_sh + 2 as i32) as isize);
    let mut k_prim: i32 = *bas
        .offset((8 as i32 * k_sh + 2 as i32) as isize);
    let mut l_prim: i32 = *bas
        .offset((8 as i32 * l_sh + 2 as i32) as isize);
    let mut rk: *mut f64 = (*envs).rk.as_mut_ptr();
    let mut rl: *mut f64 = (*envs).c2rust_unnamed_1.rl.as_mut_ptr();
    let mut ai: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
        );
    let mut aj: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
        );
    let mut ak: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
        );
    let mut al: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
        );
    let mut ci: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
        );
    let mut cj: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
        );
    let mut ck: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
        );
    let mut cl: *mut f64 = env
        .offset(
            *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
        );
    let mut expcutoff: f64 = (*envs).expcutoff;
    let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
        * (*envs).rirj[0 as i32 as usize]
        + (*envs).rirj[1 as i32 as usize]
            * (*envs).rirj[1 as i32 as usize]
        + (*envs).rirj[2 as i32 as usize]
            * (*envs).rirj[2 as i32 as usize];
    let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
        * (*envs).rkrl[0 as i32 as usize]
        + (*envs).rkrl[1 as i32 as usize]
            * (*envs).rkrl[1 as i32 as usize]
        + (*envs).rkrl[2 as i32 as usize]
            * (*envs).rkrl[2 as i32 as usize];
    let mut log_maxci: *mut f64 = 0 as *mut f64;
    let mut log_maxcj: *mut f64 = 0 as *mut f64;
    let mut log_maxck: *mut f64 = 0 as *mut f64;
    let mut log_maxcl: *mut f64 = 0 as *mut f64;
    let mut pdata_base: *mut PairData = 0 as *mut PairData;
    let mut pdata_ij: *mut PairData = 0 as *mut PairData;
    log_maxci = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut f64;
    cache = log_maxci.offset((i_prim + j_prim + k_prim + l_prim) as isize);
    pdata_base = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut PairData;
    cache = pdata_base.offset((i_prim * j_prim) as isize) as *mut f64;
    log_maxcj = log_maxci.offset(i_prim as isize);
    log_maxck = log_maxcj.offset(j_prim as isize);
    log_maxcl = log_maxck.offset(k_prim as isize);
    CINTOpt_log_max_pgto_coeff(log_maxci, ci, i_prim, i_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcj, cj, j_prim, j_ctr);
    if CINTset_pairdata(
        pdata_base,
        ai,
        aj,
        (*envs).ri.as_mut_ptr(),
        (*envs).rj.as_mut_ptr(),
        log_maxci,
        log_maxcj,
        (*envs).li_ceil,
        (*envs).lj_ceil,
        i_prim,
        j_prim,
        rr_ij,
        expcutoff,
        env,
    ) != 0
    {
        return 0 as i32;
    }
    CINTOpt_log_max_pgto_coeff(log_maxck, ck, k_prim, k_ctr);
    CINTOpt_log_max_pgto_coeff(log_maxcl, cl, l_prim, l_ctr);
    let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
        * (*envs).ncomp_tensor;
    let mut nf: size_t = (*envs).nf as size_t;
    let mut fac1i: f64 = 0.;
    let mut fac1j: f64 = 0.;
    let mut fac1k: f64 = 0.;
    let mut fac1l: f64 = 0.;
    let mut ip: i32 = 0;
    let mut jp: i32 = 0;
    let mut kp: i32 = 0;
    let mut lp: i32 = 0;
    let mut _empty: [i32; 5] = [
        1 as i32,
        1 as i32,
        1 as i32,
        1 as i32,
        1 as i32,
    ];
    let mut iempty: *mut i32 = _empty
        .as_mut_ptr()
        .offset(0 as i32 as isize);
    let mut jempty: *mut i32 = _empty
        .as_mut_ptr()
        .offset(1 as i32 as isize);
    let mut kempty: *mut i32 = _empty
        .as_mut_ptr()
        .offset(2 as i32 as isize);
    let mut lempty: *mut i32 = _empty
        .as_mut_ptr()
        .offset(3 as i32 as isize);
    let mut gempty: *mut i32 = _empty
        .as_mut_ptr()
        .offset(4 as i32 as isize);
    let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
    let mut akl: f64 = 0.;
    let mut ekl: f64 = 0.;
    let mut expijkl: f64 = 0.;
    let mut ccekl: f64 = 0.;
    let mut log_rr_kl: f64 = 0.;
    let mut eijcutoff: f64 = 0.;
    let mut cutoff: f64 = 0.;
    let mut rkl: [f64; 3] = [0.; 3];
    let mut rij: *mut f64 = 0 as *mut f64;
    akl = *ak.offset((k_prim - 1 as i32) as isize)
        + *al.offset((l_prim - 1 as i32) as isize);
    log_rr_kl = 1.7f64 - 1.5f64 * akl.ln(); //log(akl);
    let mut omega: f64 = *env.offset(8 as i32 as isize);
    if omega < 0 as i32 as f64 {
        if (*envs).rys_order > 1 as i32 {
            let mut r_guess: f64 = 8.0f64;
            let mut omega2: f64 = omega * omega;
            let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
            if lij > 0 as i32 {
                let mut aij: f64 = *ai
                    .offset((i_prim - 1 as i32) as isize)
                    + *aj.offset((j_prim - 1 as i32) as isize);
                let mut dist_ij: f64 = rr_ij.sqrt();
                let mut theta: f64 = omega2 / (omega2 + aij);
                expcutoff
                    += lij as f64
                        * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
            }
            if lkl > 0 as i32 {
                let mut theta_0: f64 = omega2 / (omega2 + akl);
                log_rr_kl
                    += lkl as f64
                        * (rr_kl.sqrt() + theta_0 * r_guess + 1.0f64).ln(); //log(sqrt(rr_kl) + theta_0 * r_guess + 1.0f64);
            }
        }
    } else if lkl > 0 as i32 {
        log_rr_kl += lkl as f64 * (rr_kl.sqrt() + 1.0f64).ln(); //log(sqrt(rr_kl) + 1.0f64);
    }
    let mut idx: *mut i32 = 0 as *mut i32;
    idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut i32;
    cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
        as *mut f64;
    CINTg2e_index_xyz(idx, envs);
    let mut non0ctri: *mut i32 = 0 as *mut i32;
    let mut non0ctrj: *mut i32 = 0 as *mut i32;
    let mut non0ctrk: *mut i32 = 0 as *mut i32;
    let mut non0ctrl: *mut i32 = 0 as *mut i32;
    let mut non0idxi: *mut i32 = 0 as *mut i32;
    let mut non0idxj: *mut i32 = 0 as *mut i32;
    let mut non0idxk: *mut i32 = 0 as *mut i32;
    let mut non0idxl: *mut i32 = 0 as *mut i32;
    non0ctri = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut i32;
    cache = non0ctri
        .offset(
            (i_prim + j_prim + k_prim + l_prim + i_prim * i_ctr + j_prim * j_ctr
                + k_prim * k_ctr + l_prim * l_ctr) as isize,
        ) as *mut f64;
    non0ctrj = non0ctri.offset(i_prim as isize);
    non0ctrk = non0ctrj.offset(j_prim as isize);
    non0ctrl = non0ctrk.offset(k_prim as isize);
    non0idxi = non0ctrl.offset(l_prim as isize);
    non0idxj = non0idxi.offset((i_prim * i_ctr) as isize);
    non0idxk = non0idxj.offset((j_prim * j_ctr) as isize);
    non0idxl = non0idxk.offset((k_prim * k_ctr) as isize);
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
    CINTOpt_non0coeff_byshell(non0idxl, non0ctrl, cl, l_prim, l_ctr);
    let mut nc: i32 = i_ctr * j_ctr * k_ctr * l_ctr;
    let mut leng: size_t = ((*envs).g_size * 3 as i32
        * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
    let mut lenl: size_t = nf
        .wrapping_mul(nc as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut lenk: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(j_ctr as libc::c_ulong)
        .wrapping_mul(k_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut lenj: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(j_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut leni: size_t = nf
        .wrapping_mul(i_ctr as libc::c_ulong)
        .wrapping_mul(n_comp as libc::c_ulong);
    let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
    let mut len: size_t = leng
        .wrapping_add(lenl)
        .wrapping_add(lenk)
        .wrapping_add(lenj)
        .wrapping_add(leni)
        .wrapping_add(len0);
    let mut g: *mut f64 = 0 as *mut f64;
    g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut f64;
    cache = g.offset(len as isize);
    let mut g1: *mut f64 = g.offset(leng as isize);
    let mut gout: *mut f64 = 0 as *mut f64;
    let mut gctri: *mut f64 = 0 as *mut f64;
    let mut gctrj: *mut f64 = 0 as *mut f64;
    let mut gctrk: *mut f64 = 0 as *mut f64;
    let mut gctrl: *mut f64 = 0 as *mut f64;
    if n_comp == 1 as i32 {
        gctrl = gctr;
        lempty = empty;
    } else {
        gctrl = g1;
        g1 = g1.offset(lenl as isize);
    }
    if l_ctr == 1 as i32 {
        gctrk = gctrl;
        kempty = lempty;
    } else {
        gctrk = g1;
        g1 = g1.offset(lenk as isize);
    }
    if k_ctr == 1 as i32 {
        gctrj = gctrk;
        jempty = kempty;
    } else {
        gctrj = g1;
        g1 = g1.offset(lenj as isize);
    }
    if j_ctr == 1 as i32 {
        gctri = gctrj;
        iempty = jempty;
    } else {
        gctri = g1;
        g1 = g1.offset(leni as isize);
    }
    if i_ctr == 1 as i32 {
        gout = gctri;
        gempty = iempty;
    } else {
        gout = g1;
        g1 = g1.offset(leng as isize);
    }
    lp = 0 as i32;
    while lp < l_prim {
        (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
        if l_ctr == 1 as i32 {
            fac1l = (*envs).common_factor * *cl.offset(lp as isize);
        } else {
            fac1l = (*envs).common_factor;
            *kempty = 1 as i32;
        }
        kp = 0 as i32;
        while kp < k_prim {
            akl = *ak.offset(kp as isize) + *al.offset(lp as isize);
            ekl = rr_kl * *ak.offset(kp as isize) * *al.offset(lp as isize) / akl;
            ccekl = ekl - log_rr_kl - *log_maxck.offset(kp as isize)
                - *log_maxcl.offset(lp as isize);
            if !(ccekl > expcutoff) {
                (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
                rkl[0 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(0 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(0 as i32 as isize))
                    / akl;
                rkl[1 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(1 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(1 as i32 as isize))
                    / akl;
                rkl[2 as i32
                    as usize] = (*ak.offset(kp as isize)
                    * *rk.offset(2 as i32 as isize)
                    + *al.offset(lp as isize) * *rl.offset(2 as i32 as isize))
                    / akl;
                eijcutoff = expcutoff - ccekl;
                ekl = (-ekl).exp(); // exp(-ekl);
                if k_ctr == 1 as i32 {
                    fac1k = fac1l * *ck.offset(kp as isize);
                } else {
                    fac1k = fac1l;
                    *jempty = 1 as i32;
                }
                pdata_ij = pdata_base;
                jp = 0 as i32;
                while jp < j_prim {
                    (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
                    if j_ctr == 1 as i32 {
                        fac1j = fac1k * *cj.offset(jp as isize);
                    } else {
                        fac1j = fac1k;
                        *iempty = 1 as i32;
                    }
                    ip = 0 as i32;
                    while ip < i_prim {
                        if !((*pdata_ij).cceij > eijcutoff) {
                            (*envs)
                                .ai[0 as i32 as usize] = *ai.offset(ip as isize);
                            rij = ((*pdata_ij).rij).as_mut_ptr();
                            cutoff = eijcutoff - (*pdata_ij).cceij;
                            expijkl = (*pdata_ij).eij * ekl;
                            if i_ctr == 1 as i32 {
                                fac1i = fac1j * *ci.offset(ip as isize) * expijkl;
                            } else {
                                fac1i = fac1j * expijkl;
                            }
                            (*envs).fac[0 as i32 as usize] = fac1i;
                            if ::core::mem::transmute::<
                                _,
                                fn(_, _, _, _, _) -> i32,
                            >(
                                (Some(
                                    ((*envs).f_g0_2e).expect("non-null function pointer"),
                                ))
                                    .expect("non-null function pointer"),
                            )(g, rij, rkl.as_mut_ptr(), cutoff, envs) != 0
                            {
                                ::core::mem::transmute::<
                                    _,
                                    fn(_, _, _, _, _),
                                >(
                                    (Some(((*envs).f_gout).expect("non-null function pointer")))
                                        .expect("non-null function pointer"),
                                )(gout, g, idx, envs, *gempty);
                                if i_ctr > 1 as i32 {
                                    if *iempty != 0 {
                                        CINTprim_to_ctr_0(
                                            gctri,
                                            gout,
                                            ci.offset(ip as isize),
                                            len0,
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
                                            len0,
                                            i_prim,
                                            i_ctr,
                                            *non0ctri.offset(ip as isize),
                                            non0idxi.offset((ip * i_ctr) as isize),
                                        );
                                    }
                                }
                                *iempty = 0 as i32;
                            }
                        }
                        ip += 1;
                        pdata_ij = pdata_ij.offset(1);
                    }
                    if *iempty == 0 {
                        if j_ctr > 1 as i32 {
                            if *jempty != 0 {
                                CINTprim_to_ctr_0(
                                    gctrj,
                                    gctri,
                                    cj.offset(jp as isize),
                                    leni,
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
                                    leni,
                                    j_prim,
                                    j_ctr,
                                    *non0ctrj.offset(jp as isize),
                                    non0idxj.offset((jp * j_ctr) as isize),
                                );
                            }
                        }
                        *jempty = 0 as i32;
                    }
                    jp += 1;
                }
                if *jempty == 0 {
                    if k_ctr > 1 as i32 {
                        if *kempty != 0 {
                            CINTprim_to_ctr_0(
                                gctrk,
                                gctrj,
                                ck.offset(kp as isize),
                                lenj,
                                k_prim,
                                k_ctr,
                                *non0ctrk.offset(kp as isize),
                                non0idxk.offset((kp * k_ctr) as isize),
                            );
                        } else {
                            CINTprim_to_ctr_1(
                                gctrk,
                                gctrj,
                                ck.offset(kp as isize),
                                lenj,
                                k_prim,
                                k_ctr,
                                *non0ctrk.offset(kp as isize),
                                non0idxk.offset((kp * k_ctr) as isize),
                            );
                        }
                    }
                    *kempty = 0 as i32;
                }
            }
            kp += 1;
        }
        if *kempty == 0 {
            if l_ctr > 1 as i32 {
                if *lempty != 0 {
                    CINTprim_to_ctr_0(
                        gctrl,
                        gctrk,
                        cl.offset(lp as isize),
                        lenk,
                        l_prim,
                        l_ctr,
                        *non0ctrl.offset(lp as isize),
                        non0idxl.offset((lp * l_ctr) as isize),
                    );
                } else {
                    CINTprim_to_ctr_1(
                        gctrl,
                        gctrk,
                        cl.offset(lp as isize),
                        lenk,
                        l_prim,
                        l_ctr,
                        *non0ctrl.offset(lp as isize),
                        non0idxl.offset((lp * l_ctr) as isize),
                    );
                }
            }
            *lempty = 0 as i32;
        }
        lp += 1;
    }
    if n_comp > 1 as i32 && *lempty == 0 {
        if *empty != 0 {
            CINTdmat_transpose(
                gctr,
                gctrl,
                nf.wrapping_mul(nc as libc::c_ulong) as i32,
                n_comp,
            );
            *empty = 0 as i32;
        } else {
            CINTdplus_transpose(
                gctr,
                gctrl,
                nf.wrapping_mul(nc as libc::c_ulong) as i32,
                n_comp,
            );
        }
    }
    return (*empty == 0) as i32;
}
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_1111_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = 1 as i32;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng.wrapping_add(len0);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     if n_comp == 1 as i32 {
//         gout = gctr;
//         gempty = empty;
//     } else {
//         gout = g.offset(leng as isize);
//     }
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         fac1l = (*envs).common_factor * *cl.offset(lp as isize);
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 fac1k = fac1l * *ck.offset(kp as isize);
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 pdata_ij = _pdata_ij;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     fac1j = fac1k * *cj.offset(jp as isize);
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             (*envs)
//                                 .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                             expij = (*pdata_ij).eij;
//                             rij = ((*pdata_ij).rij).as_mut_ptr();
//                             fac1i = fac1j * *ci.offset(ip as isize) * expij * expkl;
//                             (*envs).fac[0 as i32 as usize] = fac1i;
//                             cutoff = eijcutoff - (*pdata_ij).cceij;
//                             if ::core::mem::transmute::<
//                                 _,
//                                 fn(_, _, _, _, _) -> i32,
//                             >(
//                                 (Some(
//                                     ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                 ))
//                                     .expect("non-null function pointer"),
//                             )(g, rij, rkl, cutoff, envs) != 0
//                             {
//                                 ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _),
//                                 >(
//                                     (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                         .expect("non-null function pointer"),
//                                 )(gout, g, idx, envs, *gempty);
//                                 *gempty = 0 as i32;
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     jp += 1;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *gempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gout,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gout,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_n111_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = i_ctr;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut leni: size_t = nf
//         .wrapping_mul(i_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng.wrapping_add(leni).wrapping_add(len0);
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     let mut g1: *mut f64 = g.offset(leng as isize);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut gctri: *mut f64 = 0 as *mut f64;
//     if n_comp == 1 as i32 {
//         gctri = gctr;
//         iempty = empty;
//     } else {
//         gctri = g1;
//         g1 = g1.offset(leni as isize);
//     }
//     gout = g1;
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         fac1l = (*envs).common_factor * *cl.offset(lp as isize);
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 fac1k = fac1l * *ck.offset(kp as isize);
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 pdata_ij = _pdata_ij;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     fac1j = fac1k * *cj.offset(jp as isize);
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             if !((*pdata_ij).cceij > eijcutoff) {
//                                 (*envs)
//                                     .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                                 expij = (*pdata_ij).eij;
//                                 rij = ((*pdata_ij).rij).as_mut_ptr();
//                                 cutoff = eijcutoff - (*pdata_ij).cceij;
//                                 fac1i = fac1j * expij * expkl;
//                                 (*envs).fac[0 as i32 as usize] = fac1i;
//                                 if ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _) -> i32,
//                                 >(
//                                     (Some(
//                                         ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                     ))
//                                         .expect("non-null function pointer"),
//                                 )(g, rij, rkl, cutoff, envs) != 0
//                                 {
//                                     ::core::mem::transmute::<
//                                         _,
//                                         fn(_, _, _, _, _),
//                                     >(
//                                         (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                             .expect("non-null function pointer"),
//                                     )(gout, g, idx, envs, 1 as i32);
//                                     if i_ctr > 1 as i32 {
//                                         if *iempty != 0 {
//                                             CINTprim_to_ctr_0(
//                                                 gctri,
//                                                 gout,
//                                                 ci.offset(ip as isize),
//                                                 len0,
//                                                 i_prim,
//                                                 i_ctr,
//                                                 *non0ctri.offset(ip as isize),
//                                                 non0idxi.offset((ip * i_ctr) as isize),
//                                             );
//                                         } else {
//                                             CINTprim_to_ctr_1(
//                                                 gctri,
//                                                 gout,
//                                                 ci.offset(ip as isize),
//                                                 len0,
//                                                 i_prim,
//                                                 i_ctr,
//                                                 *non0ctri.offset(ip as isize),
//                                                 non0idxi.offset((ip * i_ctr) as isize),
//                                             );
//                                         }
//                                     }
//                                     *iempty = 0 as i32;
//                                 }
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     jp += 1;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *iempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gctri,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gctri,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_1n11_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = j_ctr;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut lenj: size_t = nf
//         .wrapping_mul(j_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng.wrapping_add(lenj).wrapping_add(len0);
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     let mut g1: *mut f64 = g.offset(leng as isize);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut gctrj: *mut f64 = 0 as *mut f64;
//     if n_comp == 1 as i32 {
//         gctrj = gctr;
//         jempty = empty;
//     } else {
//         gctrj = g1;
//         g1 = g1.offset(lenj as isize);
//     }
//     gout = g1;
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         fac1l = (*envs).common_factor * *cl.offset(lp as isize);
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 fac1k = fac1l * *ck.offset(kp as isize);
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 pdata_ij = _pdata_ij;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     fac1j = fac1k;
//                     *iempty = 1 as i32;
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             (*envs)
//                                 .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                             expij = (*pdata_ij).eij;
//                             rij = ((*pdata_ij).rij).as_mut_ptr();
//                             cutoff = eijcutoff - (*pdata_ij).cceij;
//                             fac1i = fac1j * *ci.offset(ip as isize) * expij * expkl;
//                             (*envs).fac[0 as i32 as usize] = fac1i;
//                             if ::core::mem::transmute::<
//                                 _,
//                                 fn(_, _, _, _, _) -> i32,
//                             >(
//                                 (Some(
//                                     ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                 ))
//                                     .expect("non-null function pointer"),
//                             )(g, rij, rkl, cutoff, envs) != 0
//                             {
//                                 ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _),
//                                 >(
//                                     (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                         .expect("non-null function pointer"),
//                                 )(gout, g, idx, envs, *iempty);
//                                 *iempty = 0 as i32;
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     if *iempty == 0 {
//                         if j_ctr > 1 as i32 {
//                             if *jempty != 0 {
//                                 CINTprim_to_ctr_0(
//                                     gctrj,
//                                     gout,
//                                     cj.offset(jp as isize),
//                                     len0,
//                                     j_prim,
//                                     j_ctr,
//                                     *non0ctrj.offset(jp as isize),
//                                     non0idxj.offset((jp * j_ctr) as isize),
//                                 );
//                             } else {
//                                 CINTprim_to_ctr_1(
//                                     gctrj,
//                                     gout,
//                                     cj.offset(jp as isize),
//                                     len0,
//                                     j_prim,
//                                     j_ctr,
//                                     *non0ctrj.offset(jp as isize),
//                                     non0idxj.offset((jp * j_ctr) as isize),
//                                 );
//                             }
//                         }
//                         *jempty = 0 as i32;
//                     }
//                     jp += 1;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *jempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gctrj,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gctrj,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_11n1_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = k_ctr;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut lenk: size_t = nf
//         .wrapping_mul(k_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng.wrapping_add(lenk).wrapping_add(len0);
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     let mut g1: *mut f64 = g.offset(leng as isize);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut gctrk: *mut f64 = 0 as *mut f64;
//     if n_comp == 1 as i32 {
//         gctrk = gctr;
//         kempty = empty;
//     } else {
//         gctrk = g1;
//         g1 = g1.offset(lenk as isize);
//     }
//     gout = g1;
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         fac1l = (*envs).common_factor * *cl.offset(lp as isize);
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 fac1k = fac1l;
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 pdata_ij = _pdata_ij;
//                 *jempty = 1 as i32;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     fac1j = fac1k * *cj.offset(jp as isize);
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             (*envs)
//                                 .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                             expij = (*pdata_ij).eij;
//                             rij = ((*pdata_ij).rij).as_mut_ptr();
//                             cutoff = eijcutoff - (*pdata_ij).cceij;
//                             fac1i = fac1j * *ci.offset(ip as isize) * expij * expkl;
//                             (*envs).fac[0 as i32 as usize] = fac1i;
//                             if ::core::mem::transmute::<
//                                 _,
//                                 fn(_, _, _, _, _) -> i32,
//                             >(
//                                 (Some(
//                                     ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                 ))
//                                     .expect("non-null function pointer"),
//                             )(g, rij, rkl, cutoff, envs) != 0
//                             {
//                                 ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _),
//                                 >(
//                                     (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                         .expect("non-null function pointer"),
//                                 )(gout, g, idx, envs, *jempty);
//                                 *jempty = 0 as i32;
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     jp += 1;
//                 }
//                 if *jempty == 0 {
//                     if k_ctr > 1 as i32 {
//                         if *kempty != 0 {
//                             CINTprim_to_ctr_0(
//                                 gctrk,
//                                 gout,
//                                 ck.offset(kp as isize),
//                                 len0,
//                                 k_prim,
//                                 k_ctr,
//                                 *non0ctrk.offset(kp as isize),
//                                 non0idxk.offset((kp * k_ctr) as isize),
//                             );
//                         } else {
//                             CINTprim_to_ctr_1(
//                                 gctrk,
//                                 gout,
//                                 ck.offset(kp as isize),
//                                 len0,
//                                 k_prim,
//                                 k_ctr,
//                                 *non0ctrk.offset(kp as isize),
//                                 non0idxk.offset((kp * k_ctr) as isize),
//                             );
//                         }
//                     }
//                     *kempty = 0 as i32;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *kempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gctrk,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gctrk,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_111n_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = l_ctr;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut lenl: size_t = nf
//         .wrapping_mul(l_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng.wrapping_add(lenl).wrapping_add(len0);
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     let mut g1: *mut f64 = g.offset(leng as isize);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut gctrl: *mut f64 = 0 as *mut f64;
//     if n_comp == 1 as i32 {
//         gctrl = gctr;
//         lempty = empty;
//     } else {
//         gctrl = g1;
//         g1 = g1.offset(lenl as isize);
//     }
//     gout = g1;
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         fac1l = (*envs).common_factor;
//         *kempty = 1 as i32;
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 fac1k = fac1l * *ck.offset(kp as isize);
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 pdata_ij = _pdata_ij;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     fac1j = fac1k * *cj.offset(jp as isize);
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             (*envs)
//                                 .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                             expij = (*pdata_ij).eij;
//                             rij = ((*pdata_ij).rij).as_mut_ptr();
//                             cutoff = eijcutoff - (*pdata_ij).cceij;
//                             fac1i = fac1j * *ci.offset(ip as isize) * expij * expkl;
//                             (*envs).fac[0 as i32 as usize] = fac1i;
//                             if ::core::mem::transmute::<
//                                 _,
//                                 fn(_, _, _, _, _) -> i32,
//                             >(
//                                 (Some(
//                                     ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                 ))
//                                     .expect("non-null function pointer"),
//                             )(g, rij, rkl, cutoff, envs) != 0
//                             {
//                                 ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _),
//                                 >(
//                                     (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                         .expect("non-null function pointer"),
//                                 )(gout, g, idx, envs, *kempty);
//                                 *kempty = 0 as i32;
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     jp += 1;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         if *kempty == 0 {
//             if l_ctr > 1 as i32 {
//                 if *lempty != 0 {
//                     CINTprim_to_ctr_0(
//                         gctrl,
//                         gout,
//                         cl.offset(lp as isize),
//                         len0,
//                         l_prim,
//                         l_ctr,
//                         *non0ctrl.offset(lp as isize),
//                         non0idxl.offset((lp * l_ctr) as isize),
//                     );
//                 } else {
//                     CINTprim_to_ctr_1(
//                         gctrl,
//                         gout,
//                         cl.offset(lp as isize),
//                         len0,
//                         l_prim,
//                         l_ctr,
//                         *non0ctrl.offset(lp as isize),
//                         non0idxl.offset((lp * l_ctr) as isize),
//                     );
//                 }
//             }
//             *lempty = 0 as i32;
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *lempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gctrl,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gctrl,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// #[no_mangle]
// pub unsafe extern "C" fn CINT2e_loop(
//     mut gctr: *mut f64,
//     mut envs: *mut CINTEnvVars,
//     mut cache: *mut f64,
//     mut empty: *mut i32,
// ) -> i32 {
//     let mut shls: *mut i32 = (*envs).shls;
//     let mut bas: *mut i32 = (*envs).bas;
//     let mut env: *mut f64 = (*envs).env;
//     let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
//     let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
//     let mut k_sh: i32 = *shls.offset(2 as i32 as isize);
//     let mut l_sh: i32 = *shls.offset(3 as i32 as isize);
//     let mut opt: *mut CINTOpt = (*envs).opt;
//     if !((*opt).pairdata).is_null()
//         && (*((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize)
//             == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void as *mut PairData
//             || *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize)
//                 == 0xffffffffffffffff as libc::c_ulong as *mut libc::c_void
//                     as *mut PairData)
//     {
//         return 0 as i32;
//     }
//     let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
//     let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
//     let mut k_ctr: i32 = (*envs).x_ctr[2 as i32 as usize];
//     let mut l_ctr: i32 = (*envs).x_ctr[3 as i32 as usize];
//     let mut i_prim: i32 = *bas
//         .offset((8 as i32 * i_sh + 2 as i32) as isize);
//     let mut j_prim: i32 = *bas
//         .offset((8 as i32 * j_sh + 2 as i32) as isize);
//     let mut k_prim: i32 = *bas
//         .offset((8 as i32 * k_sh + 2 as i32) as isize);
//     let mut l_prim: i32 = *bas
//         .offset((8 as i32 * l_sh + 2 as i32) as isize);
//     let mut ai: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize,
//         );
//     let mut aj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ak: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 5 as i32) as isize) as isize,
//         );
//     let mut al: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 5 as i32) as isize) as isize,
//         );
//     let mut ci: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cj: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize,
//         );
//     let mut ck: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * k_sh + 6 as i32) as isize) as isize,
//         );
//     let mut cl: *mut f64 = env
//         .offset(
//             *bas.offset((8 as i32 * l_sh + 6 as i32) as isize) as isize,
//         );
//     let mut expcutoff: f64 = (*envs).expcutoff;
//     let mut rr_ij: f64 = (*envs).rirj[0 as i32 as usize]
//         * (*envs).rirj[0 as i32 as usize]
//         + (*envs).rirj[1 as i32 as usize]
//             * (*envs).rirj[1 as i32 as usize]
//         + (*envs).rirj[2 as i32 as usize]
//             * (*envs).rirj[2 as i32 as usize];
//     let mut rr_kl: f64 = (*envs).rkrl[0 as i32 as usize]
//         * (*envs).rkrl[0 as i32 as usize]
//         + (*envs).rkrl[1 as i32 as usize]
//             * (*envs).rkrl[1 as i32 as usize]
//         + (*envs).rkrl[2 as i32 as usize]
//             * (*envs).rkrl[2 as i32 as usize];
//     let mut _pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut _pdata_kl: *mut PairData = 0 as *mut PairData;
//     let mut pdata_ij: *mut PairData = 0 as *mut PairData;
//     let mut pdata_kl: *mut PairData = 0 as *mut PairData;
//     if !((*opt).pairdata).is_null() {
//         _pdata_ij = *((*opt).pairdata).offset((i_sh * (*opt).nbas + j_sh) as isize);
//         _pdata_kl = *((*opt).pairdata).offset((k_sh * (*opt).nbas + l_sh) as isize);
//     } else {
//         let mut log_maxci: *mut f64 = *((*opt).log_max_coeff)
//             .offset(i_sh as isize);
//         let mut log_maxcj: *mut f64 = *((*opt).log_max_coeff)
//             .offset(j_sh as isize);
//         _pdata_ij = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut PairData;
//         cache = _pdata_ij.offset((i_prim * j_prim + k_prim * l_prim) as isize)
//             as *mut f64;
//         if CINTset_pairdata(
//             _pdata_ij,
//             ai,
//             aj,
//             (*envs).ri,
//             (*envs).rj,
//             log_maxci,
//             log_maxcj,
//             (*envs).li_ceil,
//             (*envs).lj_ceil,
//             i_prim,
//             j_prim,
//             rr_ij,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//         let mut log_maxck: *mut f64 = *((*opt).log_max_coeff)
//             .offset(k_sh as isize);
//         let mut log_maxcl: *mut f64 = *((*opt).log_max_coeff)
//             .offset(l_sh as isize);
//         _pdata_kl = _pdata_ij.offset((i_prim * j_prim) as isize);
//         if CINTset_pairdata(
//             _pdata_kl,
//             ak,
//             al,
//             (*envs).rk,
//             (*envs).c2rust_unnamed_1.rl,
//             log_maxck,
//             log_maxcl,
//             (*envs).lk_ceil,
//             (*envs).ll_ceil,
//             k_prim,
//             l_prim,
//             rr_kl,
//             expcutoff,
//             env,
//         ) != 0
//         {
//             return 0 as i32;
//         }
//     }
//     let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_e2
//         * (*envs).ncomp_tensor;
//     let mut nf: size_t = (*envs).nf as size_t;
//     let mut fac1i: f64 = 0.;
//     let mut fac1j: f64 = 0.;
//     let mut fac1k: f64 = 0.;
//     let mut fac1l: f64 = 0.;
//     let mut ip: i32 = 0;
//     let mut jp: i32 = 0;
//     let mut kp: i32 = 0;
//     let mut lp: i32 = 0;
//     let mut _empty: [i32; 5] = [
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     let mut iempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(0 as i32 as isize);
//     let mut jempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(1 as i32 as isize);
//     let mut kempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(2 as i32 as isize);
//     let mut lempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(3 as i32 as isize);
//     let mut gempty: *mut i32 = _empty
//         .as_mut_ptr()
//         .offset(4 as i32 as isize);
//     let mut non0ctri: *mut i32 = *((*opt).non0ctr).offset(i_sh as isize);
//     let mut non0ctrj: *mut i32 = *((*opt).non0ctr).offset(j_sh as isize);
//     let mut non0ctrk: *mut i32 = *((*opt).non0ctr).offset(k_sh as isize);
//     let mut non0ctrl: *mut i32 = *((*opt).non0ctr).offset(l_sh as isize);
//     let mut non0idxi: *mut i32 = *((*opt).sortedidx).offset(i_sh as isize);
//     let mut non0idxj: *mut i32 = *((*opt).sortedidx).offset(j_sh as isize);
//     let mut non0idxk: *mut i32 = *((*opt).sortedidx).offset(k_sh as isize);
//     let mut non0idxl: *mut i32 = *((*opt).sortedidx).offset(l_sh as isize);
//     let mut expij: f64 = 0.;
//     let mut expkl: f64 = 0.;
//     let mut eijcutoff: f64 = 0.;
//     let mut eklcutoff: f64 = 0.;
//     let mut cutoff: f64 = 0.;
//     eklcutoff = expcutoff;
//     let mut rij: *mut f64 = 0 as *mut f64;
//     let mut rkl: *mut f64 = 0 as *mut f64;
//     let mut idx: *mut i32 = *((*opt).index_xyz_array)
//         .offset(
//             ((*envs).i_l * 16 as i32 * 16 as i32 * 16 as i32
//                 + (*envs).j_l * 16 as i32 * 16 as i32
//                 + (*envs).k_l * 16 as i32 + (*envs).l_l) as isize,
//         );
//     if idx.is_null() {
//         idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//             & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//             as *mut i32;
//         cache = idx.offset(nf.wrapping_mul(3 as i32 as libc::c_ulong) as isize)
//             as *mut f64;
//         CINTg2e_index_xyz(idx, envs);
//     }
//     let mut omega: f64 = *env.offset(8 as i32 as isize);
//     if omega < 0 as i32 as f64 && (*envs).rys_order > 1 as i32
//     {
//         let mut r_guess: f64 = 8.0f64;
//         let mut omega2: f64 = omega * omega;
//         let mut lij: i32 = (*envs).li_ceil + (*envs).lj_ceil;
//         let mut lkl: i32 = (*envs).lk_ceil + (*envs).ll_ceil;
//         if lij > 0 as i32 {
//             let mut dist_ij: f64 = rr_ij.sqrt();
//             let mut aij: f64 = *ai
//                 .offset((i_prim - 1 as i32) as isize)
//                 + *aj.offset((j_prim - 1 as i32) as isize);
//             let mut theta: f64 = omega2 / (omega2 + aij);
//             expcutoff
//                 += lij as f64
//                     * ((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64)).ln(); //log((dist_ij + theta * r_guess + 1.0f64) / (dist_ij + 1.0f64));
//         }
//         if lkl > 0 as i32 {
//             let mut dist_kl: f64 = rr_kl.sqrt();
//             let mut akl: f64 = *ak
//                 .offset((k_prim - 1 as i32) as isize)
//                 + *al.offset((l_prim - 1 as i32) as isize);
//             let mut theta_0: f64 = omega2 / (omega2 + akl);
//             expcutoff
//                 += lkl as f64
//                     * ((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64)).ln(); //log((dist_kl + theta_0 * r_guess + 1.0f64) / (dist_kl + 1.0f64));
//         }
//     }
//     let mut nc: i32 = i_ctr * j_ctr * k_ctr * l_ctr;
//     let mut leng: size_t = ((*envs).g_size * 3 as i32
//         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
//     let mut lenl: size_t = nf
//         .wrapping_mul(nc as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut lenk: size_t = nf
//         .wrapping_mul(i_ctr as libc::c_ulong)
//         .wrapping_mul(j_ctr as libc::c_ulong)
//         .wrapping_mul(k_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut lenj: size_t = nf
//         .wrapping_mul(i_ctr as libc::c_ulong)
//         .wrapping_mul(j_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut leni: size_t = nf
//         .wrapping_mul(i_ctr as libc::c_ulong)
//         .wrapping_mul(n_comp as libc::c_ulong);
//     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
//     let mut len: size_t = leng
//         .wrapping_add(lenl)
//         .wrapping_add(lenk)
//         .wrapping_add(lenj)
//         .wrapping_add(leni)
//         .wrapping_add(len0);
//     let mut g: *mut f64 = 0 as *mut f64;
//     g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
//         & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
//         as *mut f64;
//     cache = g.offset(len as isize);
//     let mut g1: *mut f64 = g.offset(leng as isize);
//     let mut gout: *mut f64 = 0 as *mut f64;
//     let mut gctri: *mut f64 = 0 as *mut f64;
//     let mut gctrj: *mut f64 = 0 as *mut f64;
//     let mut gctrk: *mut f64 = 0 as *mut f64;
//     let mut gctrl: *mut f64 = 0 as *mut f64;
//     if n_comp == 1 as i32 {
//         gctrl = gctr;
//         lempty = empty;
//     } else {
//         gctrl = g1;
//         g1 = g1.offset(lenl as isize);
//     }
//     if l_ctr == 1 as i32 {
//         gctrk = gctrl;
//         kempty = lempty;
//     } else {
//         gctrk = g1;
//         g1 = g1.offset(lenk as isize);
//     }
//     if k_ctr == 1 as i32 {
//         gctrj = gctrk;
//         jempty = kempty;
//     } else {
//         gctrj = g1;
//         g1 = g1.offset(lenj as isize);
//     }
//     if j_ctr == 1 as i32 {
//         gctri = gctrj;
//         iempty = jempty;
//     } else {
//         gctri = g1;
//         g1 = g1.offset(leni as isize);
//     }
//     if i_ctr == 1 as i32 {
//         gout = gctri;
//         gempty = iempty;
//     } else {
//         gout = g1;
//         g1 = g1.offset(leng as isize);
//     }
//     pdata_kl = _pdata_kl;
//     lp = 0 as i32;
//     while lp < l_prim {
//         (*envs).al[0 as i32 as usize] = *al.offset(lp as isize);
//         if l_ctr == 1 as i32 {
//             fac1l = (*envs).common_factor * *cl.offset(lp as isize);
//         } else {
//             fac1l = (*envs).common_factor;
//             *kempty = 1 as i32;
//         }
//         kp = 0 as i32;
//         while kp < k_prim {
//             if !((*pdata_kl).cceij > eklcutoff) {
//                 (*envs).ak[0 as i32 as usize] = *ak.offset(kp as isize);
//                 expkl = (*pdata_kl).eij;
//                 rkl = ((*pdata_kl).rij).as_mut_ptr();
//                 eijcutoff = eklcutoff - (*pdata_kl).cceij;
//                 if k_ctr == 1 as i32 {
//                     fac1k = fac1l * *ck.offset(kp as isize);
//                 } else {
//                     fac1k = fac1l;
//                     *jempty = 1 as i32;
//                 }
//                 pdata_ij = _pdata_ij;
//                 jp = 0 as i32;
//                 while jp < j_prim {
//                     (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
//                     if j_ctr == 1 as i32 {
//                         fac1j = fac1k * *cj.offset(jp as isize);
//                     } else {
//                         fac1j = fac1k;
//                         *iempty = 1 as i32;
//                     }
//                     ip = 0 as i32;
//                     while ip < i_prim {
//                         if !((*pdata_ij).cceij > eijcutoff) {
//                             (*envs)
//                                 .ai[0 as i32 as usize] = *ai.offset(ip as isize);
//                             expij = (*pdata_ij).eij;
//                             rij = ((*pdata_ij).rij).as_mut_ptr();
//                             cutoff = eijcutoff - (*pdata_ij).cceij;
//                             if i_ctr == 1 as i32 {
//                                 fac1i = fac1j * *ci.offset(ip as isize) * expij * expkl;
//                             } else {
//                                 fac1i = fac1j * expij * expkl;
//                             }
//                             (*envs).fac[0 as i32 as usize] = fac1i;
//                             if ::core::mem::transmute::<
//                                 _,
//                                 fn(_, _, _, _, _) -> i32,
//                             >(
//                                 (Some(
//                                     ((*envs).f_g0_2e).expect("non-null function pointer"),
//                                 ))
//                                     .expect("non-null function pointer"),
//                             )(g, rij, rkl, cutoff, envs) != 0
//                             {
//                                 ::core::mem::transmute::<
//                                     _,
//                                     fn(_, _, _, _, _),
//                                 >(
//                                     (Some(((*envs).f_gout).expect("non-null function pointer")))
//                                         .expect("non-null function pointer"),
//                                 )(gout, g, idx, envs, *gempty);
//                                 if i_ctr > 1 as i32 {
//                                     if *iempty != 0 {
//                                         CINTprim_to_ctr_0(
//                                             gctri,
//                                             gout,
//                                             ci.offset(ip as isize),
//                                             len0,
//                                             i_prim,
//                                             i_ctr,
//                                             *non0ctri.offset(ip as isize),
//                                             non0idxi.offset((ip * i_ctr) as isize),
//                                         );
//                                     } else {
//                                         CINTprim_to_ctr_1(
//                                             gctri,
//                                             gout,
//                                             ci.offset(ip as isize),
//                                             len0,
//                                             i_prim,
//                                             i_ctr,
//                                             *non0ctri.offset(ip as isize),
//                                             non0idxi.offset((ip * i_ctr) as isize),
//                                         );
//                                     }
//                                 }
//                                 *iempty = 0 as i32;
//                             }
//                         }
//                         ip += 1;
//                         pdata_ij = pdata_ij.offset(1);
//                     }
//                     if *iempty == 0 {
//                         if j_ctr > 1 as i32 {
//                             if *jempty != 0 {
//                                 CINTprim_to_ctr_0(
//                                     gctrj,
//                                     gctri,
//                                     cj.offset(jp as isize),
//                                     leni,
//                                     j_prim,
//                                     j_ctr,
//                                     *non0ctrj.offset(jp as isize),
//                                     non0idxj.offset((jp * j_ctr) as isize),
//                                 );
//                             } else {
//                                 CINTprim_to_ctr_1(
//                                     gctrj,
//                                     gctri,
//                                     cj.offset(jp as isize),
//                                     leni,
//                                     j_prim,
//                                     j_ctr,
//                                     *non0ctrj.offset(jp as isize),
//                                     non0idxj.offset((jp * j_ctr) as isize),
//                                 );
//                             }
//                         }
//                         *jempty = 0 as i32;
//                     }
//                     jp += 1;
//                 }
//                 if *jempty == 0 {
//                     if k_ctr > 1 as i32 {
//                         if *kempty != 0 {
//                             CINTprim_to_ctr_0(
//                                 gctrk,
//                                 gctrj,
//                                 ck.offset(kp as isize),
//                                 lenj,
//                                 k_prim,
//                                 k_ctr,
//                                 *non0ctrk.offset(kp as isize),
//                                 non0idxk.offset((kp * k_ctr) as isize),
//                             );
//                         } else {
//                             CINTprim_to_ctr_1(
//                                 gctrk,
//                                 gctrj,
//                                 ck.offset(kp as isize),
//                                 lenj,
//                                 k_prim,
//                                 k_ctr,
//                                 *non0ctrk.offset(kp as isize),
//                                 non0idxk.offset((kp * k_ctr) as isize),
//                             );
//                         }
//                     }
//                     *kempty = 0 as i32;
//                 }
//             }
//             kp += 1;
//             pdata_kl = pdata_kl.offset(1);
//         }
//         if *kempty == 0 {
//             if l_ctr > 1 as i32 {
//                 if *lempty != 0 {
//                     CINTprim_to_ctr_0(
//                         gctrl,
//                         gctrk,
//                         cl.offset(lp as isize),
//                         lenk,
//                         l_prim,
//                         l_ctr,
//                         *non0ctrl.offset(lp as isize),
//                         non0idxl.offset((lp * l_ctr) as isize),
//                     );
//                 } else {
//                     CINTprim_to_ctr_1(
//                         gctrl,
//                         gctrk,
//                         cl.offset(lp as isize),
//                         lenk,
//                         l_prim,
//                         l_ctr,
//                         *non0ctrl.offset(lp as isize),
//                         non0idxl.offset((lp * l_ctr) as isize),
//                     );
//                 }
//             }
//             *lempty = 0 as i32;
//         }
//         lp += 1;
//     }
//     if n_comp > 1 as i32 && *lempty == 0 {
//         if *empty != 0 {
//             CINTdmat_transpose(
//                 gctr,
//                 gctrl,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//             *empty = 0 as i32;
//         } else {
//             CINTdplus_transpose(
//                 gctr,
//                 gctrl,
//                 nf.wrapping_mul(nc as libc::c_ulong) as i32,
//                 n_comp,
//             );
//         }
//     }
//     return (*empty == 0) as i32;
// }
// static mut CINTf_2e_loop: [Option::<
//     unsafe extern "C" fn(
//         *mut f64,
//         *mut CINTEnvVars,
//         *mut f64,
//         *mut i32,
//     ) -> i32,
// >; 16] = unsafe {
//     [
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_n111_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_1n11_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_11n1_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_111n_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//         Some(
//             CINT2e_1111_loop
//                 as unsafe extern "C" fn(
//                     *mut f64,
//                     *mut CINTEnvVars,
//                     *mut f64,
//                     *mut i32,
//                 ) -> i32,
//         ),
//     ]
// };
#[no_mangle]
pub fn CINT2e_drv(
    out: &mut [f64],
    mut dims: Vec<i32>,
    envs: &mut CINTEnvVars,
    // mut opt: *mut CINTOpt,
    cache: Vec<f64>,
    f_c2s: Option<F_FC2S>,
) -> i32 {
    let mut x_ctr: [i32; 4] = envs.x_ctr;
    let mut nf: size_t = envs.nf as size_t;
    let mut nc: size_t = nf * (x_ctr[0] * x_ctr[1] * x_ctr[2] * x_ctr[3]) as u64;
    let mut n_comp: i32 = envs.ncomp_e1 * envs.ncomp_e2 * envs.ncomp_tensor;
    // if out.is_null() {
    //     let mut bas: *mut i32 = (*envs).bas.as_mut_ptr();
    //     let mut shls: *mut i32 = (*envs).shls.as_mut_ptr();
    //     let mut i_prim: i32 = *bas
    //         .offset(
    //             (8 as i32 * *shls.offset(0 as i32 as isize)
    //                 + 2 as i32) as isize,
    //         );
    //     let mut j_prim: i32 = *bas
    //         .offset(
    //             (8 as i32 * *shls.offset(1 as i32 as isize)
    //                 + 2 as i32) as isize,
    //         );
    //     let mut k_prim: i32 = *bas
    //         .offset(
    //             (8 as i32 * *shls.offset(2 as i32 as isize)
    //                 + 2 as i32) as isize,
    //         );
    //     let mut l_prim: i32 = *bas
    //         .offset(
    //             (8 as i32 * *shls.offset(3 as i32 as isize)
    //                 + 2 as i32) as isize,
    //         );
    //     let mut pdata_size: size_t = (((i_prim * j_prim + k_prim * l_prim)
    //         * 5 as i32 + i_prim * *x_ctr.offset(0 as i32 as isize)
    //         + j_prim * *x_ctr.offset(1 as i32 as isize)
    //         + k_prim * *x_ctr.offset(2 as i32 as isize)
    //         + l_prim * *x_ctr.offset(3 as i32 as isize)
    //         + (i_prim + j_prim + k_prim + l_prim) * 2 as i32) as libc::c_ulong)
    //         .wrapping_add(nf.wrapping_mul(3 as i32 as libc::c_ulong));
    //     let mut leng: size_t = ((*envs).g_size * 3 as i32
    //         * (((1 as i32) << (*envs).gbits) + 1 as i32)) as size_t;
    //     let mut len0: size_t = nf.wrapping_mul(n_comp as libc::c_ulong);
    //     let mut cache_size: size_t = if leng
    //         .wrapping_add(len0)
    //         .wrapping_add(
    //             nc
    //                 .wrapping_mul(n_comp as libc::c_ulong)
    //                 .wrapping_mul(3 as i32 as libc::c_ulong),
    //         )
    //         .wrapping_add(pdata_size)
    //         > nc
    //             .wrapping_mul(n_comp as libc::c_ulong)
    //             .wrapping_add(nf.wrapping_mul(4 as i32 as libc::c_ulong))
    //     {
    //         leng.wrapping_add(len0)
    //             .wrapping_add(
    //                 nc
    //                     .wrapping_mul(n_comp as libc::c_ulong)
    //                     .wrapping_mul(3 as i32 as libc::c_ulong),
    //             )
    //             .wrapping_add(pdata_size)
    //     } else {
    //         nc.wrapping_mul(n_comp as libc::c_ulong)
    //             .wrapping_add(nf.wrapping_mul(4 as i32 as libc::c_ulong))
    //     };
    //     if cache_size >= 2147483647 as i32 as libc::c_ulong {
    //         println!("CINT2e_drv cache_size overflow: cache_size {} > {}, nf {}, nc {}, n_comp {}\n", 
    //             cache_size, 2147483647, nf, nc, n_comp);
    //         // fprintf(
    //         //     stderr,
    //         //     b"CINT2e_drv cache_size overflow: cache_size %zu > %d, nf %zu, nc %zu, n_comp %d\n\0"
    //         //         as *const u8 as *const libc::c_char,
    //         //     cache_size,
    //         //     2147483647 as i32,
    //         //     nf,
    //         //     nc,
    //         //     n_comp,
    //         // );
    //         cache_size = 0 as i32 as size_t;
    //     }
    //     return cache_size as i32;
    // }
    // let mut stack: &mut [f64]; // = 0 as *mut f64;
    // if cache.is_null() {
    let mut bas_0: &mut [i32] = &mut envs.bas;
    let mut shls_0: [i32; 4] = envs.shls;
    let mut i_prim_0: i32 = bas_0[8 * shls_0[0] as usize + 2];
    let mut j_prim_0: i32 = bas_0[8 * shls_0[1] as usize + 2];
    let mut k_prim_0: i32 = bas_0[8 * shls_0[2] as usize + 2];
    let mut l_prim_0: i32 = bas_0[8 * shls_0[3] as usize + 2];
    let mut pdata_size_0: size_t = (((i_prim_0 * j_prim_0 + k_prim_0 * l_prim_0) * 5 as i32 
        + i_prim_0 * x_ctr[0]
        + j_prim_0 * x_ctr[1]
        + k_prim_0 * x_ctr[2]
        + l_prim_0 * x_ctr[3]
        + (i_prim_0 + j_prim_0 + k_prim_0 + l_prim_0) * 2 as i32)
        as u64) + nf * 3;
    let mut leng_0: size_t = (envs.g_size * 3 * ((1 << envs.gbits) + 1)) as size_t;
    let mut len0_0: size_t = nf * n_comp as u64;
    let mut cache_size_0: size_t = 
        if leng_0 + len0_0 + nc * (n_comp as u64) * 3 + pdata_size_0
        > nc * (n_comp as u64) + nf * 4
        {
            leng_0 + len0_0 + nc * (n_comp as u64) * 3 + pdata_size_0
        } else {
            nc * (n_comp as u64) + nf * 4
        };
    // stack = malloc(
    //     (::core::mem::size_of::<f64>() as libc::c_ulong)
    //         .wrapping_mul(cache_size_0),
    // ) as *mut f64;
    // }
    // let gctr: &mut [f64] = &mut cache_full[..];
    // let cache: &mut [f64] = &mut cache_full[(nc * n_comp as u64) as usize..];
    let mut cache_full = vec![0.0; cache_size_0 as usize];
    let (gctr, cache) = cache_full.split_at_mut((nc * (n_comp as u64)) as usize);
    let mut n: i32 = 0;
    let mut empty: i32 = 1 as i32;
    // if !opt.is_null() {
        // (*envs).opt = opt;
        // n = (((*x_ctr.offset(0 as i32 as isize) == 1 as i32)
        //     as i32) << 3 as i32)
        //     + (((*x_ctr.offset(1 as i32 as isize) == 1 as i32)
        //         as i32) << 2 as i32)
        //     + (((*x_ctr.offset(2 as i32 as isize) == 1 as i32)
        //         as i32) << 1 as i32)
        //     + (*x_ctr.offset(3 as i32 as isize) == 1 as i32)
        //         as i32;
        // (CINTf_2e_loop[n as usize])
        //     .expect("non-null function pointer")(gctr, envs, cache, &mut empty);
    // } else {
    unsafe {
        CINT2e_loop_nopt(gctr.as_mut_ptr(), envs as *mut CINTEnvVars, cache.as_mut_ptr(), &mut empty);
    }
    // }
    let mut counts: [i32; 4] = [0; 4];
    if f_c2s == Some(c2s_sph_2e1_cpy)
        // == ::core::mem::transmute::<
        //     Option::<
        //         unsafe extern "C" fn(
        //             *mut f64,
        //             *mut f64,
        //             *mut i32,
        //             *mut CINTEnvVars,
        //             *mut f64,
        //         ) -> (),
        //     >,
        //     Option::<unsafe extern "C" fn() -> ()>,
        // >(
        //     Some(
        //         c2s_sph_2e1
        //             as unsafe extern "C" fn(
        //                 *mut f64,
        //                 *mut f64,
        //                 *mut i32,
        //                 *mut CINTEnvVars,
        //                 *mut f64,
        //             ) -> (),
        //     ),
        // )
    {
        counts[0] = (envs.i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs.j_l * 2 + 1) * x_ctr[1];
        counts[2] = (envs.k_l * 2 + 1) * x_ctr[2];
        counts[3] = (envs.l_l * 2 + 1) * x_ctr[3];
    } else {
        counts[0] = envs.nfi * x_ctr[0];
        counts[1] = envs.nfj * x_ctr[1];
        counts[2] = envs.c2rust_unnamed.nfk * x_ctr[2];
        counts[3] = envs.c2rust_unnamed_0.nfl* x_ctr[3];
    }
    if dims.len() == 0 {
        for i in 0..4 {
            dims.push(counts[i]);
        }
    }
    let mut nout: i32 = dims[0] * dims[1] * dims[2] * dims[3];
    if empty == 0 {
        n = 0 as i32;
        while n < n_comp {
            f_c2s.expect("non-null")(&mut out[(nout * n) as usize..], &mut gctr[(nc * n as u64) as usize..], &mut dims, &envs, cache);
            // unsafe {
            //     ::core::mem::transmute::<
            //         _,
            //         fn(_, _, _, _, _),
            //     >(
            //         (Some(f_c2s.expect("non-null function pointer")))
            //             .expect("non-null function pointer"),
            //     )(
            //         (&mut out[(nout * n) as usize..]).as_mut_ptr(),
            //         (&mut gctr[(nc * n as u64) as usize..]).as_mut_ptr(),
            //         dims.as_mut_ptr(),
            //         envs as *mut CINTEnvVars,
            //         cache.as_mut_ptr(),
            //     );
            // }
            n += 1;
        }
    } else {
        n = 0 as i32;
        while n < n_comp {
            c2s_dset0_cpy(&mut out[(nout * n) as usize..], &dims, &counts);
            n += 1;
        }
    }
    // if !stack.is_null() {
    //     free(stack as *mut libc::c_void);
    // }
    return (empty == 0) as i32;
}
#[no_mangle]
pub unsafe extern "C" fn CINTgout2e(
    mut gout: *mut f64,
    mut g: *mut f64,
    mut idx: *mut i32,
    mut envs: *mut CINTEnvVars,
    mut gout_empty: i32,
) {
    let mut nf: i32 = (*envs).nf;
    let mut i: i32 = 0;
    let mut ix: i32 = 0;
    let mut iy: i32 = 0;
    let mut iz: i32 = 0;
    let mut n: i32 = 0;
    let mut s: f64 = 0.;
    if gout_empty != 0 {
        match (*envs).nrys_roots {
            1 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            2 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            3 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            4 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize)
                        + *g.offset((ix + 3 as i32) as isize)
                            * *g.offset((iy + 3 as i32) as isize)
                            * *g.offset((iz + 3 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            5 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize)
                        + *g.offset((ix + 3 as i32) as isize)
                            * *g.offset((iy + 3 as i32) as isize)
                            * *g.offset((iz + 3 as i32) as isize)
                        + *g.offset((ix + 4 as i32) as isize)
                            * *g.offset((iy + 4 as i32) as isize)
                            * *g.offset((iz + 4 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            6 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize)
                        + *g.offset((ix + 3 as i32) as isize)
                            * *g.offset((iy + 3 as i32) as isize)
                            * *g.offset((iz + 3 as i32) as isize)
                        + *g.offset((ix + 4 as i32) as isize)
                            * *g.offset((iy + 4 as i32) as isize)
                            * *g.offset((iz + 4 as i32) as isize)
                        + *g.offset((ix + 5 as i32) as isize)
                            * *g.offset((iy + 5 as i32) as isize)
                            * *g.offset((iz + 5 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            7 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize)
                        + *g.offset((ix + 3 as i32) as isize)
                            * *g.offset((iy + 3 as i32) as isize)
                            * *g.offset((iz + 3 as i32) as isize)
                        + *g.offset((ix + 4 as i32) as isize)
                            * *g.offset((iy + 4 as i32) as isize)
                            * *g.offset((iz + 4 as i32) as isize)
                        + *g.offset((ix + 5 as i32) as isize)
                            * *g.offset((iy + 5 as i32) as isize)
                            * *g.offset((iz + 5 as i32) as isize)
                        + *g.offset((ix + 6 as i32) as isize)
                            * *g.offset((iy + 6 as i32) as isize)
                            * *g.offset((iz + 6 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            8 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout
                        .offset(
                            n as isize,
                        ) = *g.offset(ix as isize) * *g.offset(iy as isize)
                        * *g.offset(iz as isize)
                        + *g.offset((ix + 1 as i32) as isize)
                            * *g.offset((iy + 1 as i32) as isize)
                            * *g.offset((iz + 1 as i32) as isize)
                        + *g.offset((ix + 2 as i32) as isize)
                            * *g.offset((iy + 2 as i32) as isize)
                            * *g.offset((iz + 2 as i32) as isize)
                        + *g.offset((ix + 3 as i32) as isize)
                            * *g.offset((iy + 3 as i32) as isize)
                            * *g.offset((iz + 3 as i32) as isize)
                        + *g.offset((ix + 4 as i32) as isize)
                            * *g.offset((iy + 4 as i32) as isize)
                            * *g.offset((iz + 4 as i32) as isize)
                        + *g.offset((ix + 5 as i32) as isize)
                            * *g.offset((iy + 5 as i32) as isize)
                            * *g.offset((iz + 5 as i32) as isize)
                        + *g.offset((ix + 6 as i32) as isize)
                            * *g.offset((iy + 6 as i32) as isize)
                            * *g.offset((iz + 6 as i32) as isize)
                        + *g.offset((ix + 7 as i32) as isize)
                            * *g.offset((iy + 7 as i32) as isize)
                            * *g.offset((iz + 7 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            _ => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    s = 0 as i32 as f64;
                    i = 0 as i32;
                    while i < (*envs).nrys_roots {
                        s
                            += *g.offset((ix + i) as isize)
                                * *g.offset((iy + i) as isize)
                                * *g.offset((iz + i) as isize);
                        i += 1;
                    }
                    *gout.offset(n as isize) = s;
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
        }
    } else {
        match (*envs).nrys_roots {
            1 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            2 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            3 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            4 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize)
                            + *g.offset((ix + 3 as i32) as isize)
                                * *g.offset((iy + 3 as i32) as isize)
                                * *g.offset((iz + 3 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            5 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize)
                            + *g.offset((ix + 3 as i32) as isize)
                                * *g.offset((iy + 3 as i32) as isize)
                                * *g.offset((iz + 3 as i32) as isize)
                            + *g.offset((ix + 4 as i32) as isize)
                                * *g.offset((iy + 4 as i32) as isize)
                                * *g.offset((iz + 4 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            6 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize)
                            + *g.offset((ix + 3 as i32) as isize)
                                * *g.offset((iy + 3 as i32) as isize)
                                * *g.offset((iz + 3 as i32) as isize)
                            + *g.offset((ix + 4 as i32) as isize)
                                * *g.offset((iy + 4 as i32) as isize)
                                * *g.offset((iz + 4 as i32) as isize)
                            + *g.offset((ix + 5 as i32) as isize)
                                * *g.offset((iy + 5 as i32) as isize)
                                * *g.offset((iz + 5 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            7 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize)
                            + *g.offset((ix + 3 as i32) as isize)
                                * *g.offset((iy + 3 as i32) as isize)
                                * *g.offset((iz + 3 as i32) as isize)
                            + *g.offset((ix + 4 as i32) as isize)
                                * *g.offset((iy + 4 as i32) as isize)
                                * *g.offset((iz + 4 as i32) as isize)
                            + *g.offset((ix + 5 as i32) as isize)
                                * *g.offset((iy + 5 as i32) as isize)
                                * *g.offset((iz + 5 as i32) as isize)
                            + *g.offset((ix + 6 as i32) as isize)
                                * *g.offset((iy + 6 as i32) as isize)
                                * *g.offset((iz + 6 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            8 => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    *gout.offset(n as isize)
                        += *g.offset(ix as isize) * *g.offset(iy as isize)
                            * *g.offset(iz as isize)
                            + *g.offset((ix + 1 as i32) as isize)
                                * *g.offset((iy + 1 as i32) as isize)
                                * *g.offset((iz + 1 as i32) as isize)
                            + *g.offset((ix + 2 as i32) as isize)
                                * *g.offset((iy + 2 as i32) as isize)
                                * *g.offset((iz + 2 as i32) as isize)
                            + *g.offset((ix + 3 as i32) as isize)
                                * *g.offset((iy + 3 as i32) as isize)
                                * *g.offset((iz + 3 as i32) as isize)
                            + *g.offset((ix + 4 as i32) as isize)
                                * *g.offset((iy + 4 as i32) as isize)
                                * *g.offset((iz + 4 as i32) as isize)
                            + *g.offset((ix + 5 as i32) as isize)
                                * *g.offset((iy + 5 as i32) as isize)
                                * *g.offset((iz + 5 as i32) as isize)
                            + *g.offset((ix + 6 as i32) as isize)
                                * *g.offset((iy + 6 as i32) as isize)
                                * *g.offset((iz + 6 as i32) as isize)
                            + *g.offset((ix + 7 as i32) as isize)
                                * *g.offset((iy + 7 as i32) as isize)
                                * *g.offset((iz + 7 as i32) as isize);
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
            _ => {
                n = 0 as i32;
                while n < nf {
                    ix = *idx.offset(0 as i32 as isize);
                    iy = *idx.offset(1 as i32 as isize);
                    iz = *idx.offset(2 as i32 as isize);
                    s = 0 as i32 as f64;
                    i = 0 as i32;
                    while i < (*envs).nrys_roots {
                        s
                            += *g.offset((ix + i) as isize)
                                * *g.offset((iy + i) as isize)
                                * *g.offset((iz + i) as isize);
                        i += 1;
                    }
                    *gout.offset(n as isize) += s;
                    n += 1;
                    idx = idx.offset(3 as i32 as isize);
                }
            }
        }
    };
}
#[no_mangle]
pub unsafe fn int2e_sph(
    out: &mut [f64],
    dims: Vec<i32>,
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
    // mut opt: *mut CINTOpt,
    cache: Vec<f64>,
) -> i32 {
    let mut ng: [i32; 8] = [0, 0, 0, 0, 0, 1, 1, 1];
    let mut envs: CINTEnvVars = CINTEnvVars::new();
    CINTinit_int2e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    envs
        .f_gout = ::core::mem::transmute::<
        Option::<
            unsafe extern "C" fn(
                *mut f64,
                *mut f64,
                *mut i32,
                *mut CINTEnvVars,
                i32,
            ) -> (),
        >,
        Option::<unsafe extern "C" fn() -> ()>,
    >(
        Some(
            CINTgout2e
                as unsafe extern "C" fn(
                    *mut f64,
                    *mut f64,
                    *mut i32,
                    *mut CINTEnvVars,
                    i32,
                ) -> (),
        ),
    );
    return CINT2e_drv(
        out,
        dims,
        &mut envs,
        // opt,
        cache,
        Some(c2s_sph_2e1_cpy),
        // ::core::mem::transmute::<
        //     Option::<
        //         unsafe extern "C" fn(
        //             *mut f64,
        //             *mut f64,
        //             *mut i32,
        //             *mut CINTEnvVars,
        //             *mut f64,
        //         ) -> (),
        //     >,
        //     Option::<unsafe extern "C" fn() -> ()>,
        // >(
        //     Some(
        //         c2s_sph_2e1_cpy
        //             as unsafe extern "C" fn(
        //                 *mut f64,
        //                 *mut f64,
        //                 *mut i32,
        //                 *mut CINTEnvVars,
        //                 *mut f64,
        //             ) -> (),
        //     ),
        // ),
    );
}
// #[no_mangle]
// pub unsafe extern "C" fn int2e_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     let mut ng: [i32; 8] = [
//         0 as i32,
//         0 as i32,
//         0 as i32,
//         0 as i32,
//         0 as i32,
//         1 as i32,
//         1 as i32,
//         1 as i32,
//     ];
//     CINTall_2e_optimizer(opt, ng.as_mut_ptr(), atm, natm, bas, nbas, env);
// }
#[no_mangle]
pub unsafe fn int2e_cart(
    out: &mut [f64],
    dims: Vec<i32>,
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
    cache: Vec<f64>,
) -> i32 {
    let ng: [i32; 8] = [0, 0, 0, 0, 0, 1, 1, 1];
    let mut envs: CINTEnvVars = CINTEnvVars::new();
    CINTinit_int2e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    envs
        .f_gout = ::core::mem::transmute::<
        Option::<
            unsafe extern "C" fn(
                *mut f64,
                *mut f64,
                *mut i32,
                *mut CINTEnvVars,
                i32,
            ) -> (),
        >,
        Option::<unsafe extern "C" fn() -> ()>,
    >(
        Some(
            CINTgout2e
                as unsafe extern "C" fn(
                    *mut f64,
                    *mut f64,
                    *mut i32,
                    *mut CINTEnvVars,
                    i32,
                ) -> (),
        ),
    );
    return CINT2e_drv(
        out,
        dims,
        &mut envs,
        cache,
        Some(c2s_cart_2e1_cpy),
        // ::core::mem::transmute::<
        //     Option::<
        //         unsafe extern "C" fn(
        //             *mut f64,
        //             *mut f64,
        //             *mut i32,
        //             *mut CINTEnvVars,
        //             *mut f64,
        //         ) -> (),
        //     >,
        //     Option::<unsafe extern "C" fn() -> ()>,
        // >(
        //     Some(
        //         c2s_cart_2e1
        //             as unsafe extern "C" fn(
        //                 *mut f64,
        //                 *mut f64,
        //                 *mut i32,
        //                 *mut CINTEnvVars,
        //                 *mut f64,
        //             ) -> (),
        //     ),
        // ),
    );
}
#[no_mangle]
pub fn cint2e_sph(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
    // mut opt: *mut CINTOpt,
) -> i32 {
    let dims = vec![0;0];
    let cache = vec![0.0;0];
    unsafe {
        return int2e_sph(
            out,
            dims,
            shls,
            atm,
            natm,
            bas,
            nbas,
            env,
            cache,
        );
    }
}
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int2e_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_sph_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int2e_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_cart_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int2e_optimizer(opt, atm, natm, bas, nbas, env);
// }
#[no_mangle]
pub fn cint2e_cart(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32 {
    let dims = vec![0;0];
    let cache = vec![0.0;0];
    unsafe {
        return int2e_cart(
            out,
            dims,
            shls,
            atm,
            natm,
            bas,
            nbas,
            env,
            cache,
        );
    }
}
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_sph_(
//     mut out: *mut f64,
//     mut shls: *mut i32,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
//     mut optptr_as_integer8: size_t,
// ) -> i32 {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     todo!();
//     // return int2e_sph(
//     //     out,
//     //     0 as *mut i32,
//     //     shls,
//     //     atm,
//     //     *natm,
//     //     bas,
//     //     *nbas,
//     //     env,
//     //     *opt,
//     //     0 as *mut f64,
//     // );
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_sph_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int2e_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_cart_(
//     mut out: *mut f64,
//     mut shls: *mut i32,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
//     mut optptr_as_integer8: size_t,
// ) -> i32 {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     todo!()
//     // return int2e_cart(
//     //     out,
//     //     0 as *mut i32,
//     //     shls,
//     //     atm,
//     //     *natm,
//     //     bas,
//     //     *nbas,
//     //     env,
//     //     *opt,
//     //     0 as *mut f64,
//     // );
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_cart_optimizer_(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     int2e_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint2e_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int2e_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
