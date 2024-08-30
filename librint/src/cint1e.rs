#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments
)]

use crate::fblas::CINTdmat_transpose;
use crate::g1e::CINTcommon_fac_sp;
use crate::g1e::CINTg1e_index_xyz;
use crate::g1e::CINTg1e_nuc;
use crate::g1e::CINTg1e_ovlp;
use crate::g1e::CINTinit_int1e_EnvVars;
use crate::g1e::CINTprim_to_ctr_0;
use crate::g1e::CINTprim_to_ctr_1;
use crate::optimizer::CINTOpt_log_max_pgto_coeff;
use crate::optimizer::CINTOpt_non0coeff_byshell;
use crate::optimizer::CINTset_pairdata;

use crate::fblas::CINTdmat_transpose_cpy;

use crate::g1e::CINTg1e_nuc_cpy;
use crate::g1e::CINTg1e_ovlp_cpy;

use crate::g1e::CINTg1e_index_xyz_cpy;
use crate::g1e::CINTprim_to_ctr_0_cpy;
use crate::g1e::CINTprim_to_ctr_1_cpy;
use crate::optimizer::CINTOpt_non0coeff_byshell_cpy;

use crate::optimizer::CINTOpt_log_max_pgto_coeff_cpy;
use crate::optimizer::CINTset_pairdata_cpy;

use crate::cart2sph::c2s_cart_1e_cpy;
use crate::cart2sph::c2s_dset0_cpy;
use crate::cart2sph::c2s_sph_1e_cpy;

use crate::cint::PairData;
use crate::cint::CINTEnvVars;
use crate::cint::F_FC2S;

pub type size_t = libc::c_ulong;
pub type uintptr_t = libc::c_ulong;

extern "C" {
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(__ptr: *mut libc::c_void);
}


#[no_mangle]
pub fn CINT1e_loop_cpy(
    gctr: &mut [f64],
    envs: &mut CINTEnvVars,
    _cache: &mut [f64],
    int1e_type: i32,
) -> i32 {
    let shls: [i32; 4] = envs.shls;
    let bas: &[i32] = &envs.bas;
    let env: &[f64] = &envs.env;
    let i_sh: usize = shls[0] as usize;
    let j_sh: usize = shls[1] as usize;
    let i_ctr: usize = envs.x_ctr[0] as usize;
    let j_ctr: usize = envs.x_ctr[1] as usize;
    let i_prim: usize = bas[8 * i_sh + 2] as usize;
    let j_prim: usize = bas[8 * j_sh + 2] as usize;

    let ai: &[f64] = &env[(bas[8 * i_sh + 5] as usize)..(bas[8 * i_sh + 5] as usize + i_prim)];
    let aj: &[f64] = &env[(bas[8 * j_sh + 5] as usize)..(bas[8 * j_sh + 5] as usize + j_prim)];
    let ci: &[f64] = &env[(bas[8 * i_sh + 6] as usize)..(bas[8 * i_sh + 6] as usize + i_prim * i_ctr)];
    let cj: &[f64] = &env[(bas[8 * j_sh + 6] as usize)..(bas[8 * j_sh + 6] as usize + j_prim * j_ctr)];

    let n_comp: i32 = envs.ncomp_e1 * envs.ncomp_tensor;
    let expcutoff: f64 = envs.expcutoff;

    let mut log_maxci = vec![0.0; i_prim as usize].into_boxed_slice();
    let mut log_maxcj = vec![0.0; j_prim as usize].into_boxed_slice();

    CINTOpt_log_max_pgto_coeff_cpy(&mut log_maxci, ci, i_prim as i32, i_ctr as i32);
    CINTOpt_log_max_pgto_coeff_cpy(&mut log_maxcj, cj, j_prim as i32, j_ctr as i32);

    let mut pdata_base = vec![PairData::new(); i_prim * j_prim].into_boxed_slice();

    if CINTset_pairdata_cpy(
        &mut pdata_base,
        ai,
        aj,
        &envs.ri,
        &envs.rj,
        &log_maxci,
        &log_maxcj,
        envs.li_ceil,
        envs.lj_ceil,
        i_prim,
        j_prim,
        envs.rirj[0] * envs.rirj[0] + envs.rirj[1] * envs.rirj[1] + envs.rirj[2] * envs.rirj[2],
        expcutoff,
        env,
    ) != 0
    {
        return 0 as i32;
    }

    let mut rij: [f64; 3];

    let mut idx: Box<[i32]> = vec![0; envs.nf as usize* 3].into_boxed_slice();

    CINTg1e_index_xyz_cpy(&mut idx, envs);

    let mut non0ctri = vec![0; i_prim].into_boxed_slice();
    let mut non0ctrj = vec![0; j_prim].into_boxed_slice();
    let mut non0idxi = vec![0; i_prim * i_ctr].into_boxed_slice();
    let mut non0idxj = vec![0; j_prim * j_ctr].into_boxed_slice();

    CINTOpt_non0coeff_byshell_cpy(&mut non0idxi, &mut non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell_cpy(&mut non0idxj, &mut non0ctrj, cj, j_prim, j_ctr);

    let nc: i32 = (i_ctr * j_ctr) as i32;
    let leng: i32 = envs.g_size * 3 * ((1 << envs.gbits) + 1);
    let lenj: i32 = envs.nf * nc * n_comp;
    let leni: i32 = envs.nf * i_ctr as i32 * n_comp;
    let len0: i32 = envs.nf * n_comp;
    // let len: i32 = leng + lenj + leni + len0;

    let mut empty: [i32; 4] = [1, 1, 1, 1];
    let mut gempty_idx = 0;
    let mut iempty_idx = 1;
    let jempty_idx = 2;

    let mut g = vec![0.0; leng as usize].into_boxed_slice();
    // let mut g1 = vec![0.0; (lenj + leni + len0) as usize];

    let mut gctrj = vec![0.0; lenj as usize];
    let mut gctri = vec![0.0; leni as usize];
    let mut gout = vec![0.0; len0 as usize];

    // if n_comp == 1 {
    //     gctrj = gctr;
    // } else {
    //     (gctrj, g1) = g1.split_at_mut(lenj as usize);
    // }

    if j_ctr == 1 {
        // gctri = gctrj;
        iempty_idx = jempty_idx;
    } else {
        // (gctri, g1) = g1.split_at_mut(leni as usize);
    }

    if i_ctr == 1 {
        // gout = gctri;
        gempty_idx = iempty_idx;
    } else {
        // gout = g1;
    }

    let common_factor: f64 = envs.common_factor * CINTcommon_fac_sp(envs.i_l) * CINTcommon_fac_sp(envs.j_l);

    let pdata_ij = pdata_base;

    let mut fac1i: f64 = 0.;
    let mut fac1j: f64 = 0.;
    let mut expij: f64 = 0.;

    let mut offset = 0;

    for jp in 0..j_prim {
        envs.aj[0] = aj[jp];
        if j_ctr == 1 {
            fac1j = common_factor * cj[jp];
        } else {
            fac1j = common_factor;
            empty[iempty_idx] = 1;
        }

        for ip in 0..i_prim {
            if !(pdata_ij[offset].cceij > expcutoff) {
                envs.ai[0] = ai[ip];
                expij = pdata_ij[offset].eij;
                rij = pdata_ij[offset].rij;
                envs.rij[0] = rij[0];
                envs.rij[1] = rij[1];
                envs.rij[2] = rij[2];
                if i_ctr == 1 {
                    fac1i = fac1j * ci[ip] * expij;
                } else {
                    fac1i = fac1j * expij;
                }
                envs.fac[0] = fac1i;

                make_g1e_gout_cpy(&mut gout, &mut g, &idx, envs, empty[gempty_idx], int1e_type);

                if i_ctr > 1 {
                    if empty[iempty_idx] != 0 {
                        CINTprim_to_ctr_0_cpy(
                            &mut gctri,
                            &gout,
                            &ci[ip..],
                            (envs.nf * n_comp) as usize,
                            i_prim,
                            i_ctr,
                            non0ctri[ip] as usize,
                            &non0idxi[(ip * i_ctr)..],
                        );
                    } else {
                        CINTprim_to_ctr_1_cpy(
                            &mut gctri,
                            &gout,
                            &ci[ip..],
                            (envs.nf * n_comp) as usize,
                            i_prim,
                            i_ctr,
                            non0ctri[ip] as usize,
                            &non0idxi[(ip * i_ctr)..],
                        );
                    }
                } else {
                    for i in 0..len0 as usize {
                        gctri[i] = gout[i];
                    }
                }
                empty[iempty_idx] = 0;
            }
            offset += 1;
        }
        if empty[iempty_idx] == 0 {
            if j_ctr > 1 {
                if empty[jempty_idx] != 0 {
                    CINTprim_to_ctr_0_cpy(
                        &mut gctrj,
                        &gctri,
                        &cj[jp..],
                        (envs.nf * i_ctr as i32 * n_comp) as usize,
                        j_prim,
                        j_ctr,
                        non0ctrj[jp] as usize,
                        &non0idxj[(jp * j_ctr)..],
                    );
                } else {
                    CINTprim_to_ctr_1_cpy(
                        &mut gctrj,
                        &gctri,
                        &cj[jp..],
                        (envs.nf * i_ctr as i32 * n_comp) as usize,
                        j_prim,
                        j_ctr,
                        non0ctrj[jp] as usize,
                        &non0idxj[(jp * j_ctr)..],
                    );
                }
            } else {
                for i in 0..leni as usize { 
                    gctrj[i] = gctri[i];
                }
            }
            empty[jempty_idx] = 0;
        }
    }

    if n_comp > 1 && empty[jempty_idx] == 0 {
        CINTdmat_transpose_cpy(gctr, &gctrj, (envs.nf * nc) as usize, n_comp as usize);
    } else {
        for i in 0..lenj as usize {
            gctr[i] = gctrj[i];
        }
    }
    return (empty[jempty_idx] == 0) as i32;
}

#[no_mangle]
pub unsafe extern "C" fn CINT1e_loop(
    mut gctr: *mut f64,
    mut envs: *mut CINTEnvVars,
    mut cache: *mut f64,
    mut int1e_type: i32,
) -> i32 {
    let mut shls: *mut i32 = (*envs).shls.as_mut_ptr();
    let mut bas: *mut i32 = (*envs).bas.as_mut_ptr();
    let mut env: *mut f64 = (*envs).env.as_mut_ptr();
    let mut i_sh: i32 = *shls.offset(0 as i32 as isize);
    let mut j_sh: i32 = *shls.offset(1 as i32 as isize);
    let mut i_ctr: i32 = (*envs).x_ctr[0 as i32 as usize];
    let mut j_ctr: i32 = (*envs).x_ctr[1 as i32 as usize];
    let mut i_prim: i32 = *bas.offset((8 as i32 * i_sh + 2 as i32) as isize);
    let mut j_prim: i32 = *bas.offset((8 as i32 * j_sh + 2 as i32) as isize);
    let mut ai: *mut f64 = env.offset(*bas.offset((8 as i32 * i_sh + 5 as i32) as isize) as isize);
    let mut aj: *mut f64 = env.offset(*bas.offset((8 as i32 * j_sh + 5 as i32) as isize) as isize);
    let mut ci: *mut f64 = env.offset(*bas.offset((8 as i32 * i_sh + 6 as i32) as isize) as isize);
    let mut cj: *mut f64 = env.offset(*bas.offset((8 as i32 * j_sh + 6 as i32) as isize) as isize);
    let mut n_comp: i32 = (*envs).ncomp_e1 * (*envs).ncomp_tensor; // cannot handle (augmented) unknown intrinsic
    let mut expcutoff: f64 = (*envs).expcutoff;
    let mut log_maxci: *mut f64 = 0 as *mut f64;
    let mut log_maxcj: *mut f64 = 0 as *mut f64;
    let mut pdata_base: *mut PairData = 0 as *mut PairData;
    let mut pdata_ij: *mut PairData = 0 as *mut PairData;
    log_maxci = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void as *mut f64;
    cache = log_maxci.offset((i_prim + j_prim) as isize);
    pdata_base = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void
        as *mut PairData;
    cache = pdata_base.offset((i_prim * j_prim) as isize) as *mut f64;
    log_maxcj = log_maxci.offset(i_prim as isize);
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
        (*envs).rirj[0 as i32 as usize] * (*envs).rirj[0 as i32 as usize]
            + (*envs).rirj[1 as i32 as usize] * (*envs).rirj[1 as i32 as usize]
            + (*envs).rirj[2 as i32 as usize] * (*envs).rirj[2 as i32 as usize],
        expcutoff,
        env,
    ) != 0
    {
        return 0 as i32;
    }
    let mut fac1i: f64 = 0.;
    let mut fac1j: f64 = 0.;
    let mut expij: f64 = 0.;
    let mut ip: i32 = 0;
    let mut jp: i32 = 0;
    let mut empty: [i32; 4] = [1 as i32, 1 as i32, 1 as i32, 1 as i32];
    let mut gempty: *mut i32 = empty.as_mut_ptr().offset(0 as i32 as isize);
    let mut iempty: *mut i32 = empty.as_mut_ptr().offset(1 as i32 as isize);
    let mut jempty: *mut i32 = empty.as_mut_ptr().offset(2 as i32 as isize);
    let mut rij: *mut f64 = 0 as *mut f64;
    let mut idx: *mut i32 = 0 as *mut i32;
    idx = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void as *mut i32;
    cache = idx.offset(((*envs).nf * 3 as i32) as isize) as *mut f64;
    CINTg1e_index_xyz(idx, envs);
    let mut non0ctri: *mut i32 = 0 as *mut i32;
    let mut non0ctrj: *mut i32 = 0 as *mut i32;
    let mut non0idxi: *mut i32 = 0 as *mut i32;
    let mut non0idxj: *mut i32 = 0 as *mut i32;
    non0ctri = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void as *mut i32;
    cache =
        non0ctri.offset((i_prim + j_prim + i_prim * i_ctr + j_prim * j_ctr) as isize) as *mut f64;
    non0ctrj = non0ctri.offset(i_prim as isize);
    non0idxi = non0ctrj.offset(j_prim as isize);
    non0idxj = non0idxi.offset((i_prim * i_ctr) as isize);
    CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
    CINTOpt_non0coeff_byshell(non0idxj, non0ctrj, cj, j_prim, j_ctr);
    let nc: i32 = i_ctr * j_ctr;
    let leng: i32 = (*envs).g_size * 3 as i32 * (((1 as i32) << (*envs).gbits) + 1 as i32);
    let lenj: i32 = (*envs).nf * nc * n_comp;
    let leni: i32 = (*envs).nf * i_ctr * n_comp;
    let len0: i32 = (*envs).nf * n_comp;
    let len: i32 = leng + lenj + leni + len0;
    let mut g: *mut f64 = 0 as *mut f64;
    let mut gout: *mut f64 = 0 as *mut f64;
    let mut gctri: *mut f64 = 0 as *mut f64;
    let mut gctrj: *mut f64 = 0 as *mut f64;
    g = ((cache as uintptr_t).wrapping_add(7 as i32 as libc::c_ulong)
        & (8 as i32 as uintptr_t).wrapping_neg()) as *mut libc::c_void as *mut f64;
    cache = g.offset(len as isize);
    let mut g1: *mut f64 = g.offset(leng as isize);
    if n_comp == 1 as i32 {
        gctrj = gctr;
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
        gout = g1; // g1[leni + lenj + len0]
    }

    let mut common_factor: f64 =
        (*envs).common_factor * CINTcommon_fac_sp((*envs).i_l) * CINTcommon_fac_sp((*envs).j_l);
    pdata_ij = pdata_base;
    jp = 0 as i32;
    while jp < j_prim {
        (*envs).aj[0 as i32 as usize] = *aj.offset(jp as isize);
        if j_ctr == 1 as i32 {
            fac1j = common_factor * *cj.offset(jp as isize);
        } else {
            fac1j = common_factor;
            *iempty = 1 as i32;
        }
        ip = 0 as i32;
        while ip < i_prim {
            if !((*pdata_ij).cceij > expcutoff) {
                (*envs).ai[0 as i32 as usize] = *ai.offset(ip as isize);
                expij = (*pdata_ij).eij;
                rij = ((*pdata_ij).rij).as_mut_ptr();
                (*envs).rij[0 as i32 as usize] = *rij.offset(0 as i32 as isize);
                (*envs).rij[1 as i32 as usize] = *rij.offset(1 as i32 as isize);
                (*envs).rij[2 as i32 as usize] = *rij.offset(2 as i32 as isize);
                if i_ctr == 1 as i32 {
                    fac1i = fac1j * *ci.offset(ip as isize) * expij;
                } else {
                    fac1i = fac1j * expij;
                }
                (*envs).fac[0 as i32 as usize] = fac1i;
                make_g1e_gout(gout, g, idx, envs, *gempty, int1e_type);
                if i_ctr > 1 as i32 {
                    if *iempty != 0 {
                        CINTprim_to_ctr_0(
                            gctri,
                            gout,
                            ci.offset(ip as isize),
                            ((*envs).nf * n_comp) as size_t,
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
                            ((*envs).nf * n_comp) as size_t,
                            i_prim,
                            i_ctr,
                            *non0ctri.offset(ip as isize),
                            non0idxi.offset((ip * i_ctr) as isize),
                        );
                    }
                }
                *iempty = 0 as i32;
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
                        ((*envs).nf * i_ctr * n_comp) as size_t,
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
                        ((*envs).nf * i_ctr * n_comp) as size_t,
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
    if n_comp > 1 as i32 && *jempty == 0 {
        CINTdmat_transpose(gctr, gctrj, (*envs).nf * nc, n_comp);
    }
    return (*jempty == 0) as i32;
}
#[no_mangle]
pub fn int1e_cache_size(envs: &CINTEnvVars) -> i32 {
    let shls: [i32; 4] = envs.shls;
    let bas: &[i32] = &envs.bas;
    let i_prim: i32 = bas[8 * shls[0] as usize + 2];
    let j_prim: i32 = bas[8 * shls[1] as usize + 2];
    let x_ctr: [i32; 4] = envs.x_ctr;
    let nc: i32 = envs.nf * x_ctr[0] * x_ctr[1];
    let n_comp: i32 = envs.ncomp_e1 * envs.ncomp_tensor;
    let leng: i32 = envs.g_size * 3 * (((1 as i32) << envs.gbits) + 1);
    let lenj: i32 = envs.nf * nc * n_comp;
    let leni: i32 = envs.nf * x_ctr[0] * n_comp;
    let len0: i32 = envs.nf * n_comp;
    let pdata_size: i32 = i_prim * j_prim * 5
        + i_prim * x_ctr[0]
        + j_prim * x_ctr[1]
        + (i_prim + j_prim) * 2
        + envs.nf * 3;
    let cache_size: i32 =
        if nc * n_comp + leng + lenj + leni + len0 + pdata_size > nc * n_comp + envs.nf * 8 * 2 {
            nc * n_comp + leng + lenj + leni + len0 + pdata_size
        } else {
            nc * n_comp + envs.nf * 8 * 2
        };
    return cache_size;
}
#[no_mangle]
pub fn CINT1e_drv(
    out: &mut [f64],
    mut dims: Vec<i32>,
    envs: &mut CINTEnvVars,
    _cache: Vec<f64>,
    f_c2s: Option<F_FC2S>,
    int1e_type: i32,
) -> i32 {
    let x_ctr: [i32; 4] = envs.x_ctr;
    let nc: i32 = envs.nf * x_ctr[0] * x_ctr[1];
    let n_comp: i32 = envs.ncomp_e1 * envs.ncomp_tensor;

    let cache_size = int1e_cache_size(envs) as usize;

    let mut cache_full = vec![0.0; cache_size as usize];
    let (_gctr, cache) = cache_full.split_at_mut((nc * n_comp) as usize);

    let mut gctr = vec![0.0; (nc * n_comp) as usize];
    let has_value: i32;

    has_value = CINT1e_loop_cpy(&mut gctr, envs, cache, int1e_type);
    // unsafe {
    //     has_value = CINT1e_loop(gctr.as_mut_ptr(), envs as *mut CINTEnvVars, cache.as_mut_ptr(), int1e_type);
    // }

    let mut counts: [i32; 4] = [0; 4];

    if f_c2s == Some(c2s_sph_1e_cpy) {
        counts[0] = (envs.i_l * 2 + 1) * x_ctr[0];
        counts[1] = (envs.j_l * 2 + 1) * x_ctr[1];
    } else if f_c2s == Some(c2s_cart_1e_cpy) {
        counts[0] = envs.nfi * x_ctr[0];
        counts[1] = envs.nfj * x_ctr[1];
    }
    counts[2] = 1;
    counts[3] = 1;

    if dims.len() == 0 {
        for i in 0..4 {
            dims.push(counts[i]);
        }
    }

    let nout: i32 = dims[0] * dims[1];
    let mut n: i32 = 0;
    if has_value != 0 {
        n = 0 as i32;
        while n < n_comp {
            f_c2s.expect("non-null")(
                &mut out[(nout * n) as usize..],
                &mut gctr[(nc * n) as usize..],
                &mut dims,
                &envs,
                cache,
            );
            n += 1;
        }
    } else {
        n = 0 as i32;
        while n < n_comp {
            c2s_dset0_cpy(&mut out[(nout * n) as usize..], &dims, &counts);
            n += 1;
        }
    }
    return has_value;
}
fn make_g1e_gout_cpy(
    mut gout: &mut [f64],
    g: &mut [f64],
    idx: &[i32],
    envs: &CINTEnvVars,
    empty: i32,
    int1e_type: i32,
) {
    unsafe {
    let mut ia: i32 = 0;
    match int1e_type {
        0 => {
            CINTg1e_ovlp_cpy(g, envs);
            // ::core::mem::transmute::<
            //     _,
            //     fn(_, _, _, _, _),
            // >(
            //     (Some(((*envs).f_gout).expect("non-null function pointer")))
            //         .expect("non-null function pointer"),
            // )(gout, g, idx, envs, empty);
            CINTgout1e_cpy(gout, g, idx, envs, empty);
        }
        // 1 => {
        //     CINTg1e_nuc(g.as_mut_ptr(), envs as *mut CINTEnvVars, -(1 as i32));
        //     ::core::mem::transmute::<
        //         _,
        //         fn(_, _, _, _, _),
        //     >(
        //         (Some(((*envs).f_gout).expect("non-null function pointer")))
        //             .expect("non-null function pointer"),
        //     )(gout, g, idx, envs, empty);
        // }
        // 2 => {
        //     ia = 0 as i32;
        //     while ia < (*envs).natm {
        //         CINTg1e_nuc(g.as_mut_ptr(), envs as *mut CINTEnvVars, ia);
        //         ::core::mem::transmute::<
        //             _,
        //             fn(_, _, _, _, _),
        //         >(
        //             (Some(((*envs).f_gout).expect("non-null function pointer")))
        //                 .expect("non-null function pointer"),
        //         )(
        //             &mut gout,
        //             &g,
        //             idx,
        //             envs,
        //             (empty != 0 && ia == 0 as i32) as i32,
        //         );
        //         ia += 1;
        //     }
        // }
        _ => {}
    };
}
}
unsafe extern "C" fn make_g1e_gout(
    mut gout: *mut f64,
    mut g: *mut f64,
    mut idx: *mut i32,
    mut envs: *mut CINTEnvVars,
    mut empty: i32,
    mut int1e_type: i32,
) {
    let mut ia: i32 = 0;
    match int1e_type {
        0 => {
            CINTg1e_ovlp(g, envs);
            ::core::mem::transmute::<_, fn(_, _, _, _, _)>(
                (Some(((*envs).f_gout).expect("non-null function pointer")))
                    .expect("non-null function pointer"),
            )(gout, g, idx, envs, empty);
        }
        1 => {
            CINTg1e_nuc(g, envs, -(1 as i32));
            ::core::mem::transmute::<_, fn(_, _, _, _, _)>(
                (Some(((*envs).f_gout).expect("non-null function pointer")))
                    .expect("non-null function pointer"),
            )(gout, g, idx, envs, empty);
        }
        2 => {
            ia = 0 as i32;
            while ia < (*envs).natm {
                CINTg1e_nuc(g, envs, ia);
                ::core::mem::transmute::<_, fn(_, _, _, _, _)>(
                    (Some(((*envs).f_gout).expect("non-null function pointer")))
                        .expect("non-null function pointer"),
                )(gout, g, idx, envs, (empty != 0 && ia == 0 as i32) as i32);
                ia += 1;
            }
        }
        _ => {}
    };
}

#[no_mangle]
pub fn CINTgout1e_cpy(gout: &mut [f64], g: &[f64], idx: &[i32], envs: &CINTEnvVars, empty: i32) {
    let nf: usize = envs.nf as usize;
    let mut ix: usize = 0;
    let mut iy: usize = 0;
    let mut iz: usize = 0;
    if empty != 0 {
        for n in 0..nf {
            ix = idx[n * 3 + 0] as usize;
            iy = idx[n * 3 + 1] as usize;
            iz = idx[n * 3 + 2] as usize;
            gout[n] = g[ix] * g[iy] * g[iz];
        }
    } else {
        for n in 0..nf {
            ix = idx[n * 3 + 0] as usize;
            iy = idx[n * 3 + 1] as usize;
            iz = idx[n * 3 + 2] as usize;
            gout[n] += g[ix] * g[iy] * g[iz];
        }
    };
}
#[no_mangle]
pub fn CINTgout1e_nuc_cpy(
    gout: &mut [f64],
    g: &[f64],
    idx: &[i32],
    envs: &CINTEnvVars,
    empty: i32,
) {
    let nf: usize = envs.nf as usize;
    let nrys_roots: usize = envs.nrys_roots as usize;
    let mut gx: &[f64];
    let mut gy: &[f64];
    let mut gz: &[f64];
    let mut s: f64 = 0.;
    if empty != 0 {
        for n in 0..nf {
            gx = &g[idx[n * 3 + 0] as usize..(idx[n * 3 + 0] as usize) + nrys_roots];
            gy = &g[idx[n * 3 + 1] as usize..(idx[n * 3 + 1] as usize) + nrys_roots];
            gz = &g[idx[n * 3 + 2] as usize..(idx[n * 3 + 2] as usize) + nrys_roots];

            s = 0.0;
            for i in 0..nrys_roots {
                s += gx[i] * gy[i] * gz[i];
            }
            gout[n] = s;
        }
    } else {
        for n in 0..nf {
            gx = &g[idx[n * 3 + 0] as usize..(idx[n * 3 + 0] as usize) + nrys_roots];
            gy = &g[idx[n * 3 + 1] as usize..(idx[n * 3 + 1] as usize) + nrys_roots];
            gz = &g[idx[n * 3 + 2] as usize..(idx[n * 3 + 2] as usize) + nrys_roots];

            s = 0.0;
            for i in 0..nrys_roots {
                s += gx[i] * gy[i] * gz[i];
            }
            gout[n] += s;
        }
    };
}

#[no_mangle]
pub unsafe extern "C" fn CINTgout1e(
    mut gout: *mut f64,
    mut g: *mut f64,
    mut idx: *mut i32,
    mut envs: *mut CINTEnvVars,
    mut empty: i32,
) {
    let mut nf: i32 = (*envs).nf;
    let mut n: i32 = 0;
    let mut ix: i32 = 0;
    let mut iy: i32 = 0;
    let mut iz: i32 = 0;
    if empty != 0 {
        n = 0 as i32;
        while n < nf {
            ix = *idx.offset((n * 3 as i32 + 0 as i32) as isize);
            iy = *idx.offset((n * 3 as i32 + 1 as i32) as isize);
            iz = *idx.offset((n * 3 as i32 + 2 as i32) as isize);
            *gout.offset(n as isize) =
                *g.offset(ix as isize) * *g.offset(iy as isize) * *g.offset(iz as isize);
            n += 1;
        }
    } else {
        n = 0 as i32;
        while n < nf {
            ix = *idx.offset((n * 3 as i32 + 0 as i32) as isize);
            iy = *idx.offset((n * 3 as i32 + 1 as i32) as isize);
            iz = *idx.offset((n * 3 as i32 + 2 as i32) as isize);
            *gout.offset(n as isize) +=
                *g.offset(ix as isize) * *g.offset(iy as isize) * *g.offset(iz as isize);
            n += 1;
        }
    };
}
#[no_mangle]
pub unsafe extern "C" fn CINTgout1e_nuc(
    mut gout: *mut f64,
    mut g: *mut f64,
    mut idx: *mut i32,
    mut envs: *mut CINTEnvVars,
    mut empty: i32,
) {
    let mut nf: i32 = (*envs).nf;
    let mut nrys_roots: i32 = (*envs).nrys_roots;
    let mut n: i32 = 0;
    let mut i: i32 = 0;
    let mut gx: *mut f64 = 0 as *mut f64;
    let mut gy: *mut f64 = 0 as *mut f64;
    let mut gz: *mut f64 = 0 as *mut f64;
    let mut s: f64 = 0.;
    if empty != 0 {
        n = 0 as i32;
        while n < nf {
            gx = g.offset(*idx.offset((n * 3 as i32 + 0 as i32) as isize) as isize);
            gy = g.offset(*idx.offset((n * 3 as i32 + 1 as i32) as isize) as isize);
            gz = g.offset(*idx.offset((n * 3 as i32 + 2 as i32) as isize) as isize);
            s = 0 as i32 as f64;
            i = 0 as i32;
            while i < nrys_roots {
                s += *gx.offset(i as isize) * *gy.offset(i as isize) * *gz.offset(i as isize);
                i += 1;
            }
            *gout.offset(n as isize) = s;
            n += 1;
        }
    } else {
        n = 0 as i32;
        while n < nf {
            gx = g.offset(*idx.offset((n * 3 as i32 + 0 as i32) as isize) as isize);
            gy = g.offset(*idx.offset((n * 3 as i32 + 1 as i32) as isize) as isize);
            gz = g.offset(*idx.offset((n * 3 as i32 + 2 as i32) as isize) as isize);
            s = 0 as i32 as f64;
            i = 0 as i32;
            while i < nrys_roots {
                s += *gx.offset(i as isize) * *gy.offset(i as isize) * *gz.offset(i as isize);
                i += 1;
            }
            *gout.offset(n as isize) += s;
            n += 1;
        }
    };
}
#[no_mangle]
pub unsafe fn int1e_ovlp_sph(
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
    let ng = [0, 0, 0, 0, 0, 1, 1, 1];
    let mut envs = CINTEnvVars::new();
    CINTinit_int1e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    // envs.f_gout = ::core::mem::transmute::<
    //     Option<unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> ()>,
    //     Option<unsafe extern "C" fn() -> ()>,
    // >(Some(
    //     CINTgout1e
    //         as unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> (),
    // ));
    envs.f_gout = Some(CINTgout1e_cpy);
    return CINT1e_drv(
        out,
        dims,
        &mut envs,
        cache,
        Some(c2s_sph_1e_cpy),
        0 as i32,
    );
}
#[no_mangle]
pub unsafe fn int1e_ovlp_cart(
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
    let ng = [0, 0, 0, 0, 0, 1, 1, 1];
    let mut envs = CINTEnvVars::new();
    CINTinit_int1e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    // envs.f_gout = ::core::mem::transmute::<
    //     Option<unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> ()>,
    //     Option<unsafe extern "C" fn() -> ()>,
    // >(Some(
    //     CINTgout1e
    //         as unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> (),
    // ));
    envs.f_gout = Some(CINTgout1e_cpy);
    return CINT1e_drv(
        out,
        dims,
        &mut envs,
        cache,
        Some(c2s_cart_1e_cpy),
        0 as i32,
    );
}
#[no_mangle]
pub unsafe fn int1e_nuc_sph(
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
    let ng = [0, 0, 0, 0, 0, 1, 0, 1];
    let mut envs = CINTEnvVars::new();
    CINTinit_int1e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    // envs.f_gout = ::core::mem::transmute::<
    //     Option<unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> ()>,
    //     Option<unsafe extern "C" fn() -> ()>,
    // >(Some(
    //     CINTgout1e_nuc
    //         as unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> (),
    // ));
    envs.f_gout = Some(CINTgout1e_nuc_cpy);
    return CINT1e_drv(
        out,
        dims,
        &mut envs,
        cache,
        Some(c2s_sph_1e_cpy),
        2 as i32,
    );
}
#[no_mangle]
pub unsafe fn int1e_nuc_cart(
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
    let ng = [0, 0, 0, 0, 0, 1, 0, 1];
    let mut envs = CINTEnvVars::new();
    CINTinit_int1e_EnvVars(&mut envs, &ng, shls, atm, natm, bas, nbas, env);
    // envs.f_gout = ::core::mem::transmute::<
    //     Option<unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> ()>,
    //     Option<unsafe extern "C" fn() -> ()>,
    // >(Some(
    //     CINTgout1e_nuc
    //         as unsafe extern "C" fn(*mut f64, *mut f64, *mut i32, *mut CINTEnvVars, i32) -> (),
    // ));
    envs.f_gout = Some(CINTgout1e_nuc_cpy);
    return CINT1e_drv(
        out,
        dims,
        &mut envs,
        cache,
        Some(c2s_cart_1e_cpy),
        2 as i32,
    );
}
#[no_mangle]
pub fn cint1e_ovlp_cart(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32 {
    let dims = vec![0; 0];
    let cache = vec![0.0; 0];
    unsafe {
        return int1e_ovlp_cart(out, dims, shls, atm, natm, bas, nbas, env, cache);
    }
}
#[no_mangle]
pub fn cint1e_ovlp_sph(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32 {
    let dims = vec![0; 0];
    let cache = vec![0.0; 0];
    unsafe {
        return int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, cache);
    }
}
#[no_mangle]
pub fn cint1e_nuc_cart(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32 {
    let dims = vec![0; 0];
    let cache = vec![0.0; 0];
    unsafe {
        return int1e_nuc_cart(out, dims, shls, atm, natm, bas, nbas, env, cache);
    }
}

#[no_mangle]
pub fn cint1e_nuc_sph(
    out: &mut [f64],
    shls: [i32; 4],
    atm: &mut [i32],
    natm: i32,
    bas: &mut [i32],
    nbas: i32,
    env: &mut [f64],
) -> i32 {
    let dims = vec![0; 0];
    let cache = vec![0.0; 0];
    unsafe {
        return int1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, cache);
    }
}
// this function does nothing
// #[no_mangle]
// pub unsafe extern "C" fn int1e_ovlp_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     *opt = 0 as *mut CINTOpt;
// }
// #[no_mangle]
// pub unsafe extern "C" fn int1e_nuc_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     *opt = 0 as *mut CINTOpt;
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_cart_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_ovlp_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_sph_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_ovlp_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_ovlp_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_cart_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_nuc_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_sph_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_nuc_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_optimizer(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: i32,
//     mut bas: *mut i32,
//     mut nbas: i32,
//     mut env: *mut f64,
// ) {
//     int1e_nuc_optimizer(opt, atm, natm, bas, nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_sph_(
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
// return int1e_ovlp_sph(
//     out,
//     0 as *mut i32,
//     shls,
//     atm,
//     *natm,
//     bas,
//     *nbas,
//     env,
//     *opt,
//     0 as *mut f64,
// );
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_sph_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int1e_ovlp_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_cart_(
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
//     unimplemented!("why do we have the fn_name_ versions of fn_name functions?");
//     //return int1e_ovlp_cart(
//     //    out,
//     //    0 as *mut i32,
//     //    shls,
//     //    atm,
//     //    *natm,
//     //    bas,
//     //    *nbas,
//     //    env,
//     //    *opt,
//     //    0 as *mut f64,
//     //);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_cart_optimizer_(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     int1e_ovlp_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_ovlp_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int1e_ovlp_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_sph_(
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
// return int1e_nuc_sph(
//     out,
//     0 as *mut i32,
//     shls,
//     atm,
//     *natm,
//     bas,
//     *nbas,
//     env,
//     *opt,
//     0 as *mut f64,
// );
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_sph_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int1e_nuc_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_cart_(
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
// return int1e_nuc_cart(
//     out,
//     0 as *mut i32,
//     shls,
//     atm,
//     *natm,
//     bas,
//     *nbas,
//     env,
//     *opt,
//     0 as *mut f64,
// );
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_cart_optimizer_(
//     mut opt: *mut *mut CINTOpt,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     int1e_nuc_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
// #[no_mangle]
// pub unsafe extern "C" fn cint1e_nuc_optimizer_(
//     mut optptr_as_integer8: size_t,
//     mut atm: *mut i32,
//     mut natm: *mut i32,
//     mut bas: *mut i32,
//     mut nbas: *mut i32,
//     mut env: *mut f64,
// ) {
//     let mut opt: *mut *mut CINTOpt = optptr_as_integer8 as *mut *mut CINTOpt;
//     int1e_nuc_optimizer(opt, atm, *natm, bas, *nbas, env);
// }
