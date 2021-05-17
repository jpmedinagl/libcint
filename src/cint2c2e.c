/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * 2-center 2-electron integrals
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint_bas.h"
#include "g2e.h"
#include "optimizer.h"
#include "cint1e.h"
#include "cint2e.h"
#include "misc.h"
#include "cart2sph.h"
#include "c2f.h"

#define PRIM2CTR(ctrsymb, gp, ngp) \
        if (ctrsymb##_ctr > 1) {\
                if (*ctrsymb##empty) { \
                        CINTprim_to_ctr_0(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } else { \
                        CINTprim_to_ctr_1(gctr##ctrsymb, gp, c##ctrsymb+ctrsymb##p, \
                                          ngp, ctrsymb##_prim, ctrsymb##_ctr, \
                                          non0ctr##ctrsymb[ctrsymb##p], \
                                          non0idx##ctrsymb+ctrsymb##p*ctrsymb##_ctr); \
                } \
        } \
        *ctrsymb##empty = 0


FINT CINT2c2e_loop_nopt(double *gctr, CINTEnvVars *envs, double *cache)
{
        FINT *shls  = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT k_sh = shls[1];
        FINT i_ctr = envs->x_ctr[0];
        FINT k_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        FINT n_comp = envs->ncomp_tensor;
        double fac1i, fac1k;
        FINT ip, kp;
        FINT empty[3] = {1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *kempty = empty + 1;
        FINT *gempty = empty + 2;
        /* COMMON_ENVS_AND_DECLARE end */
        const FINT nc = i_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        FINT *idx;
        MALLOC_INSTACK(idx, envs->nf * 3);
        CINTg1e_index_xyz(idx, envs);

        FINT *non0ctri, *non0ctrk;
        FINT *non0idxi, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+k_prim+i_prim*i_ctr+k_prim*k_ctr);
        non0ctrk = non0ctri + i_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxk = non0idxi + i_prim*i_ctr;
        if (i_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        }
        if (k_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        }

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp]; // to use CINTg0_2e
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                PRIM2CTR(i, gout, envs->nf*n_comp);
                        }
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

FINT CINT2c2e_loop(double *gctr, CINTEnvVars *envs, const CINTOpt *opt, double *cache)
{
        FINT *shls = envs->shls;
        FINT *bas = envs->bas;
        double *env = envs->env;
        FINT i_sh = shls[0];
        FINT k_sh = shls[1];
        FINT i_ctr = envs->x_ctr[0];
        FINT k_ctr = envs->x_ctr[1];
        FINT i_prim = bas(NPRIM_OF, i_sh);
        FINT k_prim = bas(NPRIM_OF, k_sh);
        double *ai = env + bas(PTR_EXP, i_sh);
        double *ak = env + bas(PTR_EXP, k_sh);
        double *ci = env + bas(PTR_COEFF, i_sh);
        double *ck = env + bas(PTR_COEFF, k_sh);
        FINT n_comp = envs->ncomp_tensor;
        double fac1i, fac1k;
        FINT ip, kp;
        FINT empty[3] = {1, 1, 1};
        FINT *iempty = empty + 0;
        FINT *kempty = empty + 1;
        FINT *gempty = empty + 2;
        FINT *non0ctri, *non0ctrk;
        FINT *non0idxi, *non0idxk;
        MALLOC_INSTACK(non0ctri, i_prim+k_prim+i_prim*i_ctr+k_prim*k_ctr);
        non0ctrk = non0ctri + i_prim;
        non0idxi = non0ctrk + k_prim;
        non0idxk = non0idxi + i_prim*i_ctr;
        if (i_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxi, non0ctri, ci, i_prim, i_ctr);
        }
        if (k_ctr > 1) {
                CINTOpt_non0coeff_byshell(non0idxk, non0ctrk, ck, k_prim, k_ctr);
        }

        FINT *idx = opt->index_xyz_array[envs->i_l*LMAX1+envs->k_l];

        const FINT nc = i_ctr * k_ctr;
        const FINT leng = envs->g_size * 3 * ((1<<envs->gbits)+1);
        const FINT lenk = envs->nf * nc * n_comp; // gctrk
        const FINT leni = envs->nf * i_ctr * n_comp; // gctri
        const FINT len0 = envs->nf * n_comp; // gout
        const FINT len = leng + lenk + leni + len0;
        double *g;
        MALLOC_INSTACK(g, len);
        double *g1 = g + leng;
        double *gout, *gctri, *gctrk;

        if (n_comp == 1) {
                gctrk = gctr;
        } else {
                gctrk = g1;
                g1 += lenk;
        }
        if (k_ctr == 1) {
                gctri = gctrk;
                iempty = kempty;
        } else {
                gctri = g1;
                g1 += leni;
        }
        if (i_ctr == 1) {
                gout = gctri;
                gempty = iempty;
        } else {
                gout = g1;
        }

        *kempty = 1;
        for (kp = 0; kp < k_prim; kp++) {
                envs->ak = ak[kp];
                envs->akl = ak[kp];
                if (k_ctr == 1) {
                        fac1k = envs->common_factor * ck[kp];
                } else {
                        fac1k = envs->common_factor;
                        *iempty = 1;
                }
                for (ip = 0; ip < i_prim; ip++) {
                        envs->ai = ai[ip];
                        envs->aij = ai[ip];
                        if (i_ctr == 1) {
                                fac1i = fac1k*ci[ip];
                        } else {
                                fac1i = fac1k;
                        }
                        if ((*envs->f_g0_2e)(g, fac1i, envs)) {
                                (*envs->f_gout)(gout, g, idx, envs, *gempty);
                                PRIM2CTR(i, gout, envs->nf*n_comp);
                        }
                } // end loop i_prim
                if (!*iempty) {
                        PRIM2CTR(k, gctri, envs->nf*i_ctr*n_comp);
                }
        } // end loop k_prim

        if (n_comp > 1 && !*kempty) {
                CINTdmat_transpose(gctr, gctrk, envs->nf*nc, n_comp);
        }
        return !*kempty;
}

CACHE_SIZE_T CINT2c2e_drv(double *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                          double *cache, void (*f_c2s)())
{
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n;
        FINT has_value;

        if (opt != NULL) {
                has_value = CINT2c2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs, cache);
        }

        FINT counts[4];
        if (f_c2s == &c2s_sph_1e) {
                counts[0] = (envs->i_l*2+1) * x_ctr[0];
                counts[1] = (envs->k_l*2+1) * x_ctr[1];
        } else {
                counts[0] = envs->nfi * x_ctr[0];
                counts[1] = envs->nfk * x_ctr[1];
        }
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_c2s)(out+nout*n, gctr+nc*n, dims, envs, cache);
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_dset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}
// (spinor|spinor)
CACHE_SIZE_T CINT2c2e_spinor_drv(double complex *out, FINT *dims, CINTEnvVars *envs, CINTOpt *opt,
                        double *cache, void (*f_e1_c2s)())
{
        if (envs->ncomp_e1 > 1 || envs->ncomp_e2 > 1) {
                fprintf(stderr, "CINT2c2e_spinor_drv not implemented\n");
                exit(1);
        }
        if (out == NULL) {
                return int1e_cache_size(envs);
        }
        FINT *x_ctr = envs->x_ctr;
        FINT nc = envs->nf * x_ctr[0] * x_ctr[1];
        FINT n_comp = envs->ncomp_e1 * envs->ncomp_e2 * envs->ncomp_tensor;
        double *stack = NULL;
        if (cache == NULL) {
                size_t cache_size = int1e_cache_size(envs);
                stack = malloc(sizeof(double)*cache_size);
                cache = stack;
        }
        double *gctr;
        MALLOC_INSTACK(gctr, nc*n_comp);

        FINT n, has_value;

        if (opt != NULL) {
                has_value = CINT2c2e_loop(gctr, envs, opt, cache);
        } else {
                has_value = CINT2c2e_loop_nopt(gctr, envs, cache);
        }

        FINT counts[4];
        counts[0] = CINTcgto_spinor(envs->shls[0], envs->bas);
        counts[1] = CINTcgto_spinor(envs->shls[1], envs->bas);
        counts[2] = 1;
        counts[3] = 1;
        if (dims == NULL) {
                dims = counts;
        }
        FINT nout = dims[0] * dims[1];
        if (has_value) {
                for (n = 0; n < n_comp; n++) {
                        (*f_e1_c2s)(out+nout*n, gctr, dims, envs, cache);
                        gctr += nc;
                }
        } else {
                for (n = 0; n < n_comp; n++) {
                        c2s_zset0(out+nout*n, dims, counts);
                }
        }
        if (stack != NULL) {
                free(stack);
        }
        return has_value;
}


CACHE_SIZE_T int2c2e_sph(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_sph_1e);
}
void int2c2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                       FINT *bas, FINT nbas, double *env)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTall_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env);
}

CACHE_SIZE_T int2c2e_cart(double *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                 FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_drv(out, dims, &envs, opt, cache, &c2s_cart_1e);
}
 
CACHE_SIZE_T int2c2e_spinor(double complex *out, FINT *dims, FINT *shls, FINT *atm, FINT natm,
                   FINT *bas, FINT nbas, double *env, CINTOpt *opt, double *cache)
{
        FINT ng[] = {0, 0, 0, 0, 0, 1, 1, 1};
        CINTEnvVars envs;
        CINTinit_int2c2e_EnvVars(&envs, ng, shls, atm, natm, bas, nbas, env);
        envs.f_gout = &CINTgout2e;
        return CINT2c2e_spinor_drv(out, dims, &envs, opt, cache, &c2s_sf_1e);
}


ALL_CINT(int2c2e)
ALL_CINT_FORTRAN_(int2c2e)

