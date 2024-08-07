#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cint.h"

int cint1e_ipnuc_cart(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env);

/* general contracted DZ basis [3s1p/2s1p] for H2
    exponents    contract-coeff
S   6.0          0.7               0.4
    2.0          0.6               0.3
    0.8          0.5               0.2
P   0.9          1.
 */

int main()
{
        printf("TRYING THINGS OUT\n");
	int natm = 2;
        int nbas = 4;
        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        int *atm = malloc(sizeof(int) * natm * ATM_SLOTS);
        int *bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
        double *env = malloc(sizeof(double) * 10000);

        int i, n, off;
        off = PTR_ENV_START; // = 20

        i = 0;
        atm[CHARGE_OF + ATM_SLOTS * i] = 1;
        atm[PTR_COORD + ATM_SLOTS * i] = off;
        env[off + 0] =  0; // x (Bohr)
        env[off + 1] =  0; // y (Bohr)
        env[off + 2] =-.8; // z (Bohr)
        i++;
        off += 3;
        atm[CHARGE_OF + ATM_SLOTS * i] = 1;
        atm[PTR_COORD + ATM_SLOTS * i] = off;
        env[off + 0] = 0;
        env[off + 1] = 0;
        env[off + 2] =.8; // (Bohr)
        i++;
        off += 3;

        n = 0;
        /* basis #0, 3s -> 2s */
        bas[ATOM_OF  + BAS_SLOTS * n]  = 0;
        bas[ANG_OF   + BAS_SLOTS * n]  = 0;
        bas[NPRIM_OF + BAS_SLOTS * n]  = 3;
        bas[NCTR_OF  + BAS_SLOTS * n]  = 2;
        bas[PTR_EXP  + BAS_SLOTS * n]  = off;
        env[off + 0] = 6.;
        env[off + 1] = 2.;
        env[off + 2] = .8;
        off += 3;
        bas[PTR_COEFF+ BAS_SLOTS * n] = off;
        env[off + 0] = .7 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+0]);
        env[off + 1] = .6 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+1]);
        env[off + 2] = .5 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+2]);
        env[off + 3] = .4 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+0]);
        env[off + 4] = .3 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+1]);
        env[off + 5] = .2 * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]+2]);
        off += 6;
        n++;

        /* basis #1 */
        bas[ATOM_OF  + BAS_SLOTS * n]  = 0;
        bas[ANG_OF   + BAS_SLOTS * n]  = 1;
        bas[NPRIM_OF + BAS_SLOTS * n]  = 1;
        bas[NCTR_OF  + BAS_SLOTS * n]  = 1;
        bas[PTR_EXP  + BAS_SLOTS * n]  = off;
        env[off + 0] = .9;
        off += 1;
        bas[PTR_COEFF+ BAS_SLOTS * n] = off;
        env[off + 0] = 1. * CINTgto_norm(bas[ANG_OF+BAS_SLOTS*n], env[bas[PTR_EXP+BAS_SLOTS*n]]);
        off += 1;
        n++;

        /* basis #2 == basis #0 */
        bas[ATOM_OF  + BAS_SLOTS * n] = 1;
        bas[ANG_OF   + BAS_SLOTS * n] = bas[ANG_OF   + BAS_SLOTS * 0];
        bas[NPRIM_OF + BAS_SLOTS * n] = bas[NPRIM_OF + BAS_SLOTS * 0];
        bas[NCTR_OF  + BAS_SLOTS * n] = bas[NCTR_OF  + BAS_SLOTS * 0];
        bas[PTR_EXP  + BAS_SLOTS * n] = bas[PTR_EXP  + BAS_SLOTS * 0];
        bas[PTR_COEFF+ BAS_SLOTS * n] = bas[PTR_COEFF+ BAS_SLOTS * 0];
        n++;

        /* basis #3 == basis #1 */
        bas[ATOM_OF  + BAS_SLOTS * n] = 1;
        bas[ANG_OF   + BAS_SLOTS * n] = bas[ANG_OF   + BAS_SLOTS * 1];
        bas[NPRIM_OF + BAS_SLOTS * n] = bas[NPRIM_OF + BAS_SLOTS * 1];
        bas[NCTR_OF  + BAS_SLOTS * n] = bas[NCTR_OF  + BAS_SLOTS * 1];
        bas[PTR_EXP  + BAS_SLOTS * n] = bas[PTR_EXP  + BAS_SLOTS * 1];
        bas[PTR_COEFF+ BAS_SLOTS * n] = bas[PTR_COEFF+ BAS_SLOTS * 1];
        n++;

        /*
         * call one-electron cartesian integrals
         * the integral has 3 components, saving as
         * buf[      0:  di*dj]    for x
         * buf[  di*dj:2*di*dj]    for y
         * buf[2*di*dj:3*di*dj]    for z
         */
        int j, k, l;
        int di, dj, dk, dl;
        int shls[4];
        double *buf;

        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas);
        buf = malloc(sizeof(double) * di * dj * 3);
        if (0 != cint1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, env)) {
                printf("This gradient integral is not 0.\n");
        } else {
                printf("This integral is 0.\n");
        }
        free(buf);

        /*
         * call two-electron cartesian integrals
         */
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas);
        buf = malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL)) {
                printf("This integral is not 0.\n");
        } else {
                printf("This integral is 0.\n");
        }
        free(buf);

        CINTOpt *opt = NULL;
        cint2e_cart_optimizer(&opt, atm, natm, bas, nbas, env);
        i = 0; shls[0] = i; di = CINTcgto_cart(i, bas);
        j = 1; shls[1] = j; dj = CINTcgto_cart(j, bas);
        k = 2; shls[2] = k; dk = CINTcgto_cart(k, bas);
        l = 2; shls[3] = l; dl = CINTcgto_cart(l, bas);
        buf = malloc(sizeof(double) * di * dj * dk * dl);
        if (0 != cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt)) {
                printf("This integral is not 0.\n");
        } else {
                printf("This integral is 0.\n");
        }
        free(buf);
        CINTdel_optimizer(&opt);

        free(atm);
        free(bas);
        free(env);
}
