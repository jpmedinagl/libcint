#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

int __enzyme_autodiff(void *, ...);

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

void read_arrays(FILE * file, int natm, int nbas, int ** atm, int ** bas, double ** env)
{
	double num;
	int number;
	int count = 0;

	while (count < natm * ATM_SLOTS) {
        fscanf(file, "%d", &number);
		(*atm)[count] = number;
		count++;
	}

	count = 0;
	while (count < nbas * BAS_SLOTS) {
        fscanf(file, "%d", &number);
		(*bas)[count] = number;
		count++;
	}

	count = 0;
	while (fscanf(file, "%lf", &num) != EOF) {
		(*env)[count] = num;
		count++;
	}
	fclose(file);
}

double * c_arr(int x, int y) 
{   
    double * mat = malloc(sizeof(double) * (x * y));
    memset(mat, 0, sizeof(double) * (x * y));
    return mat;
}

void print_arr(int n, int size, double * A) {
    for (int i = 0; i < pow(n, size); i++) {
        printf("%lf ", A[i]);
        for (int p = 1; p < size; p++) {
            if ((i + 1) % (int) pow(n, p) == 0) {
                printf("\n");
            }
        }
    }
}

void int1e(int natm, int nbas, int nshells, int * atm, int * bas, double * env, double ** R) 
{   
    *R = c_arr(nshells, nshells);

    int mu, nu;
    int di, dj;
    int shls[4];

    int c;
    
    double * buf;

    mu = 0;
    for (int i = 0; i < nbas; i++) {
        nu = 0;
        for (int j = 0; j < nbas; j++) {
            shls[0] = i; di = CINTcgto_cart(i, bas);
            shls[1] = j; dj = CINTcgto_cart(j, bas);

            buf = malloc(sizeof(double) * di * dj);
            cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);

            c = 0;
            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    (*R)[mui*nshells + nuj] = buf[c];
                    c++;
                }
            }


            nu += dj;
        }
        mu += di;
    }
}

// void int2e(int natm, int nbas, int nshells, int * atm, int * bas, double * env, double ** R) 
// {   
//     *R = c_arr(nshells, nshells);

//     int mu, nu, sig, lam;
//     int di, dj, dk, dl;
//     int shls[4];

//     int c;
    
//     double * buf;

//     mu = 0;
//     for (int i = 0; i < nbas; i++) {
//         nu = 0;
//         for (int j = 0; j < nbas; j++) {
//             sig = 0;

//             shls[0] = i; di = CINTcgto_cart(i, bas);
//             shls[1] = j; dj = CINTcgto_cart(j, bas);

//             buf = malloc(sizeof(double) * di * dj);
//             cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);

//             c = 0;
//             for (int nuj = nu; nuj < nu + dj; nuj++) {
//                 for (int mui = mu; mui < mu + di; mui++) {
//                     (*S)[mui * nshells + nuj] = buf[c];
//                     c++;
//                 }
//             }

//             c = 0;
//             cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
//             for (int nuj = nu; nuj < nu + dj; nuj++) {
//                 for (int mui = mu; mui < mu + di; mui++) {
//                     T[mui * nshells + nuj] = buf[c];
//                     c++;
//                 }
//             }

//             c = 0;
//             cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
//             for (int nuj = nu; nuj < nu + dj; nuj++) {
//                 for (int mui = mu; mui < mu + di; mui++) {
//                     V[mui * nshells + nuj] = buf[c];
//                     c++;
//                 }
//             }

//             for (int nuj = nu; nuj < nu + dj; nuj++) {
//                 for (int mui = mu; mui < mu + di; mui++) {
//                     (*H)[mui*nshells + nuj] = T[mui*nshells + nuj] + V[mui*nshells + nuj];
//                 }
//             }

//             free(buf);

//             for (int k = 0; k < nbas; k++) {
//                 lam = 0;
//                 for (int l = 0; l < nbas; l++) {
//                     shls[2] = k; dk = CINTcgto_cart(k, bas);
//                     shls[3] = l; dl = CINTcgto_cart(l, bas);

//                     buf = malloc(sizeof(double) * di * dj * dk * dl);
//                     cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);

//                     c = 0;
//                     for (int laml = lam; laml < lam + dl; laml++) {
//                         for (int sigk = sig; sigk < sig + dk; sigk++) {
//                             for (int nuj = nu; nuj < nu + dj; nuj++) {
//                                 for (int mui = mu; mui < mu + di; mui++) {
//                                     if (i == 4 && j == 4 && k == 4 & l == 4) {
//                                         printf("before last\n");
//                                     }
//                                     (*two)[(int) (mui*pow(nshells, 3) + nuj*pow(nshells, 2) + sigk*nshells + laml)] = buf[c];
//                                     c++;
//                                     if (i == 4 && j == 4 && k == 4 & l == 4) {
//                                         printf("after last\n");
//                                     }
//                                 }
//                             }
//                         }
//                     }

//                     free(buf);

//                     lam += dl;
//                 }
//                 sig += dk;
//             }
//             if (i == 4 && j == 4) {
//                 printf("after last\n");
//             }
//             nu += dj;
//         }
//         mu += di;
//     }
//     free(T);
//     free(V);
// }

void cint1e_diff(double ** R, double ** dR, int natm, int nbas, int nshells, int * atm, int * bas, double * env, double * denv)  {
	__enzyme_autodiff((void *) int1e,
            enzyme_dup, R, dR,
            enzyme_const, natm, 
            enzyme_const, nbas, 
            enzyme_const, nshells, 
            enzyme_const, atm,
            enzyme_const, bas, 
            enzyme_dup, env, denv);
}

int main(int argc, char ** argv)
{	
    int natm = 2;
    int nbas = 2;
    int nshells = 2;

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);
	double * denv = malloc(sizeof(double) * 10000);

	FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/sto3g.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);
    
    double * R;
    int1e(natm, nbas, nshells, atm, bas, env, &R);
    printf("I:\n");
    print_arr(nshells, 2, R);

    double * dR = c_arr(nshells, nshells);
    dR[0] = 1.0;
	
    cint1e_diff(&R, &dR, natm, nbas, nshells, atm, bas, env, denv);
    printf("I:\n");
    print_arr(nshells, 2, R);

    printf("denv:\n");
    for (int k = 20; k < 28; k++) {
        printf("%f ", denv[k]);
    }
    printf("\n");

	free(atm);
	free(bas);
	free(env);
}
