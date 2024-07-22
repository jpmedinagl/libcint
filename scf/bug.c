#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

// #define ATM_SLOTS 6
// #define BAS_SLOTS 8

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

int __enzyme_autodiff(void *, ...);

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int main(int argc, char ** argv)
{	
	int natm = 3;
	int nbas = 5;

	int atm[natm * ATM_SLOTS] = {8, 20, 1, 23, 0, 0,
                                 1, 24, 1, 27, 0, 0, 
                                 1, 28, 1, 31, 0, 0};
	int bas[nbas * BAS_SLOTS] = {0, 0, 3, 1, 0, 32, 35, 0, 
                                 0, 0, 3, 1, 0, 38, 41, 0, 
                                 0, 1, 3, 1, 0, 44, 47, 0, 
                                 1, 0, 3, 1, 0, 50, 53, 0, 
                                 2, 0, 3, 1, 0, 50, 53, 0};
	double * env = malloc(sizeof(double) * 1000);
	double * denv;
	
	int di, dj;
	int shls[4];
	
	double * buf;
	double * dbuf;

	// SINGLE MATRIX ELEMENT

    int i = 4, j = 2;
    shls[0] = i; di = CINTcgto_cart(i, bas);
	shls[1] = j; dj = CINTcgto_cart(j, bas);

    printf("ij didj %d%d %d%d\n", i, j, di, dj);

    buf = malloc(sizeof(double) * di * dj);

    dbuf = malloc(sizeof(double) * di * dj);
	dbuf[0] = 1.0;
    denv = malloc(sizeof(double) * 1000);
	memset(denv, 0, sizeof(double) * 1000);
	
	int ret = __enzyme_autodiff((void *) cint1e_ovlp_cart,
			enzyme_dup, buf, dbuf,
			enzyme_const, shls,
			enzyme_const, atm,
			enzyme_const, natm,
			enzyme_const, bas,
			enzyme_const, nbas,
			enzyme_dup, env, denv,
			enzyme_const, NULL);

    printf("denv:\n");
    for (int k = 0; k < 56; k++) {
        printf("%f ", denv[k]);
    }
    printf("\n");

    return 0;
}