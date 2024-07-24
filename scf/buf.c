#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

int __enzyme_autodiff(void *, ...);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

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

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

void cint1e_cart_wrapper(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env) {
	cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
}

void cint1e_diff(double *buf, double * dbuf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env, double *denv) {
	__enzyme_autodiff((void *) cint1e_ovlp_cart,
			enzyme_dup, buf, dbuf,	
			enzyme_const, shls, 
			enzyme_const, atm, 
			enzyme_const, natm, 
			enzyme_const, bas, 
			enzyme_const, nbas,
			enzyme_dup, env, denv,
			enzyme_const, NULL);
}


int main(int argc, char ** argv)
{	
	int natm = 3;
	int nbas = 5;

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);
	double * denv;

	FILE * file = fopen("/u/jpmedina/libcint/molecules/h2o/sto3g.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);
	
	int di, dj;
	int shls[4];
	
	double * buf;
	double * dbuf;

	// SINGLE MATRIX ELEMENT

    int i = 2, j = 4;
    shls[0] = i; di = CINTcgto_cart(i, bas);
	shls[1] = j; dj = CINTcgto_cart(j, bas);

    printf("ij didj %d%d %d%d\n", i, j, di, dj);

    buf = malloc(sizeof(double) * di * dj);
    dbuf = malloc(sizeof(double) * di * dj);
	dbuf[0] = 1.0;
    denv = malloc(sizeof(double) * 10000);
	memset(denv, 0, sizeof(double) * 10000);
	
	printf("denv:\n");
    for (int k = 55; k < 56; k++) {
        printf("%f ", denv[k]);
    }
    cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);

	dbuf[1] = 1.0;
	memset(denv, 0, sizeof(double) * 10000);

	cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);
	for (int k = 55; k < 56; k++) {
        printf("%f ", denv[k]);
    }

	dbuf[2] = 1.0;
	memset(denv, 0, sizeof(double) * 10000);

	cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);
	for (int k = 55; k < 56; k++) {
        printf("%f ", denv[k]);
    }
    printf("\n");

	printf("finite d:\n");
    
    double grad;
    double * b1 = malloc(sizeof(double) * di * dj);
    double * b2 = malloc(sizeof(double) * di * dj);
    double h = 0.000001; 
    
    for (int k = 55; k < 56; k++) {
        env[k] += h;
		cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
		memcpy(b1, buf, sizeof(double) * di * dj);
        env[k] -= 2.0*h;
		cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
		memcpy(b2, buf, sizeof(double) * di * dj);
        env[k] += h;

		for (int l = 0; l < di * dj; l++) {
			grad = (b1[l] - b2[l])/(2.0*h);
        	printf("%lf ", grad);
		}
		printf("\n");
    }
    printf("\n");
	

	// FULL MATRIX

	// for (int i = 0; i < nbas; i++) {
	// 	for (int j = 0; j < nbas; j++) {
	// 		shls[0] = i; di = CINTcgto_cart(i, bas);
	// 		shls[1] = j; dj = CINTcgto_cart(j, bas);
    //         printf("ij didj %d%d %d%d\n", i, j, di, dj);

	// 		buf = malloc(sizeof(double) * di * dj);
	// 		dbuf = malloc(sizeof(double) * di * dj);
	// 		dbuf[0] = 1.0;

	// 		denv = malloc(sizeof(double) * 10000);
	// 		memset(denv, 0, sizeof(double) * 10000);

	// 		cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);

    //         printf("buf: ");
    //         for (int l = 0; l < di * dj; l++) {
    //             printf("%lf ", buf[l]);
    //         }
    //         printf("\n");

    //         printf("denv:");
	// 		for (int k = 24; k < 42; k++) {
	// 			printf("%f ", denv[k]);
	// 		}
	// 		printf("\n");

	// 		free(buf);
	// 		free(denv);
	// 	}
	// }
	// printf("\n");

	free(atm);
	free(bas);
	free(env);
}