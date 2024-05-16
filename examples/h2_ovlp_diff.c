#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

int cint1e_ovlp_cart(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env);

int __enzyme_autodiff(void *, ...);


void read_arrays(int natm, int nbas, int ** atm, int ** bas, double ** env)
{
	FILE * fenv = fopen("/u/jpmedina/molecules/h2/env.txt", "r");
	FILE * fatm = fopen("/u/jpmedina/molecules/h2/atm.txt", "r");
	FILE * fbas = fopen("/u/jpmedina/molecules/h2/bas.txt", "r");

	double num;
	int number;
	int count = 0;

	while (fscanf(fenv, "%lf", &num) != EOF) {
		(*env)[count] = num;
		count++;
	}
	fclose(fenv);

	count = 0;
	while (fscanf(fatm, "%d", &number) != EOF) {
		(*atm)[count] = number;
		count++;
	}
	fclose(fatm);

	count = 0;
	while (fscanf(fbas, "%d", &number) != EOF) {
		(*bas)[count] = number;
		count++;
	}
	fclose(fbas);
}

void cint1e_cart_wrapper(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env) {
	cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
}

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

void cint1e_diff(double *buf, double * dbuf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env, double *denv) {
	__enzyme_autodiff((void *) cint1e_cart_wrapper,
			enzyme_dup, buf, dbuf,	
			enzyme_const, shls, 
			enzyme_const, atm, 
			enzyme_const, natm, 
			enzyme_const, bas, 
			enzyme_const, nbas,
			enzyme_dup, env, denv);
}

int main()
{
	int natm = 2;
	int nbas = 2;

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);
	double * denv; // = malloc(sizeof(double) * 10000);

 	read_arrays(natm, nbas, &atm, &bas, &env);

	// int index = bas[0*8 + 5];
	// 0 for first exponent offset
	
	// double * env_diff = &(env[index]);
	// env_diff points to the differentiable part of env (everything after the first exponent)

	int di, dj;
	int shls[4];
	
	double * buf;
	double * dbuf;
	
	printf("denv:\n");
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			shls[0] = i; di = CINTcgto_cart(i, bas);
			shls[1] = j; dj = CINTcgto_cart(j, bas);

			buf = malloc(sizeof(double) * di * dj);
			dbuf = malloc(sizeof(double) * di * dj);
			dbuf[0] = 1.0;

			denv = malloc(sizeof(double) * 10000);
			memset(denv, 0, sizeof(double) * 10000);

			cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);
			
			// printf("%f\n", buf[0]);

			for (int k = 28; k < 34; k++) {
				printf("%f ", denv[k]);
			}
			printf("\n");

			free(buf);
			free(denv);
		}
		printf("\n");
	}
	printf("\n");

	free(atm);
	free(bas);
	free(env);
}
