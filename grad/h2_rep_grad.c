#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

int __enzyme_autodiff(void *, ...);

void cint2e_cart_wrapper(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env) {
	cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);
}

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

void cint2e_diff(double *buf, double * dbuf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env, double *denv) {
	__enzyme_autodiff((void *) cint2e_cart_wrapper,
			enzyme_dup, buf, dbuf,	
			enzyme_const, shls, 
			enzyme_const, atm, 
			enzyme_const, natm, 
			enzyme_const, bas, 
			enzyme_const, nbas,
			enzyme_dup, env, denv);
}


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

int main(int argc, char ** argv)
{	
	int natm = 2;
	int nbas = 2;

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);
	double * denv;

	FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/sto3g.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);
	
	int di, dj, dk, dl;
	int shls[4];
	
	double * buf;
	double * dbuf;
	
	printf("denv:\n");
   	for (int i = 0; i < nbas; i++) {
		for (int j = 0; j < nbas; j++) {
			for (int k = 0; k < nbas; k++) {
				for (int l = 0; l < nbas; l++) {
					// printf("%d %d %d %d\n", i, j, k, l);
					
					shls[0] = i; di = CINTcgto_cart(i, bas);
					shls[1] = j; dj = CINTcgto_cart(j, bas);
					shls[2] = k; dk = CINTcgto_cart(k, bas);
					shls[3] = l; dl = CINTcgto_cart(l, bas);

					buf = malloc(sizeof(double) * di * dj * dk * dl);
                    			dbuf = malloc(sizeof(double) * di * dj * dk * dl);
                    			dbuf[0] = 1.0;

        			        denv = malloc(sizeof(double) * 10000);
			        	memset(denv, 0, sizeof(double) * 10000);

					cint2e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);

					for (int m = 28; m < 34; m++) {
                        			printf("%f ", denv[m]);
                    			}
					printf("\n");

					free(buf);
			        	free(denv);
				}
			}
		}
	}
	printf("\n");

	free(atm);
	free(bas);
	free(env);
}
