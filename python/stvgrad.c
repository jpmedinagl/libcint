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

void cint1e_cart_wrapper(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env) {
	cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
}

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

void cint1e_diff(double *buf, double * dbuf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env, double *denv) {
    printf("before enzyme call\n");
	__enzyme_autodiff((void *) cint1e_cart_wrapper,
			enzyme_dup, buf, dbuf,	
			enzyme_const, shls, 
			enzyme_const, atm, 
			enzyme_const, natm, 
			enzyme_const, bas, 
			enzyme_const, nbas,
			enzyme_dup, env, denv);
    printf("after enzyme call\n");
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
	
	printf("denv:\n");

	for (int i = 0; i < nbas; i++) {
		for (int j = 0; j < nbas; j++) {
			shls[0] = i; di = CINTcgto_cart(i, bas);
			shls[1] = j; dj = CINTcgto_cart(j, bas);

            printf("i  j  %d %d \ndi dj %d %d\n", i, j, di, dj);

			buf = malloc(sizeof(double) * di * dj);
			dbuf = malloc(sizeof(double) * di * dj);
			dbuf[0] = 1.0;

			denv = malloc(sizeof(double) * 10000);
			memset(denv, 0, sizeof(double) * 10000);

			cint1e_diff(buf, dbuf, shls, atm, natm, bas, nbas, env, denv);
            printf("\n");
			
			printf("denv[50:56]: ");
			for (int k = 28; k < 34; k++) {
				printf("%f ", denv[k]);
			}
			printf("\n");

			free(buf);
			free(denv);
		}
	}
	printf("\n");

	free(atm);
	free(bas);
	free(env);
}
