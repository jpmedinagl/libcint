#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cint.h"

int cint1e_ovlp_cart(double *buf, int *shls,
                      int *atm, int natm, int *bas, int nbas, double *env);

void read_arrays(int natm, int nbas, int ** atm, int ** bas, double ** env)
{
    	FILE * fenv = fopen("/u/jpmedina/molecules/h2o/env.txt", "r");
    	FILE * fatm = fopen("/u/jpmedina/molecules/h2o/atm.txt", "r");
    	FILE * fbas = fopen("/u/jpmedina/molecules/h2o/bas.txt", "r");

    	double num;
    	int number;
    	int count = 0;

    	while (fscanf(fenv, "%lf", &num) != EOF) {
        	(*env)[count] = num;
        	count++;
    	}
    	fclose(fenv);

    	printf("env: (%d)\n", count);
    	for (int i = 0; i < count; i++) {
        	printf("%lf ", (*env)[i]);
    	}	
    	printf("\n");

    	count = 0;
    	while (fscanf(fatm, "%d", &number) != EOF) {
        	(*atm)[count] = number;
        	count++;
    	}
    	fclose(fatm);

    	printf("atm: (%d)\n", count);
    	for (int i = 0; i < count; i++) {
        	printf("%d ", (*atm)[i]);
    	}
    	printf("\n");

    	count = 0;
    	while (fscanf(fbas, "%d", &number) != EOF) {
        	(*bas)[count] = number;
        	count++;
    	}
    	fclose(fbas);

    	printf("bas: (%d)\n", count);
    	for (int i = 0; i < count; i++) {
        	printf("%d ", (*bas)[i]);
    	}
    	printf("\n");
}

int main()
{
	int natm = 3;
    	int nbas = 5;

    	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
    	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
    	double * env = malloc(sizeof(double) * 10000);

 	read_arrays(natm, nbas, &atm, &bas, &env);

        /*
         * call one-electron cartesian integrals
         * the integral has 3 components, saving as
         * buf[      0:  di*dj]    for x
         * buf[  di*dj:2*di*dj]    for y
         * buf[2*di*dj:3*di*dj]    for z
         */
        
        int di, dj;
        int shls[4];
        // double *buf;
	
	printf("buf:\n");
	for (int i = 0; i < nbas; i++) {
		for (int j = 0; j < nbas; j++) {
			double * buf;

         		shls[0] = i; di = CINTcgto_cart(i, bas);
         		shls[1] = j; dj = CINTcgto_cart(j, bas);

        		// printf("di: %d\n", di);
        		// printf("di: %d dj: %d\n", di, dj);

        		buf = malloc(sizeof(double) * di * dj);
        		cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);

			for (int k = 0; k < di * dj; k++) {
				printf("%f ", buf[k]);
			}
			free(buf);
		}
		printf("\n");
	}
	printf("\n");

        free(atm);
        free(bas);
        free(env);
}
